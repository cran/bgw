# Internal functions that assemble what bgw_mle returns. None is exported, and no help page is
# generated.

#' @title bgw_buildModel
#'
#' @description Builds the model object returned by bgw_mle from the final iv and v.
#'
#' @param iv BGW internal vector of integer values, as left by the estimation loop.
#' @param v BGW internal vector of numeric values, as left by the estimation loop.
#' @param beta Final parameter vector.
#' @param betaStart Starting parameter vector.
#' @param p Number of parameters.
#' @param n Number of residuals.
#' @param bgw_settings List of final settings.
#' @param hasAnalyticGrad Logical. Whether the caller supplied calcJ.
#' @param g_index_v Positions of the gradient in v.
#' @param d_index_v Positions of the scale vector in v.
#' @param silent Logical.
#'
#' @return List: the model object, without its class.
#' @noRd
bgw_buildModel <- function(iv, v, beta, betaStart, p, n, bgw_settings, hasAnalyticGrad, g_index_v,
                           d_index_v, silent) {
  # Convergence codes and messages (see bgw_constants.R)
  bgw_message_table <- bgw_messages(annotated = TRUE)
  model <- list()
  model$bgw_settings    <- bgw_settings
  model$hasAnalyticGrad <- hasAnalyticGrad
  model$betaStart       <- betaStart
  model$numParams       <- p
  model$numResids       <- n

  # Codes above 62 report failures (63, 64 and 66 at the start; 65 when the gradient could not be
  # computed); subtracting 51 gives the position of the message in the table, as in bgw_itsum. The
  # code returned is the one from the Fortran.
  iv1 <- iv[1]
  if (iv1 > 62) iv1 <- iv1 - 51
  if (iv1 >= 3 & iv1 <= 15) {
    model$code <- iv[1]
    model$message <- trimws(bgw_message_table[[iv1]])
  } else {
    model$code <- 0
    model$message <- "Problematic BGW return code. Seek expert (maintainer) help."
  }

  # DSB note: We need to consider what to do with some of these
  # values if there is a catastrophic failure, e.g., if beta
  # has garbage.
  model$betaStop        <- beta
  model$finalLL         <- -v[bgw_v_index$f]
  model$iterations      <- iv[bgw_iv_index$niter]
  model$functionEvals   <- iv[bgw_iv_index$nfcall]
  model$gradient        <- v[g_index_v]
  model$scaleVec        <- v[d_index_v]

  # 7 November 2026: estimate and maximum are always returned, whether or not
  # the search converged (see model$code and model$message for that). The
  # previous if-else had two identical branches. Returning NULL on
  # non-convergence, as the old help text said, would break Apollo, which uses
  # model$estimate and model$maximum after unsuccessful estimations too.
  model$estimate <- beta
  model$maximum  <- model$finalLL

  vcHessianRecipCond <- v[bgw_v_index$rcond]^2
  if (iv[bgw_iv_index$covmat] > 0) {
    hessianMethodAttempted <- bgw_settings$vcHessianMethod
    hessianMethodUsed      <- bgw_settings$vcHessianMethod
    varcov <- matrix(NA, nrow = length(beta), ncol = length(beta),
                     dimnames = list(names(beta), names(beta)))
    se           <- rep(0, p)
    names(se)    <- names(beta)
    tstat        <- rep(0, p)
    names(tstat) <- names(beta)
    vc_start <- iv[bgw_iv_index$covmat]   # position in v where the covariance matrix starts
    vcVec    <- v[vc_start:(vc_start + p*(p+1)/2 - 1)]
    iii <- 0
    for (ii in 1:p) {
      for (jj in 1:ii) {
        iii <- iii + 1
        varcov[ii, jj] <- vcVec[iii]
        varcov[jj, ii] <- vcVec[iii]
      }
      se[ii] <- sqrt(varcov[ii, ii])
      tstat[ii] <- beta[ii]/se[ii]
    }
  } else {
    if (iv[bgw_iv_index$rdreq] == 0) {
      # cat('       No variance-covariance matrix computation was requested. \n')
      hessianMethodUsed <- "None"
      vcHessianRecipCond <- -1
      if (!silent) cat("\n")
    } else {
      hessianMethodAttempted <- bgw_settings$vcHessianMethod
      hessianMethodUsed      <- "N/A - BHHH Hessian is singular"
    }
    vcVec  <- NULL
    varcov <- NULL
    se     <- NULL
    tstat  <- NULL
  }

  # The reciprocal of the condition number: 0 means singular, and -1 means not computed.
  model$vcHessianRecipCondBGW    <- vcHessianRecipCond
  # Deprecated name, which is misleading because the value is the reciprocal. It is kept for now
  # because Apollo tests whether this field exists.
  model$vcHessianConditionNumber <- vcHessianRecipCond
  # 7 November 2026: bug fix. hessianMethodAttempted is documented as part of
  # the result but was never stored in model.
  model$hessianMethodAttemptedBGW <- bgw_settings$vcHessianMethod
  model$hessianMethodUsedBGW      <- hessianMethodUsed
  model$vcVecBGW                  <- vcVec
  model$varcovBGW                 <- varcov
  model$seBGW                     <- se
  model$tstatBGW                  <- tstat

  return(model)
}

#' @title bgw_indLevelGradients
#'
#' @description Computes individual-level gradients at the final parameter vector.
#'
#' @param x Final parameter vector, in x order.
#' @param beta Final parameter vector, with any names.
#' @param v BGW internal vector of numeric values, as left by the estimation loop.
#' @param d_index_v Positions of the scale vector in v.
#' @param n Number of residuals.
#' @param p Number of parameters.
#' @param hasAnalyticGrad Logical. Whether the caller supplied calcJ.
#' @param calcJ_type The calcJ_type in use (NULL if calcJ was not supplied).
#' @param evalCalcR Function: calcR evaluated at a parameter vector in x order.
#' @param evalCalcJ Function: calcJ evaluated at a parameter vector in x order.
#' @param orient_calcJ Function: checks a calcJ result and returns it as p x n.
#' @param silent Logical.
#'
#' @return n x p matrix whose row i is d log(R_i) / d beta.
#' @noRd
bgw_indLevelGradients <- function(x, beta, v, d_index_v, n, p, hasAnalyticGrad, calcJ_type,
                                  evalCalcR, evalCalcJ, orient_calcJ, silent) {
  # 7 November 2026: individual-level gradients, returned as an n x p matrix
  # whose row i is d log(R_i) / d beta, at the final beta. They are computed
  # here, after the estimation loop, using local copies only, so requesting
  # them cannot change iv, v, the estimates or the covariance matrix.
  P_final <- as.numeric(evalCalcR(x))
  if (hasAnalyticGrad) {
    J_final <- evalCalcJ(x)
    J_final <- t(orient_calcJ(J_final))            # n x p
    if (calcJ_type == "probability") {
      J_final <- J_final / P_final                 # dR_i -> dlog(R_i): row i divided by R_i
    }
  } else {
    # Forward differences of R at the final x, with the same step size
    # as the finite-difference branch in the estimation loop.
    J_final <- matrix(NA_real_, nrow = n, ncol = p)
    d_final <- v[d_index_v]
    for (k in 1:p) {
      xx <- x
      hk <- v[bgw_v_index$dltfdj]*max(abs(x[k]), 1/d_final[k])
      xx[k] <- x[k] + hk
      Pk <- as.numeric(evalCalcR(xx))
      if (any(is.na(Pk))) {
        # Try a step in the other direction before giving up on this component.
        hk <- -hk
        xx[k] <- x[k] + hk
        Pk <- as.numeric(evalCalcR(xx))
      }
      if (!any(is.na(Pk))) J_final[, k] <- (Pk - P_final)/hk
    }
    if (any(is.na(J_final)) && !silent) {
      cat("\nWarning: some individual-level gradients could not be computed (NA).\n")
    }
    J_final <- J_final / P_final                   # dR_i -> dlog(R_i)
  }
  dimnames(J_final) <- list(NULL, names(beta))
  return(J_final)
}

#' @title bgw_startModel
#'
#' @description Builds the model object for a request made without estimation: everything is
#'   evaluated at betaStart, and the Fortran is not called.
#'
#' @param x Starting parameter vector, in x order.
#' @param beta Starting parameter vector, with any names.
#' @param betaStart Starting parameter vector, as supplied.
#' @param Prob Result of calcR at betaStart.
#' @param iv BGW internal vector of integer values, after setup.
#' @param v BGW internal vector of numeric values, after setup.
#' @param n Number of residuals.
#' @param p Number of parameters.
#' @param bgw_settings List of final settings.
#' @param hasAnalyticGrad Logical. Whether the caller supplied calcJ.
#' @param calcJ_type The calcJ_type in use (NULL if calcJ was not supplied).
#' @param evalCalcR Function: calcR evaluated at a parameter vector in x order.
#' @param evalCalcJ Function: calcJ evaluated at a parameter vector in x order.
#' @param orient_calcJ Function: checks a calcJ result and returns it as p x n.
#' @param silent Logical.
#'
#' @return List: the model object, without its class.
#' @noRd
bgw_startModel <- function(x, beta, betaStart, Prob, iv, v, n, p, bgw_settings, hasAnalyticGrad,
                           calcJ_type, evalCalcR, evalCalcJ, orient_calcJ, silent) {
  # The scale vector sets the finite-difference step size, as at the start of an estimation.
  d_index_v <- iv[bgw_iv_index$d]:(iv[bgw_iv_index$d] + p - 1)
  if (bgw_settings$scalingMethod == "userScaling") {
    v[d_index_v] <- bgw_settings$userScaleVector
  } else {
    v[d_index_v] <- 1
  }

  # G is n x p, and row i is d log(R_i) / d beta.
  G <- bgw_indLevelGradients(x, beta, v, d_index_v, n, p, hasAnalyticGrad, calcJ_type, evalCalcR,
                             evalCalcJ, orient_calcJ, silent)

  hessianMethodUsed <- "None"
  recipCond         <- -1
  vcVec  <- NULL
  varcov <- NULL
  se     <- NULL
  tstat  <- NULL
  if (bgw_settings$vcHessianMethod == "bhhh") {
    # The test for singularity is the Fortran's: the covariance is computed only if the reciprocal
    # condition number exceeds machine epsilon. Here that number is exact rather than estimated.
    recipCond <- 0
    singular  <- TRUE
    if (all(is.finite(G))) {
      bhhh <- crossprod(G)
      e    <- eigen(bhhh, symmetric = TRUE, only.values = TRUE)$values
      if (max(e) > 0) recipCond <- max(min(e), 0)/max(e)
      singular <- (min(e) <= 0) || (recipCond <= .Machine$double.eps)
    }
    if (!singular) {
      varcov <- tryCatch(chol2inv(chol(bhhh)), error = function(e) NULL)
    }
    if (is.null(varcov)) {
      hessianMethodUsed <- "N/A - BHHH Hessian is singular"
    } else {
      hessianMethodUsed <- "bhhh"
      dimnames(varcov)  <- list(names(beta), names(beta))
      vcVec             <- varcov[upper.tri(varcov, diag = TRUE)]   # lower triangle, by rows
      se                <- sqrt(diag(varcov))
      names(se)         <- names(beta)
      tstat             <- beta/se
    }
  }

  requestedBy <- "userRequestModelObject = TRUE"
  if (bgw_settings$maxIterations < 1) requestedBy <- "maxIterations = 0"

  model <- list()
  model$bgw_settings    <- bgw_settings
  model$hasAnalyticGrad <- hasAnalyticGrad
  model$betaStart       <- betaStart
  model$numParams       <- p
  model$numResids       <- n
  model$code            <- 20
  model$message         <- paste0("Model object returned at betaStart without estimation (",
                                  requestedBy, ")")
  model$betaStop        <- beta
  model$finalLL         <- sum(log(as.numeric(Prob)))
  model$iterations      <- 0
  model$functionEvals   <- 1
  model$gradient        <- unname(-colSums(G))
  model$scaleVec        <- v[d_index_v]
  model$estimate        <- beta
  model$maximum         <- model$finalLL
  model$vcHessianRecipCondBGW     <- recipCond
  model$vcHessianConditionNumber  <- recipCond   # deprecated name; see bgw_buildModel
  model$hessianMethodAttemptedBGW <- bgw_settings$vcHessianMethod
  model$hessianMethodUsedBGW      <- hessianMethodUsed
  model$vcVecBGW                  <- vcVec
  model$varcovBGW                 <- varcov
  model$seBGW                     <- se
  model$tstatBGW                  <- tstat
  if (bgw_settings$userRequestIndLevelGradients) model$indLevelGradients <- G

  if (!silent) {
    cat("\n", model$message, "\n", sep = "")
    cat("Log-likelihood at betaStart: ", formatC(model$finalLL, format = "f", digits = 6), "\n",
        sep = "")
  }

  return(model)
}
