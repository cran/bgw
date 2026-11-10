#' @title bgw_itsum
#'
#' @description Prints iteration summary, info on initial and final x.
#'
#' @param d  Scaling vector
#' @param g  Gradient of objective function (negative-log-likelihood) wrt x
#' @param iv BGW internal vector of integer values
#' @param v  BGW internal vector of numeric values
#' @param x  Parameter vector for which the objective function is being minimized
#' @param p  Dimension of x
#' @param betaIsNamed Logical.
#' @param betaNames Character vector. If available, has beta parameter names.
#' @param i_itsum Code from caller to select specific options
#'
#' @return iv There are iv values changed inside bgw_itsum. The iv vector is the integer workspace
#'   for Fortran BGW.
#' @noRd
# 7 November 2026: internal function (not exported; no help page is generated).
##------------------------------------------------------------------------
# November 1, 2022
# The initial version(s) of bgw_itsum were focused on replicating (to the
# degree possible) the same functionality as ditsum in Fortran BGW.
##------------------------------------------------------------------------
bgw_itsum <- function(d, g, iv, v, x, p, betaIsNamed, betaNames = NULL, i_itsum) {
  covmat <- 26
  niter <- 31
  ##------------------------------------------------------------------------
  modelS <- list()
  modelS[[1]] <- '  G  '
  modelS[[2]] <- '  S  '
  modelS[[3]] <- ' G-S '
  modelS[[4]] <- ' S-G '
  modelS[[5]] <- 'G-S-G'
  modelS[[6]] <- 'S-G-S'
  ##------------------------------------------------------------------------
  # Convergence codes and messages (see bgw_constants.R)
  bgw_message_table <- bgw_messages(annotated = FALSE)

  ##------------------------------------------------------------------------
  # if (!betaIsNamed) betaNames <- paste0("c", 1:p)

  goto120 <- FALSE
  goto390 <- FALSE
  goto430 <- FALSE
  goto460 <- FALSE
  #
  # The following is a placeholder to indicate what
  # happens in Fortran DITSUM
  # if (iv[prunit] == 0) return(iv)
  #
  iv1 <- iv[1]
  ol <- iv[bgw_iv_index$outlev] # outlev = 19
  if (iv1 > 62) iv1 <- iv1 - 51
  if ((iv1 < 2) | (iv1 > 15)) {
    # We need to document and characterize this in greater detail.
    print(" Possible internal error...")
    print(paste(" ***** iv[1] = ", iv[1]))
    return(iv)
  }

  if (iv1 < 12) {
    # Occurs: During cases after a successful start.
    # 1. During standard iterations (iv[1]-iv[2]).
    # 2. If some type of convergence has occurred (iv[3]-iv[8])).
    # 3. If the search is being stopped (iv[9]-iv[11]).
    #
    if ((iv1 == 2) && (iv[niter] == 0)) goto390 <- TRUE
    if ((ol == 0) && !goto390) goto120 <- TRUE
    if ((iv1 >= 10) && (iv[bgw_iv_index$prntit] == 0) && !goto390) goto120 <- TRUE
    if ((!goto390) && (!goto120)) {
      if (iv1 <= 2) {
        iv[bgw_iv_index$prntit] <- iv[bgw_iv_index$prntit] + 1
        if (iv[bgw_iv_index$prntit] < abs(ol)) return(iv)
      }
      # Line 10
      nf <- iv[bgw_iv_index$nfcall] - abs(iv[bgw_iv_index$nfcov])
      iv[bgw_iv_index$prntit] <- 0
      reldf <- 0.
      preldf <- 0.
      oldf <- max(abs(v[bgw_v_index$f0]), abs(v[bgw_v_index$f]))
      if (oldf > 0.) {
        reldf <- v[bgw_v_index$fdif]/oldf
        preldf <- v[bgw_v_index$preduc]/oldf
      }
      nreldf <- 0.
      if (oldf > 0.) nreldf <- v[bgw_v_index$nreduc]/oldf
      m <- iv[bgw_iv_index$sused]
      # WRITE(pu,100) stuff
      # DSB NOTE: We have not figured out why the last iteration information is printed TWICE!
      # print(iv1)
      if (iv[1] == 2) {
        if (ol > 0) cat(format(iv[niter], nsmall = 0, trim = FALSE, width = 6), format(nf,
                                                                                       nsmall = 0,
                                                                                       trim = FALSE,
                                                                                       width = 5),
                        formatC(v[bgw_v_index$f], format = "e", digits = 9),
                        formatC(reldf, format = "e", digits = 3),
                        formatC(preldf, format = "e", digits = 3),
                        formatC(v[bgw_v_index$reldx], format = "e", digits = 2),
                        modelS[[m]],
                        formatC(v[bgw_v_index$stppar], format = "e", digits = 2),
                        formatC(v[bgw_v_index$dstnrm], format = "e", digits = 2),
                        formatC(nreldf, format = "e", digits = 3))
        if (ol < 0) cat(format(iv[niter], nsmall = 0, trim = FALSE, width = 6), format(nf,
                                                                                       nsmall = 0,
                                                                                       trim = FALSE,
                                                                                       width = 5),
                        formatC(v[bgw_v_index$f], format = "e", digits = 9),
                        formatC(reldf, format = "e", digits = 3),
                        formatC(preldf, format = "e", digits = 3),
                        formatC(v[bgw_v_index$reldx], format = "e", digits = 2),
                        modelS[[m]],
                        formatC(v[bgw_v_index$stppar], format = "e", digits = 2))
        # Original FORMAT(i6,i5,e10.3,2E9.2,e8.1,a3,a4,2E8.1,e9.2)
        # FORMAT(i6,i5,e15.8,2E9.2,e8.1,a3,a4,2E8.1,e9.2)
        cat("\n")
      }
      goto120 <- TRUE
    }
  } else {
    goto120 <- TRUE
  }
  #
  if (goto120) {
    if (iv1 <= 2) return(iv)
    i <- iv[bgw_iv_index$statpr]
    if ((i == -1) | (i + iv1 < 0)) {
      goto460 <- TRUE
    } else {
      cat("\n")
      if (iv1 == 1) return(iv)
      if (iv1 == 2) return(iv)
      if ((iv1 >= 3) && (iv1 <= 15)) {
        cat(bgw_message_table[[iv1]])
      }
      if ((iv1 >= 3) && (iv1 <= 11)) goto430 <- TRUE
      if (iv1 == 12) goto390 <- TRUE
      if ((iv1 == 13) | (iv1 == 15)) return(iv)
      if (iv1 == 14) {
        if (iv[niter] > 0) { goto460 <- TRUE
        } else { goto390 <- TRUE
        }
      }
    }
  }
  #
  if (goto390) {
    if (iv[bgw_iv_index$x0prt] != 0) {
      cat("\n")
      # df <- data.frame('XName' = betaNames,
      #                  'Initial_X_i' = x, 'D_i' = d)
      # colnames(df) <- c('XName','Initial X(i)','D(i)')
      # writeLines(paste0("    ",capture.output(print(df))))
      # if (betaIsNamed) {
      #     df <- data.frame('XName' = betaNames,
      #                  'Initial_X_i' = x)
      #     colnames(df) <- c('XName','Initial X(i)')
      # } else {
      #   df <- data.frame('Initial_X_i' = x)
      #   colnames(df) <- c('Initial X(i)')
      # }
      if (betaIsNamed) {
        df <- data.frame('BetaName' = betaNames, 'Initial_X_i' = x, 'D_i' = d)
        colnames(df) <- c('BetaName', 'InitialBeta(i)', 'D(i)')
      } else {
        df <- data.frame('Initial_X_i' = x, 'D_i' = d)
        colnames(df) <- c('InitialBeta(i)', 'D(i)')
      }
      writeLines(paste0("    ", capture.output(print(df))))
      rm(df)
      cat("\n")
    }
    v[bgw_v_index$dstnrm] <- 0.
    v[bgw_v_index$fdif] <- 0.
    v[bgw_v_index$nreduc] <- 0.
    v[bgw_v_index$preduc] <- 0.
    v[bgw_v_index$reldx] <- 0.
    if (iv1 >= 12) return(iv)
    iv[bgw_iv_index$prntit] <- 0
    if (ol == 0) return(iv)
    cat("\n")
    if (ol > 0) cat(paste0("    it    nf     F            RELDF    PRELDF    RELDX    MODEL ",
                           "stppar   D*step   NPRELDF"))
    if (ol < 0) cat("    it    nf     F            RELDF    PRELDF    RELDX    MODEL stppar")
    # cat("    it    nf     F            RELDF    PRELDF    RELDX    MODEL stppar   D*step
    #          NPRELDF")
    cat("\n")
    cat(format(iv[niter], nsmall = 0, trim = FALSE, width = 6), format(iv[bgw_iv_index$nfcall],
                                                                       nsmall = 0,
                                                                       trim = FALSE, width = 5),
        formatC(v[bgw_v_index$f], format = "e", digits = 9))
    cat("\n")
  }

  if (goto430) {
    goto460 <- TRUE
    # goto430
    if (iv[bgw_iv_index$statpr] > 0) {
      oldf <- max(abs(v[bgw_v_index$f0]), abs(v[bgw_v_index$f]))
      preldf <- 0.
      nreldf <- 0.
      if (oldf > 0.) {
        preldf <- v[bgw_v_index$preduc]/oldf
        npreldf <- v[bgw_v_index$nreduc]/oldf
      }
      nf <- iv[bgw_iv_index$nfcall] - iv[bgw_iv_index$nfcov]
      ng <- iv[bgw_iv_index$ngcall] - iv[bgw_iv_index$ngcov]
      # if (i_itsum == 1) {
      cat("\n")
      cat("       FUNCTION    ", formatC(v[bgw_v_index$f], format = "e", digits = 9),
          " RELDX       ", formatC(v[bgw_v_index$reldx], format = "e", digits = 3),
          "\n       PRELDF      ", formatC(preldf, format = "e", digits = 3),
          "       NPRELDF     ", formatC(npreldf, format = "e", digits = 3),
          "\n       func. evals ", format(nf, digits = 8),
          "              grad. evals ", format(ng, digits = 8), "\n")
      #}
    }
  }

  if (goto460) {
    # goto460
    if (iv[bgw_iv_index$solprt] == 0) {
      cat("\n")
      return(iv)
    }
    # iv[covmat] = iv[26] = location of variance-covariance matrix
    if (iv[covmat] > 0) {
      cat("\n")
      if (iv[bgw_iv_index$covreq] == 3)  cat('       vcHessianMethod = Gauss-Newton/BHHH \n')
      if (iv[bgw_iv_index$covreq] == 2)  cat(paste0('       vcHessianMethod = finite differences ',
                                                    '(gradient evals) \n'))
      if (iv[bgw_iv_index$covreq] == -2) cat(paste0('       vcHessianMethod = finite differences ',
                                                    '(function evals) \n'))
      cat('       Estimated upper bound on reciprocal of Euclidean condition number of vcHessian: ',
          v[bgw_v_index$rcond]^2, ' \n')
      cat(paste0('       (Reciprocal condition number = ratio of smallest singular value to ',
                 'largest singular value.)\n'))
      cat('       Value of unit roundoff = machine epsilon for comparison purposes:               ',
          .Machine$double.eps, ' \n')
      # cat('       Value of iv[fdh] = ',iv[fdh],'\n')
      # cat('       Value of iv[hc]  = ',iv[hc],'\n')
      # cat('       Value of iv[h]   = ',iv[h],'\n')


      vc_start <- iv[bgw_iv_index$covmat]   # position in v where the covariance matrix starts
      vc_vec   <- v[vc_start:(vc_start + p*(p+1)/2 - 1)]
      vc <- matrix(0, p, p)
      iii <- 0
      for (ii in 1:p) {
        for (jj in 1:ii) {
          iii <- iii + 1
          vc[ii, jj] <- vc_vec[iii]
          vc[jj, ii] <- vc_vec[iii]
        }
      }
      cat("\n")
      se <- rep(0, p)
      tstat <- rep(0, p)
      for (i in 1:p) {
        se[i] <- sqrt(vc[i, i])
        tstat[i] <- x[i]/se[i]
      }
      if (betaIsNamed) {
        df <- data.frame('XName' = betaNames,
                         'FinalX(i)' = x, 's.e.(i)' = se, 't.rat.(0)' = tstat, 'G(i)' = g,
                         'D(i)' = d)
        colnames(df) <- c('BetaName', 'FinalBeta(i)', 's.e.(i)', 't.rat.(0)', 'G(i)', 'D(i)')

      } else {
        df <- data.frame('FinalX(i)' = x, 's.e.(i)' = se, 't.rat.(0)' = tstat, 'G(i)' = g,
                         'D(i)' = d)
        colnames(df) <- c('FinalBeta(i)', 's.e.(i)', 't.rat.(0)', 'G(i)', 'D(i)')
      }
      writeLines(paste0("    ", capture.output(print(df))))
      rm(se, tstat, df)
    } else {
      # df <- data.frame('XName' = betaNames,
      #                  'FinalX(i)' = x, 'D(i)' = d, 'G(i)' = g)
      # colnames(df) <- c('XName','FinalX(i)','D(i)','G(i)')
      # writeLines(paste0("    ",capture.output(print(df))))
      # rm(df)
      if (iv[bgw_iv_index$rdreq] == 0) {
        cat('       No variance-covariance matrix computation was requested. \n')
        cat("\n")
      } else {
        cat("\n")
        if (iv[bgw_iv_index$covreq] == 3) {
          cat('       Requested vcHessianMethod = Gauss-Newton/BHHH \n')
        }
        if (iv[bgw_iv_index$covreq] == 2) cat(paste0('       Requested vcHessianMethod = finite ',
                                                     'differences ',
                                                     '(gradient evals) \n'))
        if (iv[bgw_iv_index$covreq] == -2) cat(paste0('      Requested vcHessianMethod = finite ',
                                                      'differences ',
                                                      '(function evals) \n'))
        if (iv[covmat] == -1) {
          cat('       ***** Unable to compute a variance-covariance matrix... ***** \n')
          cat('       ***** Requested vcHessian matrix was indefinite         *****  \n')
        }
        if (iv[covmat] == -2) {
          cat(paste0('      ***** Unable to compute a variance-covariance matrix....              ',
                     '                            ***** \n'))
          cat(paste0('      ***** Too many beta values were rejected during finite-difference ',
                     'computation of vcHessian matrix *****  \n'))
        }
        cat(paste0('       Estimated upper bound on reciprocal of Euclidean condition number of ',
                   'vcHessian: '), v[bgw_v_index$rcond]^2, ' \n')
        cat(paste0('       (Reciprocal condition number = ratio of smallest singular value to ',
                   'largest singular value.)\n'))
        cat(paste0('       Value of unit roundoff = machine epsilon for comparison purposes:      ',
                   '         '), .Machine$double.eps, ' \n')
        # cat('       Value of iv[fdh] = ',iv[fdh],'\n')
        # cat('       Value of iv[hc]  = ',iv[hc],'\n')
        # cat('       Value of iv[h]   = ',iv[h],'\n')
        cat("\n")
      }

      if (betaIsNamed) {
        df <- data.frame('XName' = betaNames,
                         'FinalX(i)' = x, 'G(i)' = g, 'D(i)' = d)
        colnames(df) <- c('BetaName', 'FinalBeta(i)', 'G(i)', 'D(i)')
      } else {
        df <- data.frame('FinalX(i)' = x, 'G(i)' = g, 'D(i)' = d)
        colnames(df) <- c('FinalBeta(i)', 'G(i)', 'D(i)')
      }
      writeLines(paste0("    ", capture.output(print(df))))
      rm(df)
      # tmp <- c(x,d,g)
      # # tmp <- matrix(tmp, nrow = length(tmp), ncol = 3, dimnames=list(betaNames,c('Final
      # X(i)','D(i)','G(i)')))
      # tmp <- matrix(tmp,nrow = length(tmp), ncol = 3)
      # print(tmp)
      # rm(tmp)
      # output = cbind(c(1:p),x,d,g)
      # colnames(output) <- c('i','Final X(i)','D(i)','G(i)')
      # rownames(output) <- betaNames
      # apollo_print(output) #print(output, digits = 6)
    }
    cat("\n")
  }
  return(iv)
}
