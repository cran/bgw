#' @title bgw_mle
#'
#' @description Performs maximum likelihood estimation (MLE) for the user-provided model
#' defined in \code{calcR}. Version 0.1.5.
#'
#' @details This function has been written to provide an R-based interface to Fortran estimation
#' software published in Bunch, Gay and Welsch (1993), "Algorithm 717-Subroutines for Maximum
#' Likelihood and Quasi-Likelihood Estimation of Parameters in Nonlinear Regression Models," ACM
#' Transactions on Mathematical Software, 19 (1), March 1993, 109-130. The letters BGW will be used
#' in various ways to denote the source of the estimation functionality.
#'
#' A primary motivation was to develop a more efficient maximum likelihood estimation function for
#' use in the Apollo choice modelling package: see \url{https://www.apollochoicemodelling.com/}.
#' However, we have adopted a design whereby the BGW package is wholly independent of Apollo, and
#' can be used in a stand-alone fashion. Note also that the BGW Fortran subroutines are written to
#' support general statistical estimation for an arbitrary objective/criterion function. So,
#' although this version of the package is specifically written for MLE, the package may see future
#' updates that expand the number of estimation options (for, e.g., nonlinear least squares,
#' generalized method of moments, etc.).
#'
#' Remark: Following the convention in the numerical optimization literature, BGW minimizes the
#' objective function. That is, bgw_mle minimizes the negative-log-likelihood for the model calcR.
#'
#' \subsection{Forms of the Jacobian (calcJ)}{
#' calcJ must return a matrix (not a vector), in one of the following forms.
#'
#' \strong{Derivatives of probabilities.} Formally, the Jacobian (J) is an n x p matrix whose [i,k]
#' element is dR_i/dbeta_k. The equivalent information can be provided as a p x n matrix (J
#' transpose). Because these are derivatives of probabilities, if J is provided in either form, the
#' user should also set \code{calcJ_type = "probability"}.
#'
#' \strong{Derivatives of log-probabilities.} Alternatively, J can be provided as an n x p matrix
#' whose [i,k] element is d[log(R_i)]/dbeta_k, or, equivalently, the p x n matrix that is its
#' transpose. In either case the user should set \code{calcJ_type = "log_probability"}.
#'
#' \strong{Orientation.} bgw identifies from the dimensions whether J is n x p or p x n. In the rare
#' case where n = p, the dimensions cannot show this and the matrix is read as n x p, so a square
#' matrix must be supplied in the n x p form.
#'
#' \strong{Version 0.1.4 and earlier.} In earlier versions of bgw, the assumed form for J was the
#' n x p matrix with [i,k] element d[log(R_i)]/dbeta_k. To maintain backward compatibility, if
#' calcJ_type is missing, bgw assumes that J contains derivatives of log(R_i). Because this cannot
#' be verified from the matrix itself, a calcJ that returns derivatives of probabilities must be
#' accompanied by \code{calcJ_type = "probability"}; otherwise the results will be wrong.
#' }
#'
## Arguments for bgw_mle
##
## calcR       Function that computes an n-vector of model residuals (R) for a p-vector of
##                parameters beta.
##                In this case the residuals are likelihoods (probabilities).
##                The first argument of calcR must be the parameter vector beta.
## calcJ       User-provided function that computes the derivatives of R wrt beta, as a matrix,
##                either n x p (one row per observation) or p x n. If NULL, finite-difference
##                derivatives are used.
## calcJ_type  What calcJ returns: "probability" (derivatives of R) or "log_probability"
##                (derivatives of log(R)). If calcJ is supplied without calcJ_type,
##                "log_probability" is assumed (as in earlier versions of bgw).
## start      initial vector of starting values for parameters.
##
## bgw_settings Control parameters for BGW algorithm.
#'
#' @param calcR Function that computes an n-vector (R) of model residuals for a p-vector of
#'   (numeric) parameters beta (the first argument). (The beta vector can be named or unnamed.)
#'   Because bgw_mle performs maximum likelihood estimation (MLE), R_i is the likelihood
#'   (probability) of the ith observation, and the estimate is the value of beta that minimizes the
#'   negative log-likelihood NLL = -Sum_i [log(R_i)].
#' @param betaStart Vector of initial starting values for beta. Can be either a named or unnamed
#'   vector.
#' @param calcJ Optional function that computes analytic derivative (Jacobian) information for the
#'   n-vector R, with beta as its first argument, as for calcR. If calcJ is NULL (the default),
#'   estimation uses finite-difference derivatives. Starting with version 0.1.5, the Jacobian can be
#'   supplied in several forms (see Details), so whenever calcJ is supplied, calcJ_type should be
#'   supplied as well.
#' @param calcJ_type Character. Should be supplied whenever calcJ is supplied: \code{"probability"}
#'   if calcJ returns derivatives of R_i, or \code{"log_probability"} if it returns derivatives of
#'   log(R_i). See Details. If omitted, \code{"log_probability"} is assumed (see "Version 0.1.4 and
#'   earlier" under Details). Ignored if calcJ is NULL.
#' @param bgw_settings List. Contains control parameters for BGW estimation code. All parameters
#'   have default values, so user input is entirely optional.
#'                \itemize{
#'                  \item \strong{\code{printLevel}}: Integer (0-3). Controls the level of detail
#'                    for iteration information written to the console during estimation. 0 =
#'                    silent, 1 = print starting values and final solution only, 2 = short line
#'                    summary, 3 = long line summary (see documentation). Default = 3.
#'                  \item \strong{\code{silent}}: Logical. Suppresses all output. Default = FALSE.
#'                    (Currently redundant with printLevel = 0.)
#'
#'                  \item \strong{\code{printNonDefaultSettings}}: Logical. Echos any user-provided
#'                    bgw_settings. Default = TRUE.
#'                  \item \strong{\code{printStartingValues}}: Logical. Default = TRUE. Print
#'                    starting values for beta.
#'                  \item \strong{\code{printFinalResults}}: Logical. Default = TRUE. Print final
#'                    status of estimation (convergence/error message), summary statistics (negative
#'                    log-likelihood, iterations, function evals, etc.), beta, gradient, and (if
#'                    available) estimated standard errors and t-ratios. (See vcHessianMethod for
#'                    variance-covariance specification.)
#'                  \item \strong{\code{maxIterations}}: Numeric. Maximum number of iterations for
#'                    the estimation. Default is 200. A value of 0 requests a model object at
#'                    betaStart without estimation (see userRequestModelObject).
#'                  \item \strong{\code{maxFunctionEvals}}: Numeric. Maximum number of objective
#'                    function evaluations. Default is 300.
#'                  \item \strong{\code{modelName}}: Character. Used to create names for requested
#'                    output files. Default is "bgw_mle_model."
#'                  \item \strong{\code{outputDirectory}}: Character. Directory in which to write
#'                    output files. It is created if it does not exist. Default is the current
#'                    working directory.
#'                  \item \strong{\code{writeItSummary}}: Logical. If TRUE, iteration information is
#'                    echoed in the output file "modelName_itSummary.csv". Default is FALSE. [Not
#'                    currently implemented.]
#'                  \item \strong{\code{writeIter}}: Logical. If TRUE, parameters and log-likelihood
#'                    for each iteration are written to "modelName_iterations.csv". The columns are
#'                    the parameter names (or beta1, beta2, ... if betaStart has no names), logLike,
#'                    Iteration and ItType. An existing file is replaced, unless writeIterMode is
#'                    "append". Default is FALSE.
#'                  \item \strong{\code{writeIterMode}}: Character. "overWrite" (the default)
#'                    replaces any existing iteration files; "append" adds to them. Applies to
#'                    writeIter, writeIterGradient and writeIterScaleVec.
#'                  \item \strong{\code{writeIterGradient}}: Logical. If TRUE, the gradient of the
#'                    log-likelihood at each iteration is written to "modelName_gradients.csv", with
#'                    the same columns as for writeIter. Its sign is opposite to that of the
#'                    returned gradient, which is the gradient of the negative log-likelihood.
#'                    Default is FALSE.
#'                  \item \strong{\code{writeIterScaleVec}}: Logical. If TRUE, the scale vector at
#'                    each iteration is written to "modelName_scaleVecs.csv", with the same columns
#'                    as for writeIter. Default is FALSE.
#'                  \item \strong{\code{vcHessianMethod}}: Character. Method for computing the
#'                    Hessian approximation used for the variance-covariance matrix (VC = H^(-1)).
#'                    Options are: "none","bhhh","finiteDifferences", or "fdFunction"
#'                    ("finiteDifferences" automatically uses gradient differences if a gradient is
#'                    available, otherwise it uses objective function differences. "fdFunction"
#'                    allows the user to use objective function differences even if a gradient is
#'                    available). Default is "bhhh."
#'                  \item \strong{\code{scalingMethod}}: Character. Method used to compute a scaling
#'                    vector (scaleVec_i, i=1,..,p). Define a matrix D = diag(scaleVec_1,...,
#'                    scaleVec_p). When using scaling, values in scaleVec should be chosen so that
#'                    the elements of D*beta are roughly comparable in size. The re-scaled beta is
#'                    used when computing trial steps using trust regions, and when computing
#'                    stopping criteria. Options are: "adaptive," "none," and "userScaling": the
#'                    default is "adaptive." For a description of the "adaptive" method, which
#'                    updates scaleVec at each iteration using information from the model Jacobian,
#'                    see Bunch, Gay, and Welsch (1993). For "none," scaleVec is set to a p-vector
#'                    of ones for the entire search. The "userScaling" option indicates that the
#'                    user is supplying their own (fixed) scaleVec in
#'                    bgw_settings[["userScaleVector"]] (see next item). Both items must be properly
#'                    set or an error occurs.
#'                  \item \strong{\code{userScaleVector}}: Numeric, with dimension p = number of
#'                    free parameters. Can be either a named or unnamed vector. This is a
#'                    user-provided scaling vector that is used ONLY in conjunction with the
#'                    non-default option "userScaling" for bgw_settings[["scalingMethod"]] (see
#'                    previous item).
#'                  \item\strong{\code{userRequestIndLevelGradients}}: Logical. If TRUE, BGW returns
#'                    an n x p matrix containing individual-level gradients at the final value of
#'                    beta. Specifically, if the log-likelihood = sum_i[log(R_i)], row i of the
#'                    matrix is d log(R_i)/d beta. Available with both analytic (calcJ) and
#'                    finite-difference derivatives. Default is FALSE.
#'                  \item\strong{\code{userRequestModelObject}}: Logical. If TRUE, no estimation is
#'                    performed. BGW returns a model object evaluated at betaStart, with code 20,
#'                    iterations equal to 0, and a variance-covariance matrix based on the BHHH
#'                    matrix (NULL if that matrix is singular). Requires vcHessianMethod = "bhhh"
#'                    or "none". Default is FALSE.
#'                }
#' @return model object of class 'bgw_mle'. Output of a bgw maximum likelihood estimation procedure.
#'   A list with the following attributes:
#' \itemize{
#' \item \strong{\code{betaStart}}: Vector of initial starting values.
#' \item \strong{\code{bgw_settings}}: List. The final settings used, including any defaults.
#' \item \strong{\code{hasAnalyticGrad}}: Logical. Indicates if an analytical gradient calculation
#'   was used. If the user has not provided a calcJ function (see input parameter), it is set to
#'   FALSE.
#' \item \strong{\code{numParams}}: Numeric. Number of model parameters used in calcR.
#' \item \strong{\code{numResids}}: Numeric. Number of independent observations (model residuals) in
#'   data set = dimension of calcR output.
#' \item \strong{\code{code}}: Numeric. Numeric return code from BGW. Code 20 means that a model
#'   object at betaStart was requested, and no estimation was performed. Codes above 62 report a
#'   failure: 63, 64 and 66 at the start, and 65 when the gradient could not be computed. If code
#'   65 occurs after the start, the last iterate at which the gradient could be computed is
#'   returned, with its log-likelihood, gradient and scale vector.
#' \item \strong{\code{message}}: Character. Message statement characterizing termination of MLE
#'   search.
#' \item \strong{\code{betaStop}}: Vector. Value of parameter vector at conclusion of MLE search.
#'   See message to determine if beta is a valid estimate.
#' \item \strong{\code{finalLL}}: Numeric. Value of log-likelihood at betaStop.
#' \item \strong{\code{iterations}}: Numeric. Number of iterations used in MLE search.
#' \item \strong{\code{functionEvals}}: Numeric. Number of function evaluations used in MLE search.
#'   (This excludes function evaluations used by any finite-difference calculations for the gradient
#'   and/or the vcHessian.)
#' \item \strong{\code{gradient}}: Vector. The gradient of the negative log-likelihood, evaluated at
#'   betaStop.
#' \item \strong{\code{scaleVec}}: Vector. The scaling vector at the conclusion of the MLE search
#'   (see scalingMethod in bgw_settings).
#' \item \strong{\code{estimate}}: Vector. Parameter vector at the conclusion of the search (the
#'   same as betaStop). It is always returned; it is the MLE only if the search converged
#'   successfully (code 3, 4, 5 or 6; see code and message).
#' \item \strong{\code{maximum}}: Numeric. Log-likelihood value at estimate (the same as finalLL).
#'   It is always returned; it is the maximum of the log-likelihood only if the search converged
#'   successfully.
#' \item \strong{\code{hessianMethodAttemptedBGW}}: String. Requested method for computing vcHessian
#'   (from bgw_settings).
#' \item \strong{\code{hessianMethodUsedBGW}}: String. Method actually used for computing vcHessian.
#'   If none was requested, this is "None"; if the matrix was singular, it is "N/A - BHHH Hessian is
#'   singular".
#' \item \strong{\code{vcHessianRecipCondBGW}}: Numeric. Estimated upper bound on the reciprocal
#'   of the Euclidean condition number of vcHessian (if available); 0 corresponds to a singular
#'   matrix. Set to -1 if unavailable. For a model object at betaStart (see
#'   userRequestModelObject), this is the exact reciprocal condition number.
#' \item \strong{\code{vcHessianConditionNumber}}: Numeric. Deprecated; the same value as
#'   vcHessianRecipCondBGW (despite the name, a reciprocal). It will be removed in a future
#'   version.
#' \item \strong{\code{varcovBGW}}: Matrix. p-by-p matrix containing estimate of the
#'   variance-covariance matrix (if requested and available).
#' \item \strong{\code{vcVecBGW}}: Vector. Lower triangle of variance-covariance matrix stored in
#'   vector form (row-major order, if requested and available).
#' \item \strong{\code{seBGW}}: Vector. Estimated standard errors for parameter estimates (if
#'   requested and available).
#' \item \strong{\code{tstatBGW}}: Vector. Estimated t-statistics (versus 0, if requested and
#'   available).
#' \item\strong{\code{indLevelGradients}}: Matrix. n-by-p matrix whose row i is d log(R_i)/d beta,
#'   computed at the final value of beta (if requested via userRequestIndLevelGradients).
#' }
#' @export

bgw_mle <- function(calcR, betaStart, calcJ = NULL, bgw_settings = NULL, calcJ_type = NULL) {
  ### Parameters:
  ## calcR       function that computes an n-vector of model residuals for a p-vector of parameters
  ##             x.
  ##                In this case the residuals are likelihoods (probabilities).
  ##                The first argument must be the parameter vector x.
  ## calcJ      derivatives wrt x, n x p or p x n. If NULL, finite-difference
  ##                derivatives are used.
  ## calcJ_type "probability" (calcJ returns derivatives of R) or "log_probability"
  ##                (derivatives of log(R)). Added as the LAST argument so that existing
  ##                calls that pass arguments by position are unaffected.
  ## start      initial vector of starting values for parameters
  ##
  ## bgw_settings control parameters for BGW algorithm.
  ## Current working list of control parameters:
  #
  # printLevel
  # maxIterations
  # maxFunctionEvaluations
  # modelName
  # outputDirectory
  # writeIter
  # writeItSummary
  # vcHessianMethod.  Options are:  none, bhhh, finiteDifference, fdFunction
  #
  # scalingMethod  [Currently fixed to default d=1]
  #
  # Possible convergence tolerance parameters
  # RelativeFunction_ConvergeTol
  # Step_ConvergeTol
  # Singular_ConvergeTol
  # False_ConvergeTol
  # AbsoluteFunction_ConvergeTol

  ##------------------------------------------------------------------------
  # #################### #
  #### Loading Inputs ####
  # #################### #

  # Initial test for betaStart and calcR
  test <- is.numeric(betaStart) && is.vector(betaStart) && is.function(calcR)
  if (!test) stop(paste0('Arguments calcR (function) and betaStart (numeric vector) must be ',
                         'provided to bge_mle.'))

  # Critical inputs to BGW are p and n.
  # p, n are established below. They are needed for the call to bgw_setup.
  # The following are comments from BGW Fortran
  #
  #  ***  PARAMETER USAGE  ***
  # N....... TOTAL NUMBER OF RESIDUALS.
  # P....... NUMBER OF PARAMETERS (COMPONENTS OF X) BEING ESTIMATED.
  # PS...... NUMBER OF NON-NUISANCE PARAMETERS (THOSE INVOLVED IN S).
  # X....... PARAMETER VECTOR BEING ESTIMATED (INPUT = INITIAL GUESS,
  #             OUTPUT = BEST VALUE FOUND).
  # RHO..... SUBROUTINE FOR COMPUTING LOSS FUNCTIONS AND THEIR DERIVS.
  #             SEE  DRGLG FOR DETAILS ABOUT RHO.
  # RHOI.... PASSED WITHOUT CHANGE TO RHO.
  # RHOR.... PASSED WITHOUT CHANGE TO RHO.
  # IV...... INTEGER VALUES ARRAY.
  # LIV..... LENGTH OF IV, AT LEAST 90 + P.
  #
  # As indicated above, BGW allows parameters to be of two
  # types: nuisance and non-nuisance
  # p is the total number of parameters
  # ps is the number of non-nuisance parameters
  # Nuisance parameters are not implemented here,
  # so p = ps.

  p <- length(betaStart)
  ps <- p

  # DSB note: Could change this to be included in bgw_settings.
  # Right now it is changed manually and recompiled.
  debug <- FALSE
  if (debug) {
    cat("\nChecking inputs to bgw_mle\n")
    cat("\nChecking betaStart: \n")
    print(betaStart)
    cat("\n")
    cat("\nLength of betaStart = p\n")
    print(p)
    cat("\n")
  }

  beta <- betaStart
  # Establishing n requires a successful call to calcR
  Prob <- calcR(beta)
  if (!is.numeric(Prob)) stop(" Initial call to calcR returned a non-numeric result.")
  # MLE ONLY: R_i <= 0 is rejected here because log(R_i) is undefined. In the general BGW framework
  # R_i is a generalized residual, and zero is a valid value (for example in nonlinear least
  # squares). If bgw is extended to other estimators, this test must be changed.
  if (any(!is.finite(Prob) | Prob <= 0)) {
    stop(" Initial call to calcR failed: it returned NA, NaN, Inf, or a value <= 0.")
  }
  if (!is.vector(Prob))  stop(" Initial call to calcR returned non-vector result.")
  n <- length(Prob)

  if (debug) {
    cat("\nNumber of residuals produced by calcR = n\n")
    print(n)
    cat("\n")
    cat("\nFirst ten choice probabilities\n")
    print(Prob[1:10])
    stop('Stop due to debugging...')
  }

  betaNames <- names(betaStart)
  if (is.null(betaNames)) {
    betaIsNamed <- FALSE
    x <- beta
  } else {
    betaIsNamed <- TRUE
    copy_BetaToX <- function(beta, p) {
      x <- rep(0, p)
      for (i in 1:p) {x[i] <- beta[[i]]}
      return(x)
    }
    copy_XToBeta <- function(x, beta, p) {
      for (i in 1:p) {beta[[i]] <- x[i]}
      return(beta)
    }
    x <- copy_BetaToX(beta, p)
  }

  # Evaluate calcR and calcJ at a parameter vector in x order, passing a named copy when betaStart
  # has names. Because a copy is passed, finite-difference perturbations never alter beta.
  evalCalcR <- function(xx) {
    if (betaIsNamed) calcR(copy_XToBeta(xx, betaStart, p)) else calcR(xx)
  }
  evalCalcJ <- function(xx) {
    if (betaIsNamed) calcJ(copy_XToBeta(xx, betaStart, p)) else calcJ(xx)
  }

  # This needs to be error checked next.
  # But, echoing valid information to user is done later.
  if (!is.null(calcJ)) {
    if (is.function(calcJ)) {
      hasAnalyticGrad <- TRUE
    } else {
      stop("\nArgument calcJ provided to bgw_mle must be a function.\n")
    }
  } else {
    hasAnalyticGrad <- FALSE
  }

  # bgw_setup allocates space for iv and v, and sets up their default values
  bgw_setup <- bgw_mle_setup(p, n, hasAnalyticGrad)
  iv  <- bgw_setup$iv_r
  liv <- length(iv)
  v   <- bgw_setup$v_r
  lv  <- length(v)
  debug <- FALSE
  if (debug) write.csv(iv, "bgw_iv_r_setup.csv")
  if (debug) write.csv(v, "bgw_v_r_setup.csv")

  # 7 November 2026: calcJ_type argument. Earlier versions of bgw (and all
  # versions of Apollo to date) assume that calcJ returns an n x p matrix of
  # derivatives of log(R), and do not pass calcJ_type. To keep such callers
  # working, a missing calcJ_type means "log_probability". No warning is issued:
  # users of an existing Apollo could not act on it, and Apollo calls bgw_mle
  # repeatedly in some cases.
  if (hasAnalyticGrad) {
    if (is.null(calcJ_type)) {
      calcJ_type <- "log_probability"
    } else if (!(is.character(calcJ_type) && length(calcJ_type) == 1 &&
                 calcJ_type %in% c("probability", "log_probability"))) {
                   stop('\nArgument calcJ_type must be "probability" or "log_probability".\n')
                 }
  }

  # Settings: what can be set, what the caller asked for, and the resulting values of iv and v.
  settingDefs  <- bgw_defineSettings()
  resolved     <- bgw_resolveSettings(bgw_settings, settingDefs, hasAnalyticGrad, p)
  bgw_settings <- resolved$bgw_settings
  applied      <- bgw_applySettings(iv, v, settingDefs$default, bgw_settings,
                                    resolved$modifiedSettings, resolved$silent, hasAnalyticGrad)
  iv           <- applied$iv
  v            <- applied$v
  silent       <- applied$silent
  rm(settingDefs, resolved, applied)


  test <- is.null(bgw_settings$outputDirectory) || bgw_settings$outputDirectory == ''
  if (test) {
    bgw_settings$outputDirectory <- getwd()
  } else {
    if (!dir.exists(bgw_settings$outputDirectory)) {
      if (!silent) cat("\noutputDirectory provided by user does not exist, so will be created.\n")
      dir.create(bgw_settings$outputDirectory)
    }
  }

  tmp <- substr(bgw_settings$outputDirectory,
                nchar(bgw_settings$outputDirectory),
                nchar(bgw_settings$outputDirectory))
  if (tmp != "/") bgw_settings$outputDirectory <- paste0(bgw_settings$outputDirectory, "/")

  modelName         <- bgw_settings$modelName
  outputDirectory   <- bgw_settings$outputDirectory
  writeIter         <- bgw_settings$writeIter
  writeIterMode     <- bgw_settings$writeIterMode
  writeIterGradient <- bgw_settings$writeIterGradient
  writeIterScaleVec <- bgw_settings$writeIterScaleVec
  writeItSummary    <- bgw_settings$writeItSummary
  returnIndLevelGradients <- bgw_settings$userRequestIndLevelGradients

  # 7 November 2026: checks the matrix returned by calcJ and returns it as p x n
  # (one column per observation), the layout required by the Fortran. The layout
  # is read from the dimensions; when n = p, the matrix is read as n x p.
  orient_calcJ <- function(J) {
    if (!is.matrix(J)) {
      stop("\ncalcJ must return a matrix (n x p or p x n), not a vector.\n")
    }
    if (nrow(J) == n && ncol(J) == p) return(t(J))
    if (nrow(J) == p && ncol(J) == n) return(J)
    stop(sprintf(paste0("\ncalcJ returned a %d x %d matrix; expected %d x %d (n x p) ",
                        "or %d x %d (p x n).\n"),
                 nrow(J), ncol(J), n, p, p, n))
  }

  # A model object at betaStart, without estimation. The Fortran is not called on this path.
  if (bgw_settings$userRequestModelObject) {
    model <- bgw_startModel(x, beta, betaStart, Prob, iv, v, n, p, bgw_settings, hasAnalyticGrad,
                            calcJ_type, evalCalcR, evalCalcJ, orient_calcJ, silent)
    class(model) <- "bgw_mle"
    return(model)
  }

  if (writeIter) {
    # Remove modelName_iterations if it exists
    iterFile <- paste0(outputDirectory, modelName, "_iterations.csv")
    tmp <- iterFile
    if (writeIterMode == "overWrite") {
      txt <- paste0('\nCould not delete old ', tmp,
                    ' file. New iterations will be written after old ones.')
      if (file.exists(tmp)) tryCatch(file.remove(tmp), error = function(e) cat(txt))
      rm(txt)
      if (!silent) cat("\nIterates will be written to: \n", iterFile)
    }
    if (writeIterMode == "append") {
      if (!silent) cat("\nIterates will be appended to: \n", iterFile)
    }
  }

  # Column names for the iteration files.
  writeNames <- if (betaIsNamed) betaNames else paste0("beta", 1:p)

  gradFile  <- paste0(outputDirectory, modelName, "_gradients.csv")
  scaleFile <- paste0(outputDirectory, modelName, "_scaleVecs.csv")
  for (tmp in c(if (writeIterGradient) gradFile, if (writeIterScaleVec) scaleFile)) {
    if (writeIterMode == "overWrite") {
      txt <- paste0('\nCould not delete old ', tmp,
                    ' file. New iterations will be written after old ones.')
      if (file.exists(tmp)) tryCatch(file.remove(tmp), error = function(e) cat(txt))
      rm(txt)
      if (!silent) cat("\nIteration results will be written to: \n", tmp)
    } else {
      if (!silent) cat("\nIteration results will be appended to: \n", tmp)
    }
  }

  writeItSummary <- FALSE
  if (writeItSummary) {
    itSumFname <- paste0(outputDirectory, modelName, "_itSummary.csv")
    if (!silent) cat("\nIteration summary will be written to: \n", itSumFname, "\n")
    if (!silent) cat('  (Note: This is not yet implemented...)\n')
    if (!silent) cat("\n")
  }

  # Parameters for FD gradient calculation
  hlim        <- 0.1
  negpt5      <- -0.5

  # ***  IV AND V COMPONENTS  ***
  #--------------------------------------------------------------------------
  # BGW Fortran excerpt
  # This is included for reference to components needed inside bgw_mle,
  # that must be handled differently because of R.
  # ! D....... SCALE VECTOR.
  # ! DR...... DERIVATIVES OF R AT X.
  # ! N....... TOTAL NUMBER OF RESIDUALS.
  # ! ND...... LEADING DIMENSION OF DR -- MUST BE AT LEAST PS.
  # ! NN...... LEAD DIMENSION OF R, RD.
  # ! P....... NUMBER OF PARAMETERS (COMPONENTS OF X) BEING ESTIMATED.
  # ! PS...... NUMBER OF NON-NUISANCE PARAMETERS.
  # ! R....... RESIDUALS (OR MEANS -- FUNCTIONS OF X) WHEN  DRGLG IS CALLED
  # !          WITH IV(1) = 1.
  # ! RD...... RD(I) = HALF * (G(I)**T * H(I)**-1 * G(I)) ON OUTPUT WHEN
  # !          IV(RDREQ) IS 2, 3, 5, OR 6.   DRGLG SETS IV(REGD) = 1 IF RD
  # !          IS SUCCESSFULLY COMPUTED, TO 0 IF NO ATTEMPT WAS MADE
  # !          TO COMPUTE IT, AND TO -1 IF H (THE FINITE-DIFFERENCE HESSIAN)
  # !          WAS INDEFINITE.  BEFORE CONVERGENCE, RD IS ALSO USED AS
  # !          TEMPORARY STORAGE.
  ps <- p
  nd <- ps
  nn <- n

  # The following is expanded internal documentation to
  # support implementing the passing of vectors and arrays
  # between R and BGW Fortran.
  #
  # Managing the calls to drglg. The call looks like:
  # out = bgw_drglg(d_vec, dr_vec, iv, liv, lv, n, nd,
  #                 nn, p, ps, r_vec, rd_vec, v, x,
  #                 rhoi, rhor, i_itsum)
  #
  # Key "smaller entities" (i.e, those that are not iv, v,
  # rhoi or rhor) that are being expected on the "other side"
  # are:  d_vec, dr_vec, r_vec, and rd_vec.
  #
  # In actuality, BGW is using space in v that coincides with
  # these entities. This is a Fortran-related complexity that
  # is incompatible with other languages. Moreover, inside
  # drglg sometimes these are indexed as matrices, even though
  # they are passed as vectors. For example, here are the
  # declarations:
  #!      INTEGER IV(LIV), RHOI(*)
  #!      DOUBLE PRECISION D(P), DR(ND,N), R(*), RD(*), RHOR(*),
  #!      DIMENSION RD(N, (P-PS)*(P-PS+1)/2 + 1)

  # Internal to R, we just allocate space to these:
  # d
  # scaling vector (of length p)
  d_vec    <- rep(1, p)

  # dr
  # Derivative of R wrt x (aka the Jacobian)
  len_dr   <- nd*n
  dr_vec   <- rep(0, len_dr)

  # r
  # This is nominally the n-vector vector of residuals.
  # However, it is more complicated than that.
  # First, if BGW's capability of using p-ps nuisance parameters
  # is being implemented, an additional (p-ps)*n slots are required.
  # (Recall we are not (yet) using these). So, in this case the
  # length of r would be (p-ps+1)*n (rather than just n).
  # Moreover, the BGW Fortran implementation also uses an
  # additional 2 extra vectors of length n for extra storage
  # if finite-difference gradients are being computed. In that
  # case the length of r would be be (p-ps+3)*n (as implemented
  # in BGW Fortran). We are replicating that implementation here.
  # In that case, r actually contains 3 'sub-vectors' that are
  # concatenated, and requires indexes to be properly accessed.
  # Accordingly:
  #
  if (hasAnalyticGrad) {
    len_r_all <- (p-ps+1)*n
    len_r_resid <- len_r_all
    r_resid_index <- (1:len_r_resid)
  } else {
    len_r_all <- (p-ps+3)*n
    len_r_resid <- (p-ps+1)*n
    r_resid_index <- (1:len_r_resid)
    rs1_index <- (len_r_resid+1):(len_r_resid+n)
    r21_index <- (len_r_resid+n+1):(len_r_all)
  }
  r_all_vec <- rep(0, len_r_all)

  # rd
  # Space related to regression diagnostics
  len_rd <- n*((p-ps)*(p-ps+1)/2 + 1)
  rd_vec <- rep(0, len_rd)

  # To review, we have d_vec, dr_vec, r_all_vec and rd_vec
  # for use inside this R function. However, we must keep
  # copies of these inside the v vector that gets passed
  # back and forth to drglg.

  # # For now, we assume extra storage for r, as in dglf
  # dim2_r   <- (p - ps + 3)
  # len_r    <- nn*dim2_r
  # r_vec    <- rep(0,len_r)
  # r_ind    <- 1:n
  # rs1_ind  <- (n+1):(2*n)
  # r21_ind  <- (2*n+1):(3*n)
  # #
  # dim2_rd  <- 1
  # len_rd   <- nn*dim2_rd
  # rd_vec   <- rep(0,len_rd)
  # #
  # nf <- 0  # Why did I put this here? FIX ME
  # DSB note: These could be allocating unneeded
  # storage.  FIX ME?
  rhoi <- rep(0, 2)
  rhor <- rep(1, n)

  # Pointers/indexes for v
  d1_v        <- iv[bgw_iv_index$d]
  d_index_v   <- d1_v:(d1_v+p-1)
  # d1_v needed in FD loop
  #
  dr1_v       <- iv[bgw_iv_index$j]
  dr_index_v  <- dr1_v:(dr1_v+nd*n-1)
  # dr1_v needed in FD loop
  #
  r1_v            <- iv[bgw_iv_index$r]
  # r_all_index_v   <- r1_v:(r1_v+(p - ps + 3)*n-1)
  r_all_index_v   <- r1_v:(r1_v+len_r_all-1)
  r_resid_index_v <- r1_v:(r1_v+len_r_resid-1)
  rm(r1_v)
  #
  rd1_v         <- iv[bgw_iv_index$regd0]
  rd_index_v    <- rd1_v:(rd1_v + len_rd - 1)
  if (hasAnalyticGrad) {
    rm(rd1_v)
  } else {
    r21_v       <- rd1_v - n
    r21_index_v <- r21_v:(r21_v + n - 1)
    rs1_v       <- r21_v - n
    rs1_index_v <- rs1_v:(rs1_v + n - 1)
    # rm(rd1_v,r21_v,rs1_v)
  }
  g1_v      <- iv[bgw_iv_index$g]
  g_index_v <- g1_v:(g1_v + p - 1)
  rm(g1_v)

  estimationFinished <- 0

  if (debug) write.csv(iv, "bgw_iv_r_preloop.csv")
  if (debug) write.csv(v, "bgw_v_r_preloop.csv")

  # DSB November 5, 2026
  # This changes the initiation of d in the
  # case of analytic gradients.
  # Background: in the original BGW Fortran, David Gay chose the defaults for
  # dynamic (adaptive) scaling, and with analytic derivatives D was not
  # initialized to ones. (In the finite-difference branch D must be set to ones
  # at the start, because adaptive scaling needs scale values and no derivatives
  # are available yet to compute them.) Extensive numerical testing by the
  # Apollo developers found that starting from D = 1 also works better when
  # analytic derivatives are available, so that is now the default here.
  if (hasAnalyticGrad) {
    if ((iv[bgw_iv_index$dtype] == 1) && (v[bgw_v_index$dinit] == 0.)) {
      v[bgw_v_index$dinit]   <- 1.e0
    }
  }

  # v[bgw_v_index$lmax0] <- 0.2
  ltest <- FALSE
  i_itsum <- -1
  lastIterWritten <- -1     # the last iteration written to the iteration files
  pendingLL       <- NULL   # log-likelihood of an iterate whose gradient has not yet been written
  pendingBeta     <- NULL   # and that iterate
  requested       <- NULL   # an iterate whose gradient has been requested from calcJ or by FD
  lastGood        <- NULL   # the last iterate at which the gradient could be computed
  while (estimationFinished == 0) {
    i_itsum <- i_itsum + 1
    out <- bgw_drglg(d_vec, dr_vec, iv, liv, lv, n, nd,
                     nn, p, ps, r_all_vec, rd_vec, v, x,
                     rhoi, rhor, i_itsum)
    d_vec        <- out[[1]]
    dr_vec       <- out[[2]]
    iv           <- out[[3]]
    r_all_vec    <- out[[4]]
    rd_vec       <- out[[5]]
    v            <- out[[6]]
    x            <- out[[7]]
    rhoi         <- out[[8]]
    rhor         <- out[[9]]
    i_itsum      <- out[[10]]

    if ((i_itsum == 0) && (bgw_settings[["scalingMethod"]] == "userScaling")) {
      d_vec <- bgw_settings[["userScaleVector"]]
    }

    v[d_index_v]     <- d_vec
    v[dr_index_v]    <- dr_vec
    v[r_all_index_v] <- r_all_vec
    v[rd_index_v]    <- rd_vec
    g_vec <- v[g_index_v]

    # The gradient and scale vector for an iterate exist only after the call following the one
    # that requested its derivatives, so the iterate is recorded and written here. If the gradient
    # could not be computed (iv[1] = 65), nothing is recorded or written for it.
    if (!is.null(requested)) {
      if (iv[1] != 65) lastGood <- c(requested, list(g = g_vec, d = v[d_index_v]))
      requested <- NULL
    }
    if (!is.null(pendingLL)) {
      if (iv[1] != 65) {
        if (writeIter) {
          bgw_writeIterations(pendingBeta, pendingLL, lastIterWritten, "estimate", writeNames,
                              iterFile)
        }
        if (writeIterGradient) {
          bgw_writeIterations(-g_vec, pendingLL, lastIterWritten, "gradient", writeNames, gradFile)
        }
        if (writeIterScaleVec) {
          bgw_writeIterations(v[d_index_v], pendingLL, lastIterWritten, "scaleVec", writeNames,
                              scaleFile)
        }
      }
      pendingLL   <- NULL
      pendingBeta <- NULL
    }
    if (betaIsNamed) {
      beta <- copy_XToBeta(x, beta, p)
    } else {
      beta <- x
    }

    if (ltest) {
      if (debug) write.csv(iv, "bgw_iv_r_firstcall.csv")
      if (debug) write.csv(v, "bgw_v_r_firstcall.csv")
      ltest <- FALSE
    }

    iv1 <- iv[1]
    computeProb <- 0
    computeDeriv <- 0
    if (iv1 - 2 < 0) {
      computeProb <- 1
      if (iv1 > 0) {
        computeDeriv <- 0
      } else {
        computeDeriv <- 1
      }
    }

    if (iv1 - 2 == 0) {
      computeProb <- 0
      computeDeriv <- 1
    }

    if (iv1 - 2 > 0) {
      computeProb <- 0
      # DSB marker
      # 7 November 2026: after convergence the derivatives are not passed to the
      # Fortran again, so they are not computed here. (They were computed here
      # only for indLevelGradients, which are now computed after the loop.)
      computeDeriv <- 0
      estimationFinished <- 1
    }

    if (computeProb == 1) {
      # DSB FIX ME!  Still need to trap out bad function evaluations
      nf <- iv[bgw_iv_index$nfcall] # not using this for much at the moment
      # P <- calcR(beta[1:p],logP=FALSE) New version in Apollo fixes this.
      # P <- calcR(beta[1:p])
      # P <- calcR(beta)
      P <- evalCalcR(x)
      r_all_vec[r_resid_index] <- as.numeric(P)
      v[r_resid_index_v] <- r_all_vec[r_resid_index]

      # In BGW, calcR sets nf to 0 to indicate that it was not
      # possible to compute the residual vector R (aka probabilities),
      # and hence it is not possible to compute the objective function.
      # In Apollo calcR does not return any such error codes, so
      # R values (probabilities) must be checked here.

      # Inf and negative values are failures too: the Fortran forms -log(R_i), and would otherwise
      # take -Inf as a huge improvement. (A zero gives +Inf, and the Fortran rejects that step.)
      if (any(!is.finite(P) | P < 0)) {
        if (bgw_settings$printLevel == 3L)
          cat("             Call to calcR failed (producing NA, NaN, Inf or negative probs)\n")
        nf <- 0
      }

      if ((!hasAnalyticGrad) && (nf > 0)) {
        # r_vec[rs1_ind] <- r_vec[r_ind]
        r_all_vec[rs1_index] <- r_all_vec[r_resid_index]
        v[rs1_index_v]       <- r_all_vec[r_resid_index]
      }

      # In BGW the failure is indicated by iv[toobig] = 1.
      if (nf == 0) {
        iv[bgw_iv_index$toobig] <- 1
      }
    }

    if (computeDeriv == 1) {
      if ((iv1 == 2) && (!silent))
        iv <- bgw_itsum(d_vec, g_vec, iv, v, x, p, betaIsNamed, betaNames, i_itsum)
      # if (iv1 == 2) bgw_writeIterations(beta,-v[10],modelName,outputDirectory)
      # Each iteration is written once. (With a finite-difference Hessian, the Fortran also returns
      # with iv[1] = 2 during the covariance stage, when beta is being perturbed.)
      if ((iv1 == 2) && (iv[bgw_iv_index$niter] > lastIterWritten)) {
        lastIterWritten <- iv[bgw_iv_index$niter]
        pendingLL       <- -v[bgw_v_index$f]
        pendingBeta     <- beta
      }
      # Not in the covariance stage (iv[mode] > 0), when x is being perturbed.
      if ((iv1 == 2) && (iv[bgw_iv_index$mode] <= 0)) {
        requested <- list(x = x, beta = beta, LL = -v[bgw_v_index$f],
                          iter = iv[bgw_iv_index$niter])
      }
      if (hasAnalyticGrad) {
        # 7 November 2026: calcJ may return n x p or p x n, and the derivatives of
        # log(R) (calcJ_type = "log_probability") or of R ("probability").
        # The Fortran needs the derivatives of R, p x n. calcR is evaluated only when
        # the derivatives of log(R) have to be multiplied by R.
        dr_vec <- evalCalcJ(x)
        dr_vec <- orient_calcJ(dr_vec)   # now p x n
        if (calcJ_type == "log_probability") {
          P <- evalCalcR(x)
          # column i of the p x n matrix (observation i) times R_i
          dr_vec <- dr_vec * rep(as.numeric(P), each = p)
        }
        dr_vec <- as.numeric(dr_vec)
        # A Jacobian that cannot be computed is reported to the Fortran in the same way as a failed
        # calcR or finite-difference step, by setting toobig. (The original BGW drivers asked the
        # user's calcRJ to set nf = 0 in this case.) The Fortran then stops with code 65.
        if (!all(is.finite(dr_vec))) {
          if (bgw_settings$printLevel == 3L)
            cat("             Call to calcJ failed (producing NA, NaN or Inf)\n")
          iv[bgw_iv_index$toobig] <- 1
        }
      } else {
        # % 50 in dglf
        # %  ***  COMPUTE FINITE-DIFFERENCE APPROXIMATION TO
        # %       DR = GRAD. OF R  ***
        #   %     *** INITIALIZE D IF NECESSARY ***
        #   % display('Computing FD gradient...')
        # if(iv[mode] < 0 && v[dinit] == zero) {
        if ((iv[bgw_iv_index$mode] < 0) && (v[bgw_v_index$dinit] == 0.)) {
          # DSB FIX ME?  The following could be a loose end that
          # needs to addressed.  It is scaling related.
          # % [p, v(d1), one]=dv7scp(p, v(d1), one);
          v[d_index_v] <- rep(1, p)
          d_vec <- rep(1, p)
        }
        dk <- d1_v
        ng <- iv[bgw_iv_index$ngcall] - 1
        if (iv[1] == -1) {
          # if iv(1) == -1
          # DSB Note:  This is related to covariance calculations
          # performed after convergence has been achieved.
          # iv(ngcov) = fix(iv(ngcov) - 1);
          # Need to double check this.
          iv[bgw_iv_index$ngcov] <- iv[bgw_iv_index$ngcov] - 1
        }

        # Corresponds to code at line 70 of dglf
        # DSB Nov 4, 2026: In v 0.1.4 and earlier we did not address scaling.
        # Old line of code:
        # dk <- d1_v # This is currently unused, because we don't address scaling
        # Initialize index for the scaling vector d.
        dk <- d1_v
        j1k0 <- dr1_v
        stopLoop <- 0
        nf <- iv[bgw_iv_index$nfcall]
        if (nf != iv[bgw_iv_index$nfgcal]) {
          # If these do not match, it means we need to generate a current
          # copy of r and store it in rs1.
          # print(paste(" Updating prob at iteration ",iv[31]))
          ng <- ng + 1
          # P <- calcR(beta)
          P <- evalCalcR(x)
          r_all_vec[rs1_index] <- as.numeric(P)
          v[rs1_index_v] <- r_all_vec[rs1_index]
          # Remark: nf is set to zero if probabilities cannot be computed.
          if (any(!is.finite(P) | P < 0)) { nf <- 0}
          if (nf == 0) {
            iv[bgw_iv_index$toobig] <- 1
            iv[bgw_iv_index$ngcall] <- ng
            stopLoop <- 1
            computeDeriv <- 0
          }
        }
        # The first line of this loop corresponds to line 70 in dglf.
        for (k in 1:ps) {
          if (stopLoop == 1) break
          # DSB FIX ME: The following includes references to v[dk].
          # Investigate if this is a loose end that needs to be addressed.
          # print(paste(" FD calculation for component k: ", k))
          xk <- x[k]
          # DSB Nov 4, 2026: In v 0.1.4 and earlier, we used 1 rather than v[dk]
          # when performing these calculations.
          # hk <- v[bgw_v_index$dltfdj]*max(abs(xk),1)
          hk <- v[bgw_v_index$dltfdj]*max(abs(xk), 1/v[dk])
          # print(v[dk])
          h0 <- hk
          # DSB Nov 4, 2026: In v 0.1.4 and earlier we did not bother
          # to increment dk because we were not using it. This version
          # revives this incrementing.
          #
          # go ahead and increment dk to be used in the next iteration
          dk <- dk + 1
          compute_r21_loop <- 1
          while (compute_r21_loop == 1) {
            # Loop for computing r21
            x[k] <- xk + hk
            # print(x)
            # For dglf, ng keeps track of the function evaluations
            # required for computing FD derivatives
            ng <- ng + 1
            nf <- -ng
            # DSB FIX ME!  As elsewhere, need to set nf = 0 if there is a problem.
            # beta = copy_XToBeta(x,beta,p)
            P <- evalCalcR(x)
            if (any(!is.finite(P) | P < 0)) { nf <- 0}
            r_all_vec[r21_index] <- as.numeric(P)
            v[r21_index_v] <- r_all_vec[r21_index]
            # DSB FIXME
            # We need the following functionality:  if the probability calculation fails
            # then set nf to 0.  Otherwise, do not change nf.
            # This means that nf will be either negative or zero.
            # Negative means we are done with this loop.
            # nf = 0 means we can try again with a smaller h.
            # However, if this gets too small, we must abort completely.
            # [It might be possible to improve this...]
            if (nf < 0) { #nf check
              compute_r21_loop <- 0
            } else {
              # Since there was a problem computing P, make h smaller
              hk <- negpt5 * hk
              if (abs(hk/h0) < hlim) {
                iv[bgw_iv_index$toobig] <- 1
                compute_r21_loop <- 0
                computeDeriv <- 0
                break
              }
            } # end nf check
          } # end compute r21 loop
          # Finish up calculations for k loop, if things are still good
          x[k] <- xk
          iv[bgw_iv_index$nfgcal] <- ng
          i1 <- r21_v
          j1k <- j1k0
          j1k0 <- j1k0 + 1
          for (i in rs1_index_v) {
            v[j1k] <- (v[i1] - v[i])/hk
            # print(v[j1k])
            i1 <- i1 + 1
            j1k <- j1k + ps
          } # end rs1_index loop
        } # end k loop
        dr_vec <- v[dr_index_v]
      }
    } # end compute derivative loop
  } # estimation finished

  # If the gradient could not be computed at a later iterate (code 65), the Fortran has restored the
  # objective of the previous iterate, but x is the new point and the gradient has been set to zero.
  # The last iterate at which the gradient was computed is returned instead, with its objective,
  # gradient, scale vector and iteration number, so that all of these refer to the same point.
  if ((iv[1] == 65) && !is.null(lastGood)) {
    x    <- lastGood$x
    beta <- lastGood$beta
    v[g_index_v]           <- lastGood$g
    v[d_index_v]           <- lastGood$d
    d_vec                  <- lastGood$d
    v[bgw_v_index$f]       <- -lastGood$LL
    iv[bgw_iv_index$niter] <- lastGood$iter
  }

  g_vec <- v[g_index_v]
  if (!silent)
    iv <- bgw_itsum(d_vec, g_vec, iv, v, x, p, betaIsNamed, betaNames, i_itsum)

  model <- bgw_buildModel(iv, v, beta, betaStart, p, n, bgw_settings, hasAnalyticGrad, g_index_v,
                          d_index_v, silent)

  # If the estimation fails at the start, the Fortran restores v(f) from v(f0), which has not been
  # set yet, so the objective it reports is zero. The log-likelihood at betaStart is known here.
  if ((iv[1] > 62) && is.null(lastGood)) {
    model$finalLL <- sum(log(as.numeric(Prob)))
    model$maximum <- model$finalLL
  }

  if (returnIndLevelGradients) {
    model$indLevelGradients <- bgw_indLevelGradients(x, beta, v, d_index_v, n, p, hasAnalyticGrad,
                                                     calcJ_type, evalCalcR, evalCalcJ,
                                                     orient_calcJ, silent)
  }

  # model$bgw_iv <- iv
  # model$bgw_v <- v

  # model$LLout           <- -v[10]
  # model$nIter <- iv[31]
  # model$iv <- iv
  # model$v  <- v

  class(model) <- "bgw_mle"
  # c("\n")
  # print(attributes(model))
  # c("\n")
  # cat("\n Last line of BGW_mle")
  # cat("\n")
  return(model)
}
