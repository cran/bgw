#' @title bgw_mle
#'
#' @description Performs maximum likelihood estimation (MLE) for the user-provided model
#' defined in \code{bgw_calcR}.
#'
#' @details This function has been written to provide an R-based interface to Fortran estimation software published in
#' Bunch, Gay and Welsch (1993), "Algorithm 717-Subroutines for Maximum Likelihood and Quasi-Likelihood
#' Estimation of Parameters in Nonlinear Regression Models," ACM Transactions on Mathematical Software, 19 (1),
#' March 1993, 109-130. The letters BGW will be used in various ways to denote the source of the
#' estimation functionality.
#'
#' A primary motivation was to develop a more efficient maximum likelihood estimation function for use in
#' the Apollo choice modelling package: see \url{http://www.apollochoicemodelling.com/}. However, we
#' have adopted a design whereby the BGW package is wholly independent of Apollo, and can be
#' used in a stand-alone fashion. Note also that the BGW Fortran subroutines are written to support
#' general statistical estimation for an arbitrary objective/criterion function. So, although this
#' version of the package is specifically written for MLE, the package may see future updates that
#' expand the number of estimation options (for, e.g., nonlinear least squares, generalized method of
#' moments, etc.).
#'
#' Remark: Following the convention in the numerical optimization literature, BGW minimizes the
#' objective function. That is, bgw_mle minimizes the negative-log-likelihood for the model calcR.
#'
## Note:  The following parameters will need to be updated to final versions once we can coordinate
## with Apollo.
##
## calcR       Function that computes an n-vector of model residuals (R) for a p-vector of parameters beta.
##                In this case the residuals are likelihoods (probabilities).
##                The first argument of calcR must be the parameter vector beta.
## calcJ      User-provided function that computes the Jacobian of R. If NULL, finite-difference derivatives are used.
##                In matrix form, dim=c(p,n).  However, it could be stored as vector in column-major order.
## start      initial vector of starting values for parameters.
##
## bgw_settings control parameters for BGW algorithm.
#'
#' @param calcR Function that computes an n-vector (R) of model residuals for a p-vector of (numeric) parameters beta (the first argument). In this case the residuals are likelihoods (probabilities). (The beta vector can be named or unnamed.)
#' @param betaStart Vector of initial starting values for beta. Can be either a named or unnamed vector.
#' @param calcJ Function that computes the matrix of partial derivatives of R wrt beta (a.k.a. the Jacobian). If NULL, finite-difference derivatives are used. In matrix form, dim=c(p,n).  However, it could be stored as a vector in column-major order.
#' @param bgw_settings List. Contains control parameters for BGW estimation code. All parameters have default values, so user input is entirely optional.
#'                \itemize{
#'                  \item \strong{\code{printLevel}}: Integer (0-3). Controls the level of detail for iteration information written to the console during estimation.  0 = silent, 1 = print starting values and final solution only, 2 = short line summary, 3 = long line summary (see documentation). Default = 3.
#'                  \item \strong{\code{silent}}: Logical. Suppresses all output. Default = FALSE. (Currently redundant with printLevel = 0.)
#'
#'                  \item \strong{\code{printNonDefaultSettings}}: Logical. Echos any user-provided bgw_settings. Default = TRUE.
#'                  \item \strong{\code{printStartingValues}}: Logical. Default = TRUE. Print starting values for beta.
#'                  \item \strong{\code{printFinalResults}}:  Logical. Default = TRUE. Print final status of estimation (convergence/error message), summary statistics (negative log-likelihood, iterations, function evals, etc.), beta, gradient, and (if available) estimated standard errors and t-ratios. (See vcHessianMethod for variance-covariance specification.)
#'                  \item \strong{\code{maxIterations}}: Numeric. Maximum number of iterations for the estimation. BGW default is 150.
#'                  \item \strong{\code{maxFunctionEvals}}: Numeric. Maximum number of objective function evaluations. BGW default is 200.
#'                  \item \strong{\code{modelName}}: Character. Used to create names for requested output files. Default is "bgw_mle_model."
#'                  \item \strong{\code{outputDirectory}}: Character. Name of sub-directory (~/outputDirectory) to write output files. Default is NULL (current directory).
#'                  \item \strong{\code{writeItSummary}}: Logical. If TRUE, iteration information is echoed in the output file "modelName_itSummary.csv". Default is FALSE. [Not currently implemented.]
#'                  \item \strong{\code{writeIter}}: Logical. If TRUE, parameters and log-likelihood for each iteration are written to "modelName_iterations.csv". Default is FALSE.
#'                  \item \strong{\code{vcHessianMethod}}: Character. Method for computing the Hessian approximation used for the variance-covariance matrix (VC = H^(-1)). Options are: "none","bhhh","finiteDifferences", or "fdFunction" ("finiteDifferences" automatically uses gradient differences if a gradient is available, otherwise it uses objective function differences. "fdFunction" allows the user to use objective function differences even if a gradient is available). Default is "bhhh."
#'                  \item \strong{\code{scalingMethod}}: Character. Method used to compute a scaling vector (scaleVec_i, i=1,..,p). Define a matrix D = diag(scaleVec_1,..., scaleVec_p). When using scaling, values in scaleVec should be chosen so that the elements of D*beta are roughly comparable in size. The re-scaled beta is used when computing trial steps using trust regions, and when computing stopping criteria. Options are: "adaptive," "none," and "userScaling": the default is "adaptive."  For a description of the "adaptive" method, which updates scaleVec at each iteration using information from the model Jacobian, see Bunch, Gay, and Welsch (1993). For "none," scaleVec is set to a p-vector of ones for the entire search. The "userScaling" option indicates that the user is supplying their own (fixed) scaleVec in bgw_settings[["userScaleVector"]] (see next item). Both items must be properly set or an error occurs.
#'                  \item \strong{\code{userScaleVector}}: Numeric, with dimension p = number of free parameters. Can be either a named or unnamed vector. This is a user-provided scaling vector that is used ONLY in conjunction with the non-default option "userScaling" for bgw_settings[["scalingMethod"]] (see previous item).
#'                }
#' @return model object of class 'bgw_mle'. Output of a bgw maximum likelihood estimation procedure. A list with the following attributes:
#' \itemize{
#' \item \strong{\code{betaStart}}: Vector of initial starting values.
#' \item \strong{\code{bgw_settings}}: List. The same as the input argument.
#' \item \strong{\code{hasAnalyticGrad}}: Logical. Indicates in an analytical gradient calculation was used. If the user has not provided a calcJ function (see input parameter), it is set to FALSE.
#' \item \strong{\code{numParams}}: Numeric. Number of model parameters used in calcR.
#' \item \strong{\code{numResids}}: Numeric. Number of independent observations (model residuals) in data set = dimension of calcR output.
#' \item \strong{\code{code}}: Numeric. Numeric return code from BGW.
#' \item \strong{\code{message}}: Character. Message statement characterizing termination of MLE search.
#' \item \strong{\code{betaStop}}: Vector. Value of parameter vector at conclusion of MLE search. See message to determine if beta is a valid estimate.
#' \item \strong{\code{finalLL}}: Numeric. Value of log-likelihood at betaStop.
#' \item \strong{\code{iterations}}: Numeric. Number of iterations used in MLE search. In BGW, this is the same as the number of gradient evaluations.
#' \item \strong{\code{functionEvals}}: Numeric. Number of function evaluations used in MLE search. (This excludes function evaluations used by any finite-difference calculations for the gradient and/or the vcHessian.
#' \item \strong{\code{gradient}}: Vector. The gradient evaluated at betaStop.
#' \item \strong{\code{scaleVec}}: Vector. The scaling vector at the conclusion of the MLE search. Note: In the current version, this will be a p-vector of 1's (used throughout the search). In future versions, additional scaling options may be implemented.
#' \item \strong{\code{estimate}}: Vector. MLE parameter vector obtained by BGW. The same as betaStop if a valid convergence condition is achieved. Null otherwise.
#' \item \strong{\code{maximum}}: Numeric. Final log-likelihood value for a (successful) MLE search.
#' \item \strong{\code{hessianMethodAttempted}}: String. Requested method for computing vcHessian (from bgw_settings).
#' \item \strong{\code{hessianMethodUsed}}: String. Method actually used for computing vcHessian (if vcHessian was requested, and if the computation was successful. If not, a message indicating 'no request' or 'singular vcHessian' is provided.)
#' \item \strong{\code{vcHessianConditionNumber}}: Numeric. Estimated upper bound on reciprocal of Euclidean condition number of vcHessian (if available). Set to -1 if unavailable.
#' \item \strong{\code{varcovBGW}}: Matrix. p-by-p matrix containing estimate of the variance-covariance matrix (if requested and available).
#' \item \strong{\code{vcVec}}: Vector. Lower triangle of variance-covariance matrix stored in vector form (row-major order, if requested and available).
#' \item \strong{\code{seBGW}}: Vector. Estimated standard errors for parameter estimates (if requested and available).
#' \item \strong{\code{tstatBGW}}: Vector. Estimated t-statistics (versus 0, if requested and available).
#' }
#' @export
bgw_mle <- function(calcR, betaStart, calcJ=NULL, bgw_settings=NULL) {
  ### Parameters:
  ## calcR       function that computes an n-vector of model residuals for a p-vector of parameters x.
  ##                In this case the residuals are likelihoods (probabilities).
  ##                The first argument must be the parameter vector x.
  ## calcJ      derivative of R wrt x. If NULL, finite-difference derivatives are used.
  ##                In matrix form, dim=c(p,n).  However, it could be stored as vector in
  ##                  column-major order.
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
  # scalingMethod	[Currently fixed to default d=1]
  #
  # Possible convergence tolerance parameters
  # RelativeFunction_ConvergeTol
  # Step_ConvergeTol
  # Singular_ConvergeTol
  # False_ConvergeTol
  # AbsoluteFunction_ConvergeTol
  # ##------------------------------------------------------------------------
  # ## This perhaps should be moved to a function somewhere
  # ## Convergence codes and messages
  bgw_message_table = list()
  bgw_message_table[[1]]  = "       ***** Evaluate function ***** \n"
  bgw_message_table[[2]]  = "       ***** Evaluate gradient ***** \n"
  bgw_message_table[[3]]  = "       ***** X-convergence ***** \n"
  bgw_message_table[[4]]  = "       ***** Relative function convergence ***** \n"
  bgw_message_table[[5]]  = "       ***** X- and relative function convergence ***** \n"
  bgw_message_table[[6]]  = "       ***** Absolute function convergence ***** \n"
  bgw_message_table[[7]]  = "       ***** Singular convergence ***** \n"
  bgw_message_table[[8]]  = "       ***** False convergence ***** \n"
  bgw_message_table[[9]]  = "       ***** Function evaluation limit ***** \n"
  bgw_message_table[[10]] = "       ***** Iteration limit ***** \n"
  bgw_message_table[[11]] = "       ***** stopx ***** \n"
  # Note:  The following occurs on entry to ditsum if iv[1] > 62
  # If iv[1] > 61, iv1 = iv[1] - 51
  # iv[1] = 63 => iv1 = 12
  bgw_message_table[[12]] = "       ***** Initial f(x) cannot be computed ***** \n"
  # iv[1] = 64 => iv1 = 13
  bgw_message_table[[13]] = "       ***** Bad parameters to assess ***** \n"
  # iv[1] = 65 => iv1 = 14
  bgw_message_table[[14]] = "       ***** Gradient could not be computed ***** \n"
  # iv[1] = 66 => iv1 = 15
  bgw_message_table[[15]] = "       ***** Inconsistent dimensions (Could be N, ND, P, or LIV. )***** \n"
  # iv1 >= 16 (i.e., iv[1] >= 67) simply writes out iv[1]
  # These values would typically occur when the user overwrites
  # a default value with an illegal value. We will fill in the details.
  # iv[1] = 67 => iv1 = 16
  # This can be set two different ways in dparck. The following should cover it.
  bgw_message_table[[16]] = " ***** Internal error specifying alg in divset and/or dparck. It should be 1. ***** \n"

  bgw_success = list()
  bgw_success[[1]] = 0
  bgw_success[[2]] = 0
  bgw_success[[3]] = 0
  bgw_success[[4]] = 0
  bgw_success[[5]] = 0
  bgw_success[[6]] = 0
  bgw_success[[7]] = 1
  bgw_success[[8]] = 1
  bgw_success[[9]] = 1
  bgw_success[[10]] = 1
  bgw_success[[11]] = 1
  bgw_success[[12]] = 1
  bgw_success[[13]] = 1
  bgw_success[[14]] = 1
  bgw_success[[15]] = 1
  # In DITSUM IV(1) < 2 or > 15	Write out value of IV(1)
  ##------------------------------------------------------------------------
  # #################### #
  #### Loading Inputs ####
  # #################### #

  # Initial test for betaStart and calcR
  test <- is.numeric(betaStart) && is.vector(betaStart) && is.function(calcR)
  if(!test) stop('Arguments calcR (function) and betaStart (numeric vector) must be provided to bge_mle.')

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
  if (any(is.na(Prob)))  stop(" Initial call to calcR failed, producing NaNs.")
  if (!is.numeric(Prob)) stop(" Initial call to calcR returned a non-numeric result.")
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
      x <- rep(0,p)
      for (i in 1:p) {x[i] <- beta[[i]]}
      return(x)
    }
    copy_XToBeta <- function(x, beta, p) {
      for (i in 1:p) {beta[[i]] <- x[i]}
      return(beta)
    }
    x = copy_BetaToX(beta,p)
  }

  # This needs to be error checked next.
  # But, echoing valid information to user is done later.
  if (!is.null(calcJ)) {
    if (is.function(calcJ)) {
      hasAnalyticGrad = TRUE
      # cat("BGW using analytic model derivatives supplied by caller...\n")
    } else {
      stop("\nArgument calcJ provided to bgw_mle must be a function.\n")
    }
  } else {
    hasAnalyticGrad = FALSE
    # cat("BGW is using FD derivatives for model Jacobian. (Caller did not provide derivatives.)\n")
  }

  # bgw_setup allocates space for iv and v, and sets up default values
  bgw_setup = bgw_mle_setup(p, n, hasAnalyticGrad)
  iv  <- bgw_setup$iv_r
  liv <- length(iv)
  v   <- bgw_setup$v_r
  lv  <- length(v)
  debug <- FALSE
  if (debug) write.csv(iv,"bgw_iv_r_setup.csv")
  if (debug) write.csv(v,"bgw_v_r_setup.csv")

  # First, default bgw_settings (for those values that can be altered by the user).
    bgw_settings_default        <- list()
    bgw_settings_type           <- list()
    bgw_settings_validDiscrete  <-list()
    bgw_settings_continuousLB   <- list()
    bgw_settings_continuousUB   <- list()
  #
  # Printing/writing options
  # Fortran BGW has a relatively complex subroutine DITSUM that produces output.
  # A variety of control parameters in iv are available to select options.
  # [Depending on the options, DITSUM writes out the initial beta and scaling vector,
  # an iteration summary line for each iteration, final status, and final beta
  # vector and gradient (which may or may not correspond to a converged solution).
  #
  # Because R does not allow any Fortran I/O statements, this routine had to be
  # replaced with an R function (bgw_itsum.R), which writes all DITSUM output
  # to the R console. One possible future option would be to also save the
  # same output to a file (see
  #
  # DITSUM standard output contains two iteration summary formats: a long line
  # (the default) and a short line. Patterned after Apollo, the current
  # standard output (long line, the most verbose) is assigned printLevel = 3L.
  # We hae also implemented the short line (printLevel = 2L.) The option
  # printLevel = 1L is also assigned to the short line, but in the future could
  # be modified to provide an even less verbose option. The option printLevel = 0L
  # suppresses all DITSUM-related output.
    outlev_iv <- 19
    iv[outlev_iv] <- 1
    # We have some redundancy. "silent" is the same as "printLevel = 0L".
    silent <- FALSE
    bgw_settings_default[["silent"]]             <- silent
    bgw_settings_type[["silent"]]                <- "discrete"
    bgw_settings_validDiscrete[["silent"]]       <- c(FALSE,TRUE)
    bgw_settings_continuousLB[["silent"]]        <- 0
    bgw_settings_continuousUB[["silent"]]        <- 0
    #
    bgw_settings_default[["printLevel"]]         <- 3L
    bgw_settings_type[["printLevel"]]            <- "discrete"
    bgw_settings_validDiscrete[["printLevel"]]   <- c(0L,1L,2L,3L)
    bgw_settings_continuousLB[["printLevel"]]    <- 0
    bgw_settings_continuousUB[["printLevel"]]    <- 0
  #
    printNDS <- TRUE
    bgw_settings_default[["printNonDefaultSettings"]]         <- printNDS
    bgw_settings_type[["printNonDefaultSettings"]]            <- "discrete"
    bgw_settings_validDiscrete[["printNonDefaultSettings"]]   <- c(FALSE,TRUE)
    bgw_settings_continuousLB[["printNonDefaultSettings"]]    <- 0
    bgw_settings_continuousUB[["printNonDefaultSettings"]]    <- 0
  #
    bgw_settings_default[["printStartingValues"]]         <- TRUE
    bgw_settings_type[["printStartingValues"]]            <- "discrete"
    bgw_settings_validDiscrete[["printStartingValues"]]   <- c(FALSE,TRUE)
    bgw_settings_continuousLB[["printStartingValues"]]    <- 0
    bgw_settings_continuousUB[["printStartingValues"]]    <- 0
  #
    bgw_settings_default[["printFinalResults"]]         <- TRUE
    bgw_settings_type[["printFinalResults"]]            <- "discrete"
    bgw_settings_validDiscrete[["printFinalResults"]]   <- c(FALSE,TRUE)
    bgw_settings_continuousLB[["printFinalResults"]]    <- 0
    bgw_settings_continuousUB[["printFinalResults"]]    <- 0
  # Note that BGW uses multiple flags to control different aspects of this.
  # Some of these could be added to bgw_itsum at a later time. Either way,
  # additional documentation to explain the iteration summary will be required.
  # (DITSUM output can also be written to a file. See below.)
  #
  # Stopping tests
  # maximum number of function evaluations. BGW (Fortran) default is 200.
  # Re-defining default for BGW R package to 300.
    mxfcal_iv <- 17
    bgw_settings_default[["maxFunctionEvals"]]        <- iv[mxfcal_iv]
    bgw_settings_default[["maxFunctionEvals"]]        <- 300
    iv[mxfcal_iv]                                     <- bgw_settings_default[["maxFunctionEvals"]]
    bgw_settings_type[["maxFunctionEvals"]]           <- "continuous"
    bgw_settings_validDiscrete[["maxFunctionEvals"]]  <- list()
    bgw_settings_continuousLB[["maxFunctionEvals"]]   <- 1
    bgw_settings_continuousUB[["maxFunctionEvals"]]   <- Inf
  # maximum number of iterations. BGW (Fortran) default is 150.
    mxiter_iv <- 18
    bgw_settings_default[["maxIterations"]]        <- iv[mxiter_iv]
    # Re-defining default for BGW R package to 200
    bgw_settings_default[["maxIterations"]]        <- 200
    iv[mxiter_iv]                                  <- bgw_settings_default[["maxIterations"]]
    bgw_settings_type[["maxIterations"]]           <- "continuous"
    bgw_settings_validDiscrete[["maxIterations"]]  <- list()
    bgw_settings_continuousLB[["maxIterations"]]   <- 1
    bgw_settings_continuousUB[["maxIterations"]]   <- Inf
  #
  # To write output to files, we need filenames. Output filenames are
  # constructed from a modelName. To ensure there
  # is always a modelName, we adopt a default.
    bgw_settings_default[["modelName"]]        <- "bgw_mle_model"
    bgw_settings_type[["modelName"]]           <- "NA"
    bgw_settings_validDiscrete[["modelName"]]  <- list()
    bgw_settings_continuousLB[["modelName"]]   <- 0
    bgw_settings_continuousUB[["modelName"]]   <- 0
  #
  # We also require an output subdirectory. Default is to use current directory.
    bgw_settings_default[["outputDirectory"]]        <- ''
    bgw_settings_type[["outputDirectory"]]           <- "NA"
    bgw_settings_validDiscrete[["outputDirectory"]]  <- list()
    bgw_settings_continuousLB[["outputDirectory"]]   <- 0
    bgw_settings_continuousUB[["outputDirectory"]]   <- 0
  #
  # In Fortran BGW, DITSUM output is written to Fortran unit iv[prunit].
  # The current implementation in R only writes DITSUM output to the console.
  # The following is a placeholder (not yet implemented) to write the
  # same output to a file.
    bgw_settings_default[["writeItSummary"]] <- FALSE
    bgw_settings_type[["writeItSummary"]]    <- "discrete"
    # bgw_settings_valid[["writeItSummary"]]   <- c(FALSE,TRUE)
    # Right now the only valid option is FALSE.
    # The TRUE option is yet to be implemented.
    bgw_settings_validDiscrete[["writeItSummary"]]  <- c(FALSE)
    bgw_settings_continuousLB[["writeItSummary"]]   <- 0
    bgw_settings_continuousUB[["writeItSummary"]]   <- 0
  # If the option is true, DITSUM output goes to
  # "outputDirectory/modelName_itSummary.csv".
  #
  # Another potentially important output is a more detailed record of iterates.
  # The writeIter option writes details of each iteration to a file
  # (specifically, the current beta vector and negative log-likelihood).
    bgw_settings_type[["writeIter"]]           <- "discrete"
    bgw_settings_validDiscrete[["writeIter"]]  <- c(FALSE,TRUE)
    bgw_settings_default[["writeIter"]]        <- FALSE
    # If true, the file is "modelName_iterations.csv."
    bgw_settings_continuousLB[["writeIter"]]   <- 0
    bgw_settings_continuousUB[["writeIter"]]   <- 0

  # Additional option related to writeIter:  overWrite or append?
    bgw_settings_type[["writeIterMode"]]           <- "discrete"
    bgw_settings_validDiscrete[["writeIterMode"]]  <- c("overWrite","append")
    bgw_settings_default[["writeIterMode"]]        <- "overWrite"
    bgw_settings_continuousLB[["writeIterMode"]]   <- 0
    bgw_settings_continuousUB[["writeIterMode"]]   <- 0

  # Variance-covariance (vc) matrix calculation
  # First, indicate whether or not to compute a vc matrix.
  # Note that BGW has capabilities for both vc matrices and more advanced
  # regression diagnostics. We only allow for vc matrices in this package.
  # This is controlled by iv(rdreq_iv) = iv(57)
    rdreq_iv <- 57
  # iv(rdreq_iv) = 0 means do not compute anything.
  #              = 1 means compute VC matrix.
  # BGW default is 3 because  of the regression diagnostic
  # options.  We immediately override so that it is 1.
    iv[rdreq_iv] <- 1
  # This means that we will compute a VC matrix, using
  # the default method (BHHH) discussed next.
  # bgw_settings_default[["vcHessianMethod"]] <= "bhhh" (see below).
  #
  # Assuming a vc matrix is to be calculated, BGW needs to
  # know the type of Hessian approximation to use.
  # This is controlled by iv[covreq_iv] = iv[15].
  # The options are:
  # iv[covreq_iv] =  0, 1, 2 => a FD Hessian from gradient differences.
  #               = -1, -2   => a FD Hessian from function differences,
  #               = -3 or 3  => a BHHH Hessian (a.k.a. Gauss-Newton Hessian)
  #
  # Currently, if |iv[covreq_iv]| = 1 or 2 then VC = H^(-1);
  #                               = 3      then VC = BHHH^(-1)
  # [An alternate version of the code computes robust VC's, but this is
  # not available in the 1993 BGW]
  # The BGW default value is 1, but if an analytic gradient is not available
  # it must be changed to -1. We are overriding this and changing it to 3.
  # However, note that the iv[rdreq_iv] precludes calculating it.
    covreq_iv <- 15
    iv[covreq_iv] <- 3
    bgw_settings_default[["vcHessianMethod"]]       <- "bhhh"
    bgw_settings_type[["vcHessianMethod"]]          <- "discrete"
    bgw_settings_validDiscrete[["vcHessianMethod"]] <- c("none","bhhh",
                                                         "finiteDifferences","fdFunction")
    bgw_settings_continuousLB[["vcHessianMethod"]]  <- 0
    bgw_settings_continuousUB[["vcHessianMethod"]]  <- 0
  # Alternative settings that can be used are:
  #  bgw_settings[["vcHessianMethod"]] <= "none"
  #  bgw_settings[["vcHessianMethod"]] <= "bhhh"
  #  bgw_settings[["vcHessianMethod"]] <= "finiteDifferences"
  #  bgw_settings[["vcHessianMethod"]] <= "fdFunction"
  # Note that "fdGradient" is only an option if analytic model derivatives
  # are available.
  # Apollo will typically set this to "none" because it does its own
  # VC calculations.
  #
  # Scaling options.
  # BGW supports multiple scaling options.
  # The simplest two are: dynamic scaling, and no scaling (with scale vector d=1).
  # The type of scaling is controlled by iv(dtype_iv) = iv(16).
  # The options are:
  # iv[dtype_iv] = 0 => do not update d (see dinit discussion below).
  #              = 1 => update d on every iteration
  #              = 2 => update d on first iteration only.
  # BGW default is iv[dtype_iv] = 1, which indicates adaptive scaling.
  #
  # Remark: Our experience is that many problems are reasonably well
  # scaled, and setting the scaling vector d = 1 often works well.
  # The vector d can be set to 1 by setting
  #     v[dinit_v] = v[38] = 1 and iv[dtype_iv] = 0 (see above).
  # This is our preferred default.
    # dtype_iv     <- 16
    # iv[dtype_iv] <- 0
    # dinit_v      <- 38
    # v[dinit_v]   <- 1.e0
    # iv[dtype_iv] <- 1
    dtype_iv     <- 16
    dinit_v      <- 38
  #
    bgw_settings_default[["scalingMethod"]] <- "adaptive"
    # bgw_settings_default[["scalingMethod"]] <- "Fixed scaling p-vector d = 1."
    bgw_settings_type[["scalingMethod"]]           <- "discrete"
    bgw_settings_validDiscrete[["scalingMethod"]]  <- c("adaptive","none","userScaling")
    bgw_settings_continuousLB[["scalingMethod"]]   <- 0
    bgw_settings_continuousUB[["scalingMethod"]]   <- 0
    #
    bgw_settings_default[["userScaleVector"]]        <- "NA"
    bgw_settings_type[["userScaleVector"]]           <- "NA"
    bgw_settings_validDiscrete[["userScaleVector"]]  <- list()
    bgw_settings_continuousLB[["userScaleVector"]]   <- 0
    bgw_settings_continuousUB[["userScaleVector"]]   <- 0
  #
  # Stopping tolerances.
  # The next parameters involve stopping tolerances for convergence criteria.
  # For now, we assume these cannot be changed by the user, and default values
  # will be used. However, we will go ahead
  #
  # absolute function convergence tolerance.  Default is max{10e-20,machep^2}
  # Remark: This is most relevant for the zero-residual case in nonlinear least squares.
  # But note: If for some reason a model can fit choice modeling data perfectly in
  # MLE, the results would be f(x) = sum_i(ln(prob_i)) = 0.
    # afctol_v <- 31;
    # bgw_settings_default["absFunctionConvTol"] <- v[afctol_v]
  # [Add more here...]
  #
  #
  # Other options could be added, but it is possible to write other functions
  # for displaying/summarizing results based on the contents of the model object
  # AFTER the estimation has been completed.
  # -------------------------------------------------------------------------------------

  # Establishing final bgw_settings
  # If no bgw_settings are in argument list (or if a passed bgw_settings is NULL)
  #    Use default values for bgw_settings
  # Otherwise, bgw_settings are user-provided and must be checked.
  if ( (length(bgw_settings)==0) || (length(bgw_settings)==1 && is.na(bgw_settings)) ||
    (is.null(bgw_settings)) ) {
      bgw_settings <- bgw_settings_default
      if (!is.null(calcJ)) {
        if (is.function(calcJ)) {
          # hasAnalyticGrad = TRUE
          cat("\nBGW using analytic model derivatives supplied by caller...\n")
        }
      } else {
        # hasAnalyticGrad = FALSE
        cat("\nBGW is using FD derivatives for model Jacobian. (Caller did not provide derivatives.)\n")
      }
      cat("\nThere are no user-provided BGW settings. \n")
      cat("BGW settings have been set to default values...\n")
      check_bgw_settings <- FALSE
  } else {
    check_bgw_settings <- TRUE
  }

  # User has provided bgw_settings that must be checked.
  if (check_bgw_settings) {
    # We need to check for "silent" first, before proceeding.
    if ( !is.null(bgw_settings$printLevel) ) {
      if (bgw_settings$printLevel == 0L) {
        bgw_settings$silent <- TRUE
        silent <- TRUE
      }
    }
    if ( !is.null(bgw_settings$printNonDefaultSettings) ) {
      if (bgw_settings$printNonDefaultSettings == FALSE) {
        printNDS <- FALSE
      }
    }
    if ( !is.null(bgw_settings$silent) ) {
      if (bgw_settings$silent) {
        bgw_settings$printLevel <- 0L
        silent <- TRUE
        bgw_settings$printNonDefaultSettings <- FALSE
        printNDS <- FALSE
        bgw_settings$printStartingValues <- FALSE
        bgw_settings$printFinalResults <- FALSE
      }
    }
    if (!is.null(calcJ)) {
      if (is.function(calcJ)) {
        # hasAnalyticGrad = TRUE
        if (!silent) cat("\nBGW using analytic model derivatives supplied by caller...\n")
      }
    } else {
      # hasAnalyticGrad = FALSE
      if (!silent) cat("\nBGW is using FD derivatives for model Jacobian. (Caller did not provide derivatives.)\n")
    }

    userSettings          <- bgw_settings
    userSettingNames      <- names(bgw_settings)
    numUserSettingNames   <- length(userSettingNames)
    # A list of validUserSettingNames is created, in the same order as they would be in bgw_settings_default
    validUserSettingNames <- names(bgw_settings_default)[names(bgw_settings_default) %in% userSettingNames]
    numValidUserSettingNames <- length(validUserSettingNames)
    # print(validUserSettingNames)
    # print(numValidUserSettingNames)

    # A completely new version of bgw_settings is created.
    # It is initialized with bgw_settings_default, and then systematically updated
    # with valid user-provided settings.
    # Any valid user-provided settings that are different from default settings are
    # stored in modifiedSettings.  modifiedSettings require changes to BGW default values.
    bgw_settings <- bgw_settings_default
    if (numValidUserSettingNames > 0) {
      if (printNDS) cat("\nValid user-provided bgw_setting names(s) detected...\n")
      modifiedSettings <- list()
      for (i in validUserSettingNames) {
        if ( !setequal(bgw_settings[[i]],userSettings[[i]]) ) {
            if (bgw_settings_type[[i]] != "NA") {
            settingOk <- bgw_checkSetting(userSettings[[i]],
                               bgw_settings_type[[i]],bgw_settings_validDiscrete[[i]],
                               bgw_settings_continuousLB[[i]],bgw_settings_continuousUB[[i]])
            } else {
              settingOk <- TRUE
            }
            if (settingOk) {
               bgw_settings[[i]] <- userSettings[[i]]
               if (is.character(bgw_settings[[i]])) {
                  if (printNDS) cat('\nbgw_settings[[\"',i,'\"]] is set to user-provided value (\"',userSettings[[i]],'\").',sep='')
               } else {
                  if (printNDS) cat('\nbgw_settings[[\"',i,'\"]] is set to user-provided value (',userSettings[[i]],').',sep='')
               }
              modifiedSettings <- c(modifiedSettings,i)
             } else {
                if (is.character(bgw_settings[[i]])) {
                  if (printNDS) {
                    cat('\nUser-provided value for bgw_settings[[\"',i,'\"]]  (\"',userSettings[[i]],'\") is not valid.\n',sep='')
                    cat('    Default value (\"',bgw_settings[[i]],'\") will be used instead.\n',sep='')
                  } else {
                    stop('\nUser-provided value for bgw_settings[[\"',i,'\"]]  (\"',userSettings[[i]],'\") is not valid.\n',sep='')
                  }
                } else {
                  if (printNDS) {
                    cat('\nUser-provided value for bgw_settings[[\"',i,'\"]] is not valid.\n')
                    cat('    Default value (',userSettings[[i]],') will be used instead.',sep='')
                  } else {
                    stop('\nUser-provided value for bgw_settings[[\"',i,'\"]] is not valid.\n')
                  }
                }
              }
          } else {
              if (is.character(bgw_settings[[i]])) {
                if (printNDS) cat('\nUser-provided value for bgw_settings[[\"',i,']]\" is the same as default (\"',userSettings[[i]],'\").',sep='')
              } else {
                if (printNDS) cat('\nUser-provided value for bgw_settings[[\"',i,']]\" is the same as default (',userSettings[[i]],').',sep='')
              }
           }
        }
      if (!silent || printNDS) cat("\n")
      }

      if ( (numUserSettingNames > 0)&&(numUserSettingNames != numValidUserSettingNames) ) {
        badUserSettingNames <- userSettingNames[!(userSettingNames %in% validUserSettingNames)]
        for (i in badUserSettingNames) {
          if (printNDS) {
            cat('\nNon-valid user-provided bgw_setting setting name detected: \"',i,'\".\n',sep="")
          } else {
            stop('\nNon-valid user-provided bgw_setting setting name detected: \"',i,'\".\n',sep="")
          }
        }
        rm(badUserSettingNames)
        if (!silent || printNDS) cat("\n")
      }
    rm(userSettings,userSettingNames,numUserSettingNames,validUserSettingNames,numValidUserSettingNames)

    lms <- length(modifiedSettings)
    if (lms > 0) {
      # First, handle more complex situations.
      # Attempts to use a user-provided scale vector have potential issues.
      if (any(modifiedSettings=="scalingMethod") ) {
        if ( !any(modifiedSettings=="userScaleVector") ) {
          stop("Error. User has set bgw_settings to userScaling, but no userScaleVector has been provided.")
        }
      }
      if (any(modifiedSettings=="userScaleVector") ) {
        userScaleVector <- bgw_settings[["userScaleVector"]]
        if (!is.numeric(userScaleVector)) stop("Error in userScaleVector. Must be a numeric vector.")
        if (length(userScaleVector)!=p) stop("Error in userScaleVector. It is the wrong length.")
        if (any(userScaleVector <= 0)) stop("Error in userScaleVector. One or more values is non-positive.")
        if (bgw_settings[["scalingMethod"]] != "userScaling" ) {
            stop("\nUser-provided scale vector detected, but user has not set scalingMethod to userScaling \n")
          # bgw_settings[["scalingMethod"]] <- "userScaling"
          # cat("\nUser-provided scale vector detected. bgw_settings updated to reflect scalingMethod = userScaling\n")
        } else {
        # bgw_settings[["scalingMethod"]] <- "userScaling"
        iv[dtype_iv] <- 0
        v[dinit_v]   <- -1.e0
        }
      }
      #
      for (i in 1:lms) {
        if (modifiedSettings[[i]] == "scalingMethod") {
          if (bgw_settings[["scalingMethod"]] == "none") {
            iv[dtype_iv] <- 0
            v[dinit_v]   <- 1.e0
          }
        }
        if (modifiedSettings[[i]] == "printStartingValues") {
          if (bgw_settings[["printStartingValues"]] == FALSE) {
            x0prt_iv <- 24
            iv[x0prt_iv] <- 0
          }
        }
        if (modifiedSettings[[i]] == "printFinalResults") {
          if (bgw_settings$printFinalResults == FALSE) {
            solprt_iv <- 22
            iv[solprt_iv] <- 0
            statpr_iv <- 23
            iv[statpr_iv] <- -1
          }
        }
        if (modifiedSettings[[i]] == "vcHessianMethod") {
          if (bgw_settings[["vcHessianMethod"]] == "none") {
            iv[rdreq_iv] <- 0
          }
          if (bgw_settings[["vcHessianMethod"]] == "finiteDifferences") {
            if (hasAnalyticGrad) {
              iv[covreq_iv] <- 2
            } else {
              iv[covreq_iv] <- -2
            }
          }
          if (bgw_settings[["vcHessianMethod"]] == "fdFunction") {
            iv[covreq_iv] <- -2
          }
        }
        if (modifiedSettings[[i]] == "maxIterations") {
          iv[mxiter_iv] <- bgw_settings[["maxIterations"]]
        }
        if (modifiedSettings[[i]] == "maxFunctionEvals") {
          iv[mxfcal_iv] <- bgw_settings[["maxFunctionEvals"]]
        }
        # Some of the following is redundant
        if (modifiedSettings[[i]] == "printLevel") {
          if (bgw_settings[["printLevel"]] == 0L) {
            silent <- TRUE
            iv[outlev_iv] <- 0
          }
          if (bgw_settings[["printLevel"]] == 1L) {
            silent <- FALSE
            iv[outlev_iv] <- 0
          }
          if (bgw_settings[["printLevel"]] == 2L) {
            silent <- FALSE
            iv[outlev_iv] <- -1
          }
          if (bgw_settings[["printLevel"]] == 3L) {
            silent <- FALSE
            iv[outlev_iv] <- 1
          }
        }
      }
    }
  }


  test <- is.null(bgw_settings$outputDirectory) || bgw_settings$outputDirectory==''
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
  if (tmp!="/") bgw_settings$outputDirectory <- paste0(bgw_settings$outputDirectory,"/")

  modelName       <- bgw_settings$modelName
  outputDirectory <- bgw_settings$outputDirectory
  writeIter       <- bgw_settings$writeIter
  writeIterMode   <- bgw_settings$writeIterMode
  writeItSummary  <- bgw_settings$writeItSummary

  if (writeIter) {
    # Remove modelName_iterations if it exists
    iterFile <- paste0(outputDirectory,modelName,"_iterations.csv")
    tmp <- iterFile
    if (writeIterMode == "overWrite") {
      txt <- paste0('\nCould not delete old ', tmp, ' file. New iterations will be written after old ones.')
      if(file.exists(tmp)) tryCatch(file.remove(tmp), error=function(e) cat(txt))
      rm(txt)
      if (!silent) cat("\nIterates will be written to: \n",iterFile)
    }
    if (writeIterMode == "append") {
      if (!silent) cat("\nIterates will be appended to: \n",iterFile)
    }
  }

  writeItSummary <- FALSE
  if (writeItSummary) {
    itSumFname <- paste0(outputDirectory,modelName,"_itSummary.csv")
    if (!silent) cat("\nIteration summary will be written to: \n",itSumFname,"\n")
    if (!silent) cat('  (Note: This is not yet implemented...)\n')
    if (!silent) cat("\n")
  }

  # Parameters for FD gradient calculation
  hlim        <- 0.1
  negpt5      <- -0.5

  # ***  IV AND V COMPONENTS  ***
  lmax0_v     <- 35
  mode_iv     <- 35
  nfcall_iv   <- 6
  nfgcal_iv   <- 7
  ngcall_iv   <- 30
  ngcov_iv    <- 53
  toobig_iv   <- 2
  dltfdj_v    <- 43
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
  d_vec    <- rep(1,p)

  # dr
  # Derivative of R wrt x (aka the Jacobian)
  len_dr   <- nd*n
  dr_vec   <- rep(0,len_dr)

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
  r_all_vec <- rep(0,len_r_all)

  # rd
  # Space related to regression diagnostics
  len_rd <- n*((p-ps)*(p-ps+1)/2 + 1)
  rd_vec <- rep(0,len_rd)

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
  rhoi = rep(0,2)
  rhor = rep(1,n)

  # Pointers/indexes for v
  d_iv        <- 27
  d1_v        <- iv[d_iv]
  d_index_v   <- d1_v:(d1_v+p-1)
  # d1_v needed in FD loop
  rm(d_iv)
  #
  j_iv        <- 70
  dr1_v       <- iv[j_iv]
  dr_index_v  <- dr1_v:(dr1_v+nd*n-1)
  # dr1_v needed in FD loop
  rm(j_iv)
  #
  r_iv <- 61
  r1_v            <- iv[r_iv]
  # r_all_index_v   <- r1_v:(r1_v+(p - ps + 3)*n-1)
  r_all_index_v   <- r1_v:(r1_v+len_r_all-1)
  r_resid_index_v <- r1_v:(r1_v+len_r_resid-1)
  rm(r_iv,r1_v)
  #
  regd0_iv      <- 82
  rd1_v         <- iv[regd0_iv]
  rm(regd0_iv)
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
  g_iv      <- 28
  g1_v      <- iv[g_iv]
  g_index_v <- g1_v:(g1_v + p - 1)
  rm(g_iv,g1_v)

  estimationFinished <- 0

  if (debug) write.csv(iv,"bgw_iv_r_preloop.csv")
  if (debug) write.csv(v,"bgw_v_r_preloop.csv")

  # v[lmax0_v] <- 0.2
  ltest <- FALSE
  i_itsum <- -1
  while (estimationFinished == 0) {
    i_itsum <- i_itsum + 1
    out = bgw_drglg(d_vec, dr_vec, iv, liv, lv, n, nd,
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

    if ((i_itsum == 0) && (bgw_settings[["scalingMethod"]] == "userScaling") ){
      d_vec <- userScaleVector
    }

    v[d_index_v]     <- d_vec
    v[dr_index_v]    <- dr_vec
    v[r_all_index_v] <- r_all_vec
    v[rd_index_v]    <- rd_vec
    g_vec <- v[g_index_v]
    if (betaIsNamed) {
      beta = copy_XToBeta(x,beta,p)
    } else {
      beta <- x
    }

    if (ltest) {
      if (debug) write.csv(iv,"bgw_iv_r_firstcall.csv")
      if (debug) write.csv(v,"bgw_v_r_firstcall.csv")
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
      computeDeriv <- 0
      estimationFinished <- 1
    }

    if (computeProb == 1) {
      # DSB FIX ME!  Still need to trap out bad function evaluations
      nf <- iv[nfcall_iv] # not using this for much at the moment
      # P <- calcR(beta[1:p],logP=FALSE) New version in Apollo fixes this.
      # P <- calcR(beta[1:p])
      # P <- calcR(beta)
      if (betaIsNamed) {
         P <- calcR(beta)
      } else {
        P <- calcR(x)
      }
      r_all_vec[r_resid_index] <- as.numeric(P)
      v[r_resid_index_v] <- r_all_vec[r_resid_index]

      # In BGW, calcR sets nf to 0 to indicate that it was not
      # possible to compute the residual vector R (aka probabilities),
      # and hence it is not possible to compute the objective function.
      # In Apollo calcR does not return any such error codes, so
      # R values (probabilities) must be checked here.

      if (any(is.na(P))) {
        if (bgw_settings$printLevel == 3L)
          cat("             Call to calcR failed (producing NA probs)\n")
        nf <- 0
        }

      if ( (!hasAnalyticGrad)&&(nf > 0) ) {
        # r_vec[rs1_ind] <- r_vec[r_ind]
        r_all_vec[rs1_index] <- r_all_vec[r_resid_index]
        v[rs1_index_v]       <- r_all_vec[r_resid_index]
      }

      # In BGW the failure is indicated by iv[toobig_iv] = 1.
      if (nf == 0) {
        iv[toobig_iv] <- 1
      }
    }

    if (computeDeriv == 1) {
      if ( (iv1 == 2)&&(!silent) )
         iv <- bgw_itsum(d_vec, g_vec, iv, v, x, p, betaIsNamed, betaNames, i_itsum)
      # if (iv1 == 2) bgw_writeIterations(beta,-v[10],modelName,outputDirectory)
      if ( (iv1 == 2) && writeIter ) bgw_writeIterations(beta,-v[10],iterFile)
      if (hasAnalyticGrad) {
        # DSB note: It is actually possible to add code that saves
        # prior values of R to avoid recomputing it (for slightly greater
        # efficiency). Instead, we just compute it.
        # P <- calcR(beta)
        if (betaIsNamed) {
          P <- calcR(beta)
        } else {
          P <- calcR(x)
        }
        # DSB note: We do not do any error checking here on calcJ values.
        # Problems typically occur with calcR values, which are handled
        # prior to arriving here. We should consider whether to add
        # something here.
        # dr_vec <- calcJ(beta)
        if (betaIsNamed) {
          dr_vec <- calcJ(beta)
        } else {
          dr_vec <- calcJ(x)
        }
        dr_vec <- as.numeric(dr_vec)
        dim(dr_vec) <- c(n,p)
        # PP <- rep(r_vec[r_ind],p)
        P <- rep(P,p)
        dim(P) <- c(n,p)
        # print(PP)
        dr_vec <- t(dr_vec)*t(P)
      } else {
        # % 50 in dglf
        # %  ***  COMPUTE FINITE-DIFFERENCE APPROXIMATION TO
        # %       DR = GRAD. OF R  ***
        #   %     *** INITIALIZE D IF NECESSARY ***
        #   % display('Computing FD gradient...')
        # if(iv[mode] < 0 && v[dinit] == zero) {
        if ( (iv[mode_iv] < 0) && (v[dinit_v] == 0.) ) {
          # DSB FIX ME?  The following could be a loose end that
          # needs to addressed.  It is scaling related.
          # % [p, v(d1), one]=dv7scp(p, v(d1), one);
          # v[d_index] = rep(1,p)
          d_vec = rep(1,p)
        }
        ng <- iv[ngcall_iv] - 1
        if (iv[1] == -1) {
          # if iv(1) == -1
          # DSB Note:  This is related to covariance calculations
          # performed after convergence has been achieved.
          # iv(ngcov) = fix(iv(ngcov) - 1);
          # Need to double check this.
          iv[ngcov_iv] = iv[ngcov_iv] - 1
        }

        # Corresponds to code at line 70 of dglf
        dk <- d1_v # This is currently unused, because we don't address scaling
        j1k0 <- dr1_v
        stopLoop <- 0
        nf <- iv[nfcall_iv]
        if (nf != iv[nfgcal_iv]) {
          # If these do not match, it means we need to generate a current
          # copy of r and store it in rs1.
          # print(paste(" Updating prob at iteration ",iv[31]))
          ng <- ng + 1
          # P <- calcR(beta)
          if (betaIsNamed) {
            P <- calcR(beta)
          } else {
            P <- calcR(x)
          }
          r_all_vec[rs1_index] <- as.numeric(P)
          v[rs1_index_v] <- r_all_vec[rs1_index]
          # Remark: nf is set to zero if probabilities cannot be computed.
          if (any(is.na(P))) { nf <- 0}
          if (nf == 0) {
            iv[toobig_iv] <- 1
            iv[ngcall_iv] <- ng
            stopLoop <- 1
            computeDeriv <- 0
          }
        }
        for (k in 1:ps) {
          if (stopLoop == 1) break
          # DSB FIX ME: The following includes references to v[dk].
          # Investigate if this is a loose end that needs to be addressed.
          # print(paste(" FD calculation for component k: ", k))
          xk <- x[k]
          # h <- v[dltfdj]*max(abs(xk),1/v[dk])
          hk <- v[dltfdj_v]*max(abs(xk),1)
          # print(v[dk])
          h0 <- hk
          # go ahead and increment dk to be used in the next iteration
          # dk <- dk + 1
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
            if (betaIsNamed) {
              beta = copy_XToBeta(x,beta,p)
              P <- calcR(beta)
            } else {
              P <- calcR(x)
            }
            if (any(is.na(P))) { nf <- 0}
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
              if ( abs(hk/h0) < hlim) {
                iv[toobig_iv] <- 1
                compute_r21_loop <- 0
                computeDeriv <- 0
                break
              }
            } # end nf check
          } # end compute r21 loop
          # Finish up calculations for k loop, if things are still good
          x[k] <- xk
          iv[nfgcal_iv] <- ng
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

  g_vec <- v[g_index_v]
  if (!silent)
    iv <- bgw_itsum(d_vec, g_vec, iv, v, x, p, betaIsNamed, betaNames, i_itsum)

  model <- list()
  model$bgw_settings    <- bgw_settings
  model$hasAnalyticGrad <- hasAnalyticGrad
  model$betaStart       <- betaStart
  model$numParams       <- p
  model$numResids       <- n

  if(iv[1] >= 3 & iv[1] <= 15) {
    model$code <- iv[1]
    model$message <- trimws(bgw_message_table[[iv[1]]])
  } else {
    model$code <- 0
    model$message <-"Problematic BGW return code. Seek expert (maintainer) help."
  }

  # DSB note: We need to consider what to do with some of these
  # values if there is a catastrophic failure, e.g., if beta
  # has garbage.
  model$betaStop        <- beta
  f_iv                  <- 10
  model$finalLL         <- -v[f_iv]
  # model$LLOut           <- -v[f_iv]
  niter_iv              <- 31
  model$iterations      <- iv[niter_iv]
  nfcall_iv             <- 6
  model$functionEvals   <- iv[nfcall_iv]
  model$gradient        <- v[g_index_v]
  model$scaleVec        <- v[d_index_v]

  if ( (iv[1] >= 3) && (iv[1] <= 6) ) {
    model$estimate <- beta
    model$maximum  <- model$finalLL
  } else {
    model$estimate <- beta
    model$maximum  <- model$finalLL
  }

  covmat_iv <- 26
  rcond_iv <- 53
  vcHessianConditionNumber <- v[rcond_iv]^2
  if (iv[covmat_iv] > 0) {
    hessianMethodAttempted <- bgw_settings$vcHessianMethod
    hessianMethodUsed      <- bgw_settings$vcHessianMethod
    varcov <- matrix(NA, nrow=length(beta), ncol=length(beta),
                     dimnames=list(names(beta), names(beta)))
    se           <- rep(0,p)
    names(se)    <- names(beta)
    tstat        <- rep(0,p)
    names(tstat) <- names(beta)
    vcVec <- v[iv[26]:(iv[26] + p*(p+1)/2 - 1)]
    iii <- 0
    for (ii in 1:p) {
      for (jj in 1:ii) {
        iii <- iii + 1
        varcov[ii,jj] = vcVec[iii];
        varcov[jj,ii] = vcVec[iii];
      }
      se[ii] <- sqrt(varcov[ii,ii])
      tstat[ii] <- beta[ii]/se[ii]
    }
  } else {
    rdreq <- 57
    if (iv[rdreq] == 0) {
      # cat('       No variance-covariance matrix computation was requested. \n')
      hessianMethodUsed <- "None"
      vcHessianConditionNumber <- -1
      cat("\n")
    } else {
    hessianMethodAttempted <- bgw_settings$vcHessianMethod
    hessianMethodUsed      <- "N/A - BHHH Hessian is singular"
    }
    vcVec  <- NULL
    varcov <- NULL
    se     <- NULL
    tstat  <- NULL
  }

  model$vcHessianConditionNumber <- vcHessianConditionNumber
  # model$hessianMethodUsedBGW     <- hessianMethodUsed
  model$hessianMethodUsed        <- hessianMethodUsed
  model$vcVec                    <- vcVec
  model$varcovBGW                <- varcov
  model$seBGW                    <- se
  model$tstatBGW                 <- tstat

  # model$bgw_iv <- iv
  # model$bgw_v <- v

  # model$LLout           <- -v[10]
  # model$nIter <- iv[31]
  # model$iv <- iv
  # model$v  <- v

  class(model) = "bgw_mle"
  # c("\n")
  # print(attributes(model))
  # c("\n")
  return(model)
}
