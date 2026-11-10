# Internal functions that handle bgw_settings, in the order bgw_mle calls them:
#   bgw_defineSettings   what can be set: defaults, types, valid values and bounds
#   bgw_resolveSettings  combines those with what the caller supplied, and checks the result
#   bgw_applySettings    writes the final settings into the BGW work vectors iv and v
# None of these is exported, and no help page is generated.

#' @title bgw_defineSettings
#'
#' @description Defines the settings a user can change.
#'
#' @return List with elements default, type, validDiscrete, continuousLB and continuousUB, each a
#'   list indexed by setting name.
#' @noRd
bgw_defineSettings <- function() {
  # First, default bgw_settings (for those values that can be altered by the user).
  bgw_settings_default        <- list()
  bgw_settings_type           <- list()
  bgw_settings_validDiscrete  <- list()
  bgw_settings_continuousLB   <- list()
  bgw_settings_continuousUB   <- list()
  #
  # Printing/writing options
  # Fortran BGW has a relatively complex subroutine DITSUM that produces output.
  # A variety of control parameters in iv are available to select options.
  # Depending on the options, DITSUM writes out the initial beta and scaling vector,
  # an iteration summary line for each iteration, final status, and final beta
  # vector and gradient (which may or may not correspond to a converged solution).
  #
  # Because R does not allow any Fortran I/O statements, this routine had to be
  # replaced with an R function (bgw_itsum.R), which writes all DITSUM output
  # to the R console. One possible future option would be to also save the
  # same output to a file.
  #
  # DITSUM standard output contains two iteration summary formats: a long line
  # (the default) and a short line. Patterned after Apollo, the current
  # standard output (long line, the most verbose) is assigned printLevel = 3L.
  # We have also implemented the short line (printLevel = 2L.) The option
  # printLevel = 1L is also assigned to the short line, but in the future could
  # be modified to provide an even less verbose option. The option printLevel = 0L
  # suppresses all DITSUM-related output.
  # We have some redundancy. "silent" is the same as "printLevel = 0L".
  bgw_settings_default[["silent"]]             <- FALSE
  bgw_settings_type[["silent"]]                <- "discrete"
  bgw_settings_validDiscrete[["silent"]]       <- c(FALSE, TRUE)
  bgw_settings_continuousLB[["silent"]]        <- 0
  bgw_settings_continuousUB[["silent"]]        <- 0
  #
  bgw_settings_default[["printLevel"]]         <- 3L
  bgw_settings_type[["printLevel"]]            <- "discrete"
  bgw_settings_validDiscrete[["printLevel"]]   <- c(0L, 1L, 2L, 3L)
  bgw_settings_continuousLB[["printLevel"]]    <- 0
  bgw_settings_continuousUB[["printLevel"]]    <- 0
  #
  bgw_settings_default[["printNonDefaultSettings"]]         <- TRUE
  bgw_settings_type[["printNonDefaultSettings"]]            <- "discrete"
  bgw_settings_validDiscrete[["printNonDefaultSettings"]]   <- c(FALSE, TRUE)
  bgw_settings_continuousLB[["printNonDefaultSettings"]]    <- 0
  bgw_settings_continuousUB[["printNonDefaultSettings"]]    <- 0
  #
  bgw_settings_default[["printStartingValues"]]         <- TRUE
  bgw_settings_type[["printStartingValues"]]            <- "discrete"
  bgw_settings_validDiscrete[["printStartingValues"]]   <- c(FALSE, TRUE)
  bgw_settings_continuousLB[["printStartingValues"]]    <- 0
  bgw_settings_continuousUB[["printStartingValues"]]    <- 0
  #
  bgw_settings_default[["printFinalResults"]]         <- TRUE
  bgw_settings_type[["printFinalResults"]]            <- "discrete"
  bgw_settings_validDiscrete[["printFinalResults"]]   <- c(FALSE, TRUE)
  bgw_settings_continuousLB[["printFinalResults"]]    <- 0
  bgw_settings_continuousUB[["printFinalResults"]]    <- 0

  # Note that BGW uses multiple flags to control different aspects of printing results.
  # Some of these could be added to bgw_itsum at a later time. Either way,
  # additional documentation to explain the iteration summary will be required.
  # (DITSUM output can also be written to a file. See below.)
  #
  # Stopping tests
  # maximum number of function evaluations. BGW (Fortran) default is 200.
  # Re-defining default for BGW R package to 300.
  bgw_settings_default[["maxFunctionEvals"]]        <- 300
  bgw_settings_type[["maxFunctionEvals"]]           <- "continuous"
  bgw_settings_validDiscrete[["maxFunctionEvals"]]  <- list()
  bgw_settings_continuousLB[["maxFunctionEvals"]]   <- 1
  # Capped at the largest integer, because iv holds this value as a Fortran integer.
  bgw_settings_continuousUB[["maxFunctionEvals"]]   <- .Machine$integer.max
  # maximum number of iterations. BGW (Fortran) default is 150.
  # Re-defining default for BGW R package to 200
  bgw_settings_default[["maxIterations"]]        <- 200
  bgw_settings_type[["maxIterations"]]           <- "continuous"
  bgw_settings_validDiscrete[["maxIterations"]]  <- list()
  # Zero requests a model object at betaStart (see userRequestModelObject).
  bgw_settings_continuousLB[["maxIterations"]]   <- 0
  # Capped at the largest integer, because iv holds this value as a Fortran integer.
  bgw_settings_continuousUB[["maxIterations"]]   <- .Machine$integer.max
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
  bgw_settings_default[["writeIter"]]        <- FALSE
  bgw_settings_type[["writeIter"]]           <- "discrete"
  bgw_settings_validDiscrete[["writeIter"]]  <- c(FALSE, TRUE)
  # If true, the file is "modelName_iterations.csv."
  bgw_settings_continuousLB[["writeIter"]]   <- 0
  bgw_settings_continuousUB[["writeIter"]]   <- 0

  # Additional option related to writeIter:  overWrite or append?
  bgw_settings_type[["writeIterMode"]]           <- "discrete"
  bgw_settings_validDiscrete[["writeIterMode"]]  <- c("overWrite", "append")
  bgw_settings_default[["writeIterMode"]]        <- "overWrite"
  bgw_settings_continuousLB[["writeIterMode"]]   <- 0
  bgw_settings_continuousUB[["writeIterMode"]]   <- 0
  # writeIterGradient and writeIterScaleVec write the gradient of the log-likelihood and the scale
  # vector at each iteration, in the same layout as writeIter and following writeIterMode.
  bgw_settings_default[["writeIterGradient"]]        <- FALSE
  bgw_settings_type[["writeIterGradient"]]           <- "discrete"
  bgw_settings_validDiscrete[["writeIterGradient"]]  <- c(FALSE, TRUE)
  bgw_settings_continuousLB[["writeIterGradient"]]   <- 0
  bgw_settings_continuousUB[["writeIterGradient"]]   <- 0
  bgw_settings_default[["writeIterScaleVec"]]        <- FALSE
  bgw_settings_type[["writeIterScaleVec"]]           <- "discrete"
  bgw_settings_validDiscrete[["writeIterScaleVec"]]  <- c(FALSE, TRUE)
  bgw_settings_continuousLB[["writeIterScaleVec"]]   <- 0
  bgw_settings_continuousUB[["writeIterScaleVec"]]   <- 0

  # Variance-covariance (vc) matrix calculation
  # How iv[rdreq] and iv[covreq] are set is documented in bgw_applySettings.
  bgw_settings_default[["vcHessianMethod"]]       <- "bhhh"
  bgw_settings_type[["vcHessianMethod"]]          <- "discrete"
  bgw_settings_validDiscrete[["vcHessianMethod"]] <- c("none", "bhhh",
                                                       "finiteDifferences", "fdFunction")
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
  # The type of scaling is controlled by iv[dtype] = iv[16].
  # The options are:
  # iv[dtype] = 0 => do not update d (see dinit discussion below).
  #           = 1 => update d on every iteration
  #           = 2 => update d on first iteration only.
  # BGW default is iv[dtype] = 1, which indicates adaptive scaling.
  #
  # Remark: Our experience is that many problems are reasonably well
  # scaled, and setting the scaling vector d = 1 often works well.
  # The vector d can be set to 1 by setting
  #     v[dinit] = v[38] = 1 and iv[dtype] = 0 (see above).
  # Default settings for adaptive scaling were already set by bgw_mle_setup.
  bgw_settings_default[["scalingMethod"]]        <- "adaptive"
  bgw_settings_type[["scalingMethod"]]           <- "discrete"
  bgw_settings_validDiscrete[["scalingMethod"]]  <- c("adaptive", "none", "userScaling")
  bgw_settings_continuousLB[["scalingMethod"]]   <- 0
  bgw_settings_continuousUB[["scalingMethod"]]   <- 0
  #
  bgw_settings_default[["userScaleVector"]]        <- "NA"
  bgw_settings_type[["userScaleVector"]]           <- "NA"
  bgw_settings_validDiscrete[["userScaleVector"]]  <- list()
  bgw_settings_continuousLB[["userScaleVector"]]   <- 0
  bgw_settings_continuousUB[["userScaleVector"]]   <- 0
  #
  # For bgw 0.1.4, Apollo developers requested a capability for bgw to
  # write out "individual level gradients" at the request of the user.
  # Individual level gradients are terms of the form d[log(R_i)]/dbeta, that
  # when summed yields the gradient of LL.
  bgw_settings_default[["userRequestIndLevelGradients"]]        <- FALSE
  bgw_settings_type[["userRequestIndLevelGradients"]]           <- "discrete"
  bgw_settings_validDiscrete[["userRequestIndLevelGradients"]]  <- c(FALSE, TRUE)
  bgw_settings_continuousLB[["userRequestIndLevelGradients"]]   <- 0
  bgw_settings_continuousUB[["userRequestIndLevelGradients"]]   <- 0
  # Apollo developers requested a way to obtain a model object at betaStart without estimation.
  # bgw_resolveSettings also turns this on when maxIterations is 0.
  bgw_settings_default[["userRequestModelObject"]]        <- FALSE
  bgw_settings_type[["userRequestModelObject"]]           <- "discrete"
  bgw_settings_validDiscrete[["userRequestModelObject"]]  <- c(FALSE, TRUE)
  bgw_settings_continuousLB[["userRequestModelObject"]]   <- 0
  bgw_settings_continuousUB[["userRequestModelObject"]]   <- 0
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
  # bgw_settings_default["absFunctionConvTol"] <- v[bgw_v_index$afctol]
  # [Add more here...]
  #
  #
  # Other options could be added, but it is possible to write other functions
  # for displaying/summarizing results based on the contents of the model object
  # AFTER the estimation has been completed.
  # -------------------------------------------------------------------------------------

  return(list(default       = bgw_settings_default,
              type          = bgw_settings_type,
              validDiscrete = bgw_settings_validDiscrete,
              continuousLB  = bgw_settings_continuousLB,
              continuousUB  = bgw_settings_continuousUB))
}

#' @title bgw_resolveSettings
#'
#' @description Combines the default settings with those supplied by the caller, reporting or
#'   rejecting anything that is not valid.
#'
#' @param bgw_settings List of settings passed to bgw_mle (may be NULL).
#' @param settingDefs List returned by bgw_defineSettings.
#' @param hasAnalyticGrad Logical. Whether the caller supplied calcJ (only used to report which
#'   derivatives are in use).
#' @param p Number of parameters.
#'
#' @return List with elements bgw_settings (final settings), modifiedSettings (names of settings
#'   that differ from the defaults) and silent.
#' @noRd
bgw_resolveSettings <- function(bgw_settings, settingDefs, hasAnalyticGrad, p) {
  bgw_settings_default       <- settingDefs$default
  bgw_settings_type          <- settingDefs$type
  bgw_settings_validDiscrete <- settingDefs$validDiscrete
  bgw_settings_continuousLB  <- settingDefs$continuousLB
  bgw_settings_continuousUB  <- settingDefs$continuousUB
  silent           <- bgw_settings_default[["silent"]]
  printNDS         <- bgw_settings_default[["printNonDefaultSettings"]]
  modifiedSettings <- list()

  # Establishing final bgw_settings
  #
  # If no bgw_settings are in argument list (or if a passed bgw_settings is NULL)
  #    => use default values for bgw_settings
  # Otherwise, bgw_settings are user-provided and must be checked.
  check_bgw_settings <- !((length(bgw_settings) == 0) ||
                          (length(bgw_settings) == 1 && is.na(bgw_settings)) ||
                          (is.null(bgw_settings)))

  # silent must be known before anything is printed, so look for it first.
  if (check_bgw_settings) {
    # We need to check for "silent" first, before proceeding.
    if (!is.null(bgw_settings$printLevel)) {
      if (bgw_settings$printLevel == 0L) {
        bgw_settings$silent <- TRUE
        silent <- TRUE
      }
    }
    if (!is.null(bgw_settings$printNonDefaultSettings)) {
      if (bgw_settings$printNonDefaultSettings == FALSE) {
        printNDS <- FALSE
      }
    }
    if (!is.null(bgw_settings$silent)) {
      if (bgw_settings$silent) {
        bgw_settings$printLevel <- 0L
        silent <- TRUE
        bgw_settings$printNonDefaultSettings <- FALSE
        printNDS <- FALSE
        bgw_settings$printStartingValues <- FALSE
        bgw_settings$printFinalResults <- FALSE
      }
    }
  }

  if (!silent) {
    if (hasAnalyticGrad) {
      cat("\nBGW using analytic model derivatives supplied by caller...\n")
    } else {
      cat(paste0("\nBGW is using FD derivatives for model Jacobian. (Caller did not provide ",
                 "derivatives.)\n"))
    }
  }

  if (!check_bgw_settings) {
    bgw_settings <- bgw_settings_default
    cat("\nThere are no user-provided BGW settings. \n")
    cat("BGW settings have been set to default values...\n")
  } else {
    # User has provided bgw_settings that must be checked.
    userSettings          <- bgw_settings
    userSettingNames      <- names(bgw_settings)
    numUserSettingNames   <- length(userSettingNames)
    # A list of validUserSettingNames is created, in the same order as they would be in
    # bgw_settings_default
    validUserSettingNames <-
      names(bgw_settings_default)[names(bgw_settings_default) %in% userSettingNames]
    numValidUserSettingNames <- length(validUserSettingNames)
    # print(validUserSettingNames)
    # print(numValidUserSettingNames)

    # A completely new version of bgw_settings is created.
    # It is initialized with bgw_settings_default, and then systematically updated
    # with valid user-provided settings.
    # Any valid user-provided settings that are different from default settings are
    # stored in modifiedSettings.  modifiedSettings require changes to BGW default values.
    bgw_settings <- bgw_settings_default
    # 7 November 2026: bug fix. modifiedSettings is now initialized here rather
    # than inside the if-block below. Previously, if every user-provided setting
    # name was invalid (numValidUserSettingNames == 0), modifiedSettings was never
    # created and length(modifiedSettings) further down failed with an error.
    modifiedSettings <- list()
    if (numValidUserSettingNames > 0) {
      if (printNDS) cat("\nValid user-provided bgw_setting names(s) detected...\n")
      for (i in validUserSettingNames) {
        if (!setequal(bgw_settings[[i]], userSettings[[i]])) {
          if (bgw_settings_type[[i]] != "NA") {
            settingOk <- bgw_checkSetting(userSettings[[i]],
                                          bgw_settings_type[[i]], bgw_settings_validDiscrete[[i]],
                                          bgw_settings_continuousLB[[i]],
                                          bgw_settings_continuousUB[[i]])
          } else {
            settingOk <- TRUE
          }
          if (settingOk) {
            bgw_settings[[i]] <- userSettings[[i]]
            if (is.character(bgw_settings[[i]])) {
              if (printNDS) cat('\nbgw_settings[[\"', i, '\"]] is set to user-provided value (\"',
                                userSettings[[i]], '\").', sep = '')
            } else {
              if (printNDS) cat('\nbgw_settings[[\"', i, '\"]] is set to user-provided value (',
                                userSettings[[i]], ').', sep = '')
            }
            modifiedSettings <- c(modifiedSettings, i)
          } else {
            if (is.character(bgw_settings[[i]])) {
              if (printNDS) {
                cat('\nUser-provided value for bgw_settings[[\"', i, '\"]]  (\"', userSettings[[i]],
                    '\") is not valid.\n', sep = '')
                cat('    Default value (\"', bgw_settings[[i]], '\") will be used instead.\n',
                    sep = '')
              } else {
                stop('\nUser-provided value for bgw_settings[[\"', i, '\"]]  (\"',
                     userSettings[[i]], '\") is not valid.\n', sep = '')
              }
            } else {
              if (printNDS) {
                cat('\nUser-provided value for bgw_settings[[\"', i, '\"]] is not valid.\n')
                # 7 November 2026: bug fix. Previously printed userSettings[[i]] (the
                # rejected value) here instead of the default value.
                cat('    Default value (', bgw_settings[[i]], ') will be used instead.', sep = '')
              } else {
                stop('\nUser-provided value for bgw_settings[[\"', i, '\"]] is not valid.\n')
              }
            }
          }
        } else {
          if (is.character(bgw_settings[[i]])) {
            if (printNDS) cat('\nUser-provided value for bgw_settings[[\"', i,
                              ']]\" is the same as default (\"', userSettings[[i]], '\").',
                              sep = '')
          } else {
            if (printNDS) cat('\nUser-provided value for bgw_settings[[\"', i,
                              ']]\" is the same as default (', userSettings[[i]], ').', sep = '')
          }
        }
      }
      if (!silent || printNDS) cat("\n")
    }

    if ((numUserSettingNames > 0) && (numUserSettingNames != numValidUserSettingNames)) {
      badUserSettingNames <- userSettingNames[!(userSettingNames %in% validUserSettingNames)]
      for (i in badUserSettingNames) {
        if (printNDS) {
          cat('\nNon-valid user-provided bgw_setting setting name detected: \"', i, '\".\n',
              sep = "")
        } else {
          stop('\nNon-valid user-provided bgw_setting setting name detected: \"', i, '\".\n',
               sep = "")
        }
      }
      rm(badUserSettingNames)
      if (!silent || printNDS) cat("\n")
    }
    rm(userSettings, userSettingNames, numUserSettingNames, validUserSettingNames,
       numValidUserSettingNames)

    lms <- length(modifiedSettings)
    if (lms > 0) {
      # First, handle more complex situations.
      # Attempts to use a user-provided scale vector have potential issues.
      if (any(modifiedSettings == "scalingMethod")) {
        if ((bgw_settings[["scalingMethod"]] == "userScaling") &&
            (!any(modifiedSettings == "userScaleVector"))) {
              stop(paste0("Error. User has set bgw_settings to userScaling, but no ",
                          "userScaleVector has been provided."))
            }
      }
      if (any(modifiedSettings == "userScaleVector")) {
        userScaleVector <- bgw_settings[["userScaleVector"]]
        if (!is.numeric(userScaleVector))
          stop("Error in userScaleVector. Must be a numeric vector.")
        if (length(userScaleVector) != p) stop("Error in userScaleVector. It is the wrong length.")
        if (any(userScaleVector <= 0)) stop(paste0("Error in userScaleVector. One or more values ",
                                                   "is non-positive."))
        if (bgw_settings[["scalingMethod"]] != "userScaling") {
          stop(paste0("\nUser-provided scale vector detected, but user has not set scalingMethod ",
                      "to userScaling \n"))
          # bgw_settings[["scalingMethod"]] <- "userScaling"
          # cat("\nUser-provided scale vector detected. bgw_settings updated to reflect
          # scalingMethod = userScaling\n")
        }
      }
    }
  }

  # maxIterations = 0 is another way to request a model object at betaStart.
  if (bgw_settings[["maxIterations"]] < 1) bgw_settings[["userRequestModelObject"]] <- TRUE
  if (bgw_settings[["userRequestModelObject"]] &&
      !(bgw_settings[["vcHessianMethod"]] %in% c("bhhh", "none"))) {
    stop(paste0("A model object at betaStart (userRequestModelObject = TRUE or maxIterations = 0) ",
                'requires vcHessianMethod = "bhhh" or "none", but vcHessianMethod is "',
                bgw_settings[["vcHessianMethod"]], '".'))
  }

  return(list(bgw_settings = bgw_settings, modifiedSettings = modifiedSettings, silent = silent))
}

#' @title bgw_applySettings
#'
#' @description Writes the final settings into the BGW work vectors.
#'
#' @param iv BGW internal vector of integer values.
#' @param v BGW internal vector of numeric values.
#' @param bgw_settings_default List of default settings (from bgw_defineSettings).
#' @param bgw_settings List of final settings (from bgw_resolveSettings).
#' @param modifiedSettings Names of settings that differ from the defaults.
#' @param silent Logical (from bgw_resolveSettings).
#' @param hasAnalyticGrad Logical. Whether the caller supplied calcJ.
#'
#' @return List with elements iv, v and silent.
#' @noRd
bgw_applySettings <- function(iv, v, bgw_settings_default, bgw_settings, modifiedSettings, silent,
                              hasAnalyticGrad) {
  # Package defaults that replace the BGW Fortran defaults.
  iv[bgw_iv_index$outlev] <- 1
  iv[bgw_iv_index$mxfcal] <- bgw_settings_default[["maxFunctionEvals"]]
  iv[bgw_iv_index$mxiter] <- bgw_settings_default[["maxIterations"]]
  #
  # First, indicate whether or not to compute a vc matrix.
  # Note that BGW has capabilities for both vc matrices and more advanced
  # regression diagnostics. We only allow for vc matrices in this package.
  # This is controlled by iv[rdreq] = iv[57]
  # iv[rdreq] = 0 means do not compute anything.
  #           = 1 means compute VC matrix.
  # BGW default is 3 because  of the regression diagnostic
  # options.  We immediately override so that it is 1.
  iv[bgw_iv_index$rdreq] <- 1
  # This means that we will compute a VC matrix, using
  # the default method (BHHH) discussed next.
  # bgw_settings_default[["vcHessianMethod"]] <= "bhhh" (see below).
  #
  # Assuming a vc matrix is to be calculated, BGW needs to
  # know the type of Hessian approximation to use.
  # This is controlled by iv[covreq] = iv[15].
  # The options are:
  # iv[covreq] =  0, 1, 2 => a FD Hessian from gradient differences.
  #            = -1, -2   => a FD Hessian from function differences,
  #            = -3 or 3  => a BHHH Hessian (a.k.a. Gauss-Newton Hessian)
  #
  # Currently, if |iv[covreq]| = 1 or 2 then VC = H^(-1);
  #                            = 3      then VC = BHHH^(-1)
  # [An alternate version of the code computes robust VC's, but this is
  # not available in the 1993 BGW]
  # The BGW default value is 1, but if an analytic gradient is not available
  # it must be changed to -1. We are overriding this and changing it to 3.
  # However, note that the iv[rdreq] precludes calculating it.
  iv[bgw_iv_index$covreq] <- 3

  # Settings the caller changed.
  lms <- length(modifiedSettings)
  if (lms > 0) {
    # bgw_resolveSettings has confirmed that scalingMethod is "userScaling" in this case.
    if (any(modifiedSettings == "userScaleVector")) {
      iv[bgw_iv_index$dtype] <- 0
      v[bgw_v_index$dinit]   <- -1.e0
    }
    #
    for (i in 1:lms) {
      if (modifiedSettings[[i]] == "scalingMethod") {
        if (bgw_settings[["scalingMethod"]] == "none") {
          iv[bgw_iv_index$dtype] <- 0
          v[bgw_v_index$dinit]   <- 1.e0
        }
      }
      if (modifiedSettings[[i]] == "printStartingValues") {
        if (bgw_settings[["printStartingValues"]] == FALSE) {
          iv[bgw_iv_index$x0prt] <- 0
        }
      }
      if (modifiedSettings[[i]] == "printFinalResults") {
        if (bgw_settings$printFinalResults == FALSE) {
          iv[bgw_iv_index$solprt] <- 0
          iv[bgw_iv_index$statpr] <- -1
        }
      }
      if (modifiedSettings[[i]] == "vcHessianMethod") {
        if (bgw_settings[["vcHessianMethod"]] == "none") {
          iv[bgw_iv_index$rdreq] <- 0
        }
        if (bgw_settings[["vcHessianMethod"]] == "finiteDifferences") {
          if (hasAnalyticGrad) {
            iv[bgw_iv_index$covreq] <- 2
          } else {
            iv[bgw_iv_index$covreq] <- -2
          }
        }
        if (bgw_settings[["vcHessianMethod"]] == "fdFunction") {
          iv[bgw_iv_index$covreq] <- -2
        }
      }
      if (modifiedSettings[[i]] == "maxIterations") {
        iv[bgw_iv_index$mxiter] <- bgw_settings[["maxIterations"]]
      }
      if (modifiedSettings[[i]] == "maxFunctionEvals") {
        iv[bgw_iv_index$mxfcal] <- bgw_settings[["maxFunctionEvals"]]
      }
      # Some of the following is redundant
      if (modifiedSettings[[i]] == "printLevel") {
        if (bgw_settings[["printLevel"]] == 0L) {
          silent <- TRUE
          iv[bgw_iv_index$outlev] <- 0
        }
        if (bgw_settings[["printLevel"]] == 1L) {
          silent <- FALSE
          iv[bgw_iv_index$outlev] <- 0
        }
        if (bgw_settings[["printLevel"]] == 2L) {
          silent <- FALSE
          iv[bgw_iv_index$outlev] <- -1
        }
        if (bgw_settings[["printLevel"]] == 3L) {
          silent <- FALSE
          iv[bgw_iv_index$outlev] <- 1
        }
      }
    }
  }

  return(list(iv = iv, v = v, silent = silent))
}
