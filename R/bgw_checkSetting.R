#' @title bgw_checkSetting
#'
#' @description Checks to see if user-provided value for a bgw_setting is valid.
#'
#' @param settingValue Setting value (submitted by user).
#' @param bgw_setting_type Type of setting being checked. Possible values are "discrete" and "continuous".
#' @param bgw_validDiscrete List. Contains valid values for a discrete setting.
#' @param bgw_contLB Numerical value. The lower bound for a valid continuous setting.
#' @param bgw_contUB Numerical value. The upper bound for a valid continuous setting.
#'
#' @return Logical. Indicates if the setting is okay or not.
#' @export
# NOTE:  I need to check and see if @export is only included for user-accessible functions.
##------------------------------------------------------------------------
# October 14, 2022
bgw_checkSetting <- function(settingValue, bgw_setting_type, bgw_validDiscrete,
                             bgw_contLB, bgw_contUB) {
  settingOk <- TRUE
  if (bgw_setting_type == 'discrete') {
    settingOk <- FALSE
    for (i in bgw_validDiscrete) {
      if (settingValue == i) {
        settingOk <- TRUE
      }
    }
  }

  if (bgw_setting_type == 'continuous') {
    settingOk <- FALSE
    if ( (settingValue >= bgw_contLB)&&(settingValue <= bgw_contUB) ) {
      settingOk <- TRUE
    }
  }
  return(settingOk)
}
