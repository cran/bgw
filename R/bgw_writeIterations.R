#' Writes one row, [values, ll, iteration, itType], to one of the iteration files. All of them
#' have the same layout, so that they can be stacked.
#' #
#' Was created using apollo_writeTheta as a starting point...
#' Because this is an internal function, the inputs will be assumed to be clean.
#' @param values vector to be written: the parameters, the gradient or the scale vector.
#' @param ll scalar representing the log-likelihood of the whole model.
#' @param iteration iteration number.
#' @param itType Character. What values contains: "estimate", "gradient" or "scaleVec".
#' @param valueNames Character vector. Column names for values: the parameter names, or beta1,
#'   beta2, ... if the parameters are not named.
#' @param outputFile Character. Name of the output file.
#' @return Nothing.
#' @noRd
# 7 November 2026: internal function (not exported; no help page is generated).
bgw_writeIterations <- function(values, ll, iteration, itType, valueNames, outputFile) {

  if (is.null(outputFile)) {
    stop("Attempt to call bgw_writeIterations with no file name.  Internal error.")
  }

  if (file.exists(outputFile)) {
    # If file already exists, append. A file from an earlier version of bgw, or written by Apollo,
    # has no ItType column; the rows must then match the layout it already has.
    header <- tryCatch(readLines(outputFile, n = 1), error = function(e) character(0))
    tmp <- matrix(c(values, ll), nrow = 1)
    if (any(grepl("ItType", header, fixed = TRUE))) {
      tmp <- as.data.frame(tmp)
      tmp$Iteration <- iteration
      tmp$ItType    <- itType
    }
    tryCatch(utils::write.table(tmp, file = outputFile, append = TRUE, sep = ',', col.names = FALSE,
                                row.names = FALSE),
             error = function(e) cat('Current iteration could not be written to ', outputFile,
                                     '.\n', sep = ''))
  } else {
    # If file does not exist, write afresh with headings
    tmp <- as.data.frame(matrix(c(values, ll), nrow = 1))
    tmp$Iteration <- iteration
    tmp$ItType    <- itType
    colnames(tmp) <- c(valueNames, 'logLike', 'Iteration', 'ItType')
    tryCatch(utils::write.table(tmp, file = outputFile, sep = ',', row.names = FALSE,
                                append = FALSE),
             error = function(e) cat('Initial iteration could not be written to ', outputFile,
                                     '.\n', sep = ''))
  }
}
