#' Writes the vector [beta,ll] to a file called \code{modelname_iterations.csv}
#' #
#' Was created using apollo_writeTheta as a starting point...
#' Because this is an internal function, the inputs will be assumed to be clean.
#' @param beta vector of parameters to be written (for now, no fixed betas).
#' @param ll scalar representing the log-likelihood of the whole model.
#' @param outputFile Character. Name of the output file.
#' @return Nothing.
#' @export
bgw_writeIterations <- function(beta, ll, outputFile){

  # One main issue (for now) is: Does beta have names, or is it only a numeric vector?
  # Answer: Beta will always have names.  If it starts with no names, we will assign names.
  if (is.null(outputFile)) {
    stop("Attempt to call bgw_writeIterations with no file name.  Internal error.")
  }

  ### Initialise
  tmp <- matrix(c(beta,ll ),nrow=1)

  if(file.exists(outputFile)){
    # If file already exists, append
    tryCatch( utils::write.table(tmp,file=outputFile, append=TRUE, sep=',', col.names=FALSE, row.names=FALSE),
              error=function(e) cat('Current iteration could not be written to ',outputFile,'.\n', sep='') )
  } else {
    # If file does not exist, write afresh with headings
    colnames(tmp) <- c(names(beta),'logLike')
    tryCatch( utils::write.table(tmp,file=outputFile, sep=',', row.names=FALSE, append=FALSE),
              error=function(e) cat('Initial iteration could not be written to ',outputFile,'.\n', sep='') )
  }
}
