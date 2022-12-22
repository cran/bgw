test_that("Simple mnl example 1", {
  #' Simple MNL log-likelihood function assuming linear utility
  #' @param b Numeric K-long vector of parameters to be estimated.
  #' @param Y Numeric NxJ matrix. \code{Y[n,j]=1} if alt j selected in obs n.
  #' @param X List of J NxK matrices with explanatory variables.
  #' @return an N-vector of choice probabilities.

  mnl <- function(b, Y, X){
    ### Check input
    test <- is.vector(b) && is.matrix(Y) && is.list(X) && all(sapply(X, is.matrix))
    if(!test) stop("Arguments 'b', 'Y', 'X' must be a vector, matrix, ",
                   "and list of matrices, respectively.")
    N <- nrow(Y)
    K <- length(b)
    J <- length(X)
    test <- all(dim(Y)==c(N,J)) && all(sapply(X, nrow)==N) && all(sapply(X, ncol)==K)
    if(!test) stop("Dimensions of arguments 'b', 'Y', 'X' do not match.")

    ### Calculate and return MNL loglikelihood
    eV <- sapply(X, function(x) exp(x%*%b))
    p  <- rowSums(Y*eV)/rowSums(eV)
    return( p )
  }

  ### Generate synthetic data
  N <- 1000
  K <- 3
  J <- 4
  set.seed(27)
  X <- list()
  for(j in 1:J) X[[j]] <- matrix(runif(N*K), nrow=N, ncol=K)
  b<- round(runif(K, min=-1, max=1), 2)
  U <- sapply(X, function(x) x%*%b + -log(-log(runif(N))))
  Y <- U==apply(U, MARGIN=1, FUN=max)

  ### Create starting values for estimation
  # b0<- setNames(rep(0, K), paste0("b", 1:K))
  b0 <- rep(0, K)

  ### Estimate using bgw
  mnl2 <- function(b) mnl(b, Y, X) # necessary so Y and X do not need to be given

  # First test example.  CRAN reviewers prefer non-verbose output, so
  # we impose this on the initial example.
  # We do not prefer this as the default...
  test_settings = list()
  test_settings[["modelName"]] <- "Simple mnl example 1 - non-verbose"
  test_settings[["silent"]] <- TRUE
  model <- bgw_mle(calcR=mnl2, betaStart=b0, bgw_settings=test_settings)
  #
  # Print information indicating what happened.
  cat(model$bgw_settings[["modelName"]])
  cat("\n")
  cat("\nOutcome of estimation search:",model$message,"\n")
  cat("Number of iterations: ",model$iterations,"\n")
  cat("\nMLE parameter estimate:\n")
  cat(model$estimate)
  cat("\n")
  cat("\nEstimated t-ratios:\n")
  cat(model$tstatBGW)
  cat("\n")

  expect_equal(model$maximum, -1361.97, tolerance = 0.01)
})

test_that("Simple mnl model 1 - non-verbose", {
  #' Simple MNL loglikelihood function assuming linear utility
  #' @param b Numeric K-long vector of parameters to be estimated.
  #' @param Y Numeric NxJ matrix. \code{Y[n,j]=1} if alt j selected in obs n.
  #' @param X List of J NxK matrices with explanatory variables.
  #' @return an N-vector of choice probabilities.

  mnl <- function(b, Y, X){
    ### Check input
    test <- is.vector(b) && is.matrix(Y) && is.list(X) && all(sapply(X, is.matrix))
    if(!test) stop("Arguments 'b', 'Y', 'X' must be a vector, matrix, ",
                   "and list of matrices, respectively.")
    N <- nrow(Y)
    K <- length(b)
    J <- length(X)
    test <- all(dim(Y)==c(N,J)) && all(sapply(X, nrow)==N) && all(sapply(X, ncol)==K)
    if(!test) stop("Dimensions of arguments 'b', 'Y', 'X' do not match.")

    ### Calculate and return MNL loglikelihood
    eV <- sapply(X, function(x) exp(x%*%b))
    p  <- rowSums(Y*eV)/rowSums(eV)
    return( p )
  }

  ### Generate synthetic data
  N <- 1000
  K <- 3
  J <- 4
  set.seed(27)
  X <- list()
  for(j in 1:J) X[[j]] <- matrix(runif(N*K), nrow=N, ncol=K)
  b<- round(runif(K, min=-1, max=1), 2)
  U <- sapply(X, function(x) x%*%b + -log(-log(runif(N))))
  Y <- U==apply(U, MARGIN=1, FUN=max)

  ### Create starting values for estimation
  b0<- setNames(rep(0, K), paste0("b", 1:K))

  ### Estimate using bgw
  mnl2 <- function(b) mnl(b, Y, X) # necessary so Y and X do not need to be given

  test_settings = list()
  test_settings[["printLevel"]] <- 2L
  model <- bgw_mle(calcR=mnl2, betaStart=b0, bgw_settings = test_settings)

  expect_equal(model$maximum, -1361.97, tolerance = 0.01)
})

test_that("Using settings to turn off VC calculation for mnl", {
  #' Simple MNL loglikelihood function assuming linear utility
  #' @param b Numeric K-long vector of parameters to be estimated.
  #' @param Y Numeric NxJ matrix. \code{Y[n,j]=1} if alt j selected in obs n.
  #' @param X List of J NxK matrices with explanatory variables.

  mnl <- function(b, Y, X){
    ### Check input
    test <- is.vector(b) && is.matrix(Y) && is.list(X) && all(sapply(X, is.matrix))
    if(!test) stop("Arguments 'b', 'Y', 'X' must be a vector, matrix, ",
                   "and list of matrices, respectively.")
    N <- nrow(Y)
    K <- length(b)
    J <- length(X)
    test <- all(dim(Y)==c(N,J)) && all(sapply(X, nrow)==N) && all(sapply(X, ncol)==K)
    if(!test) stop("Dimensions of arguments 'b', 'Y', 'X' do not match.")

    ### Calculate and return MNL loglikelihood
    eV <- sapply(X, function(x) exp(x%*%b))
    p  <- rowSums(Y*eV)/rowSums(eV)
    return( p )
  }

  ### Generate synthetic data
  N <- 1000
  K <- 3
  J <- 4
  set.seed(27)
  X <- list()
  for(j in 1:J) X[[j]] <- matrix(runif(N*K), nrow=N, ncol=K)
  b<- round(runif(K, min=-1, max=1), 2)
  U <- sapply(X, function(x) x%*%b + -log(-log(runif(N))))
  Y <- U==apply(U, MARGIN=1, FUN=max)

  ### Create starting values for estimation
  b0<- setNames(rep(0, K), paste0("b", 1:K))

  ### Estimate using bgw
  mnl2 <- function(b) mnl(b, Y, X) # necessary so Y and X do not need to be given

  test_settings = list()
  test_settings[["printLevel"]] <- 2L
  test_settings[["vcHessianMethod"]] <- "none"
  model <- bgw_mle(calcR=mnl2, betaStart=b0, bgw_settings = test_settings)

  expect_equal(model$maximum, -1361.97, tolerance = 0.01)
})
