## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bgw)

## -----------------------------------------------------------------------------
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
  rm(U, j)

  ### Create starting values for estimation
  b0<- setNames(rep(0, K), paste0("b", 1:K))

  ### Estimate using bgw
  mnl2 <- function(b) mnl(b, Y, X) # necessary so Y and X do not need to be given
  model <- bgw_mle(calcR=mnl2, betaStart=b0)

