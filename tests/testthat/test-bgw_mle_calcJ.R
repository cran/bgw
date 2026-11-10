# Tests for analytic derivatives (calcJ, calcJ_type),
# individual-level gradients, and related fixes in bgw 0.1.5.
# All tests use a small simulated multinomial logit (MNL) model, for which the
# exact derivatives are known. Added 7 November 2026.

# ---- Helpers ---------------------------------------------------------------

make_mnl_data <- function(N, K = 3, J = 4, seed = 27) {
  set.seed(seed)
  X <- lapply(1:J, function(j) matrix(runif(N*K), nrow = N, ncol = K))
  b <- round(runif(K, min = -1, max = 1), 2)
  U <- sapply(X, function(x) x %*% b - log(-log(runif(N))))
  Y <- U == apply(U, MARGIN = 1, FUN = max)
  list(X = X, Y = Y)
}

# n-vector of choice probabilities R_i
mnl_prob <- function(b, D) {
  eV <- sapply(D$X, function(x) exp(x %*% b))
  as.numeric(rowSums(D$Y*eV)/rowSums(eV))
}

# n x p matrix of d log(R_i) / d b  (the MNL scores)
mnl_score <- function(b, D) {
  eV <- sapply(D$X, function(x) exp(x %*% b))
  pr <- eV/rowSums(eV)
  xChosen <- Reduce(`+`, Map(function(x, j) x * D$Y[, j], D$X, seq_along(D$X)))
  xMean   <- Reduce(`+`, Map(function(x, j) x * pr[, j],  D$X, seq_along(D$X)))
  xChosen - xMean
}

D  <- make_mnl_data(1000)
b0 <- c(b1 = 0, b2 = 0, b3 = 0)
fR <- function(b) mnl_prob(b, D)
fS <- function(b) mnl_score(b, D)            # n x p, derivatives of log(R)
fJ <- function(b) mnl_score(b, D) * fR(b)    # n x p, derivatives of R

quiet <- function(...) c(list(silent = TRUE), list(...))

fit <- function(calcJ = NULL, calcJ_type = NULL, ..., calcR = fR, start = b0) {
  bgw_mle(calcR = calcR, betaStart = start, calcJ = calcJ,
          bgw_settings = quiet(...), calcJ_type = calcJ_type)
}

m_fd <- fit()   # finite-difference reference solution

# ---- All forms of calcJ give the same solution ----------------------------

test_that("analytic derivatives in every form agree with finite differences", {
  expect_true(m_fd$code %in% 3:6)
  forms <- list(
    list(J = fS,                    type = "log_probability"),
    list(J = function(b) t(fS(b)),  type = "log_probability"),
    list(J = fJ,                    type = "probability"),
    list(J = function(b) t(fJ(b)),  type = "probability")
  )
  for (f in forms) {
    m <- fit(f$J, f$type)
    expect_true(m$code %in% 3:6)
    expect_equal(m$maximum, m_fd$maximum, tolerance = 1e-8)
    expect_equal(m$estimate, m_fd$estimate, tolerance = 1e-5)
  }
})

# ---- Backward compatibility: calcJ without calcJ_type ---------------------

test_that("calcJ without calcJ_type assumes log_probability, without a warning", {
  expect_warning(m <- fit(fS), NA)
  expect_equal(m$maximum, m_fd$maximum, tolerance = 1e-8)
})

test_that("no warning when calcJ_type is given, or when calcJ is NULL", {
  expect_warning(fit(fS, "log_probability"), NA)
  expect_warning(fit(NULL, "probability"), NA)   # calcJ_type ignored without calcJ
})

# ---- Errors ----------------------------------------------------------------

test_that("invalid calcJ_type is an error", {
  expect_error(fit(fS, "probabilities"), "calcJ_type")
  expect_error(fit(fS, c("probability", "log_probability")), "calcJ_type")
})

test_that("calcJ returning a vector or wrong dimensions is an error", {
  expect_error(fit(function(b) as.numeric(fS(b)), "log_probability"), "not a vector")
  expect_error(fit(function(b) fS(b)[, -1], "log_probability"), "expected")
})

test_that("n = p: a square matrix is read as n x p", {
  D3  <- make_mnl_data(3, seed = 5)
  fR3 <- function(b) mnl_prob(b, D3)
  fS3 <- function(b) mnl_score(b, D3)
  m_nxp <- fit(fS3, "log_probability", calcR = fR3, start = b0, maxIterations = 3)
  m_def <- fit(fS3, calcR = fR3, start = b0, maxIterations = 3)
  expect_identical(m_nxp$estimate, m_def$estimate)
})

# ---- Individual-level gradients ---------------------------------------------

test_that("indLevelGradients are n x p derivatives of log(R) at the estimate", {
  m_a <- fit(fS, "log_probability", userRequestIndLevelGradients = TRUE)
  G   <- m_a$indLevelGradients
  expect_equal(dim(G), c(1000, 3))
  expect_equal(colnames(G), names(b0))
  expect_equal(unname(G), unname(fS(m_a$estimate)), tolerance = 1e-12)

  m_p <- fit(fJ, "probability", userRequestIndLevelGradients = TRUE)
  expect_equal(m_p$indLevelGradients, G, tolerance = 1e-10)

  m_f <- fit(userRequestIndLevelGradients = TRUE)
  expect_equal(unname(m_f$indLevelGradients), unname(fS(m_f$estimate)), tolerance = 1e-5)
})

test_that("indLevelGradients are not returned unless requested", {
  expect_null(m_fd$indLevelGradients)
})

# ---- Fixes -----------------------------------------------------------------

test_that("requesting indLevelGradients does not change the estimate (named beta, FD)", {
  m <- fit(userRequestIndLevelGradients = TRUE)
  expect_identical(m$estimate, m_fd$estimate)
  expect_identical(m$varcovBGW, m_fd$varcovBGW)
})

test_that("only invalid setting names no longer causes an internal error", {
  expect_error(capture.output(m <- bgw_mle(fR, b0, bgw_settings = list(notASetting = 1))), NA)
  expect_equal(m$maximum, m_fd$maximum, tolerance = 1e-8)
})

test_that("hessianMethodAttemptedBGW is returned", {
  expect_equal(m_fd$hessianMethodAttemptedBGW, "bhhh")
  m <- fit(vcHessianMethod = "none")
  expect_equal(m$hessianMethodAttemptedBGW, "none")
})
