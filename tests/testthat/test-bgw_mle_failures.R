# Tests of failed evaluations of calcR and calcJ.

set.seed(4)
fl_n  <- 400
fl_x1 <- rnorm(fl_n)
fl_x2 <- rnorm(fl_n)
fl_x3 <- rnorm(fl_n)
fl_y  <- rbinom(fl_n, 1, plogis(0.8*fl_x1 - 0.5*fl_x2 + 0.2*fl_x3))
fl_P  <- function(b) {
  b <- unname(b)
  plogis(b[1]*fl_x1 + b[2]*fl_x2 + b[3]*fl_x3)
}
fl_R  <- function(b) ifelse(fl_y == 1, fl_P(b), 1 - fl_P(b))
fl_J  <- function(b) cbind(fl_x1, fl_x2, fl_x3)*(fl_y - fl_P(b))
fl_LL <- function(b) sum(log(fl_R(b)))
fl_b0 <- c(a = 0.1, b = 0.2, c = -0.1)
fl_Q  <- list(silent = TRUE)

# calcJ that fails (one element NaN, NA or Inf) from its k-th call onwards
fl_badJ <- function(k, bad = NaN, type = "log_probability") {
  calls <- 0
  function(b) {
    calls <<- calls + 1
    G <- fl_J(b)
    if (type == "probability") G <- G*fl_R(b)
    if (calls >= k) G[1, 2] <- bad
    G
  }
}

test_that("calcJ failing at betaStart stops at once with code 65", {
  for (type in c("log_probability", "probability")) {
    for (bad in list(NaN, NA, Inf)) {
      nR <- 0
      R <- function(b) { nR <<- nR + 1; fl_R(b) }
      m <- bgw_mle(R, fl_b0, calcJ = fl_badJ(1, bad, type), bgw_settings = fl_Q, calcJ_type = type)
      expect_equal(m$code, 65)
      expect_match(m$message, "Gradient could not be computed", fixed = TRUE)
      expect_equal(m$iterations, 0)
      expect_identical(m$estimate, fl_b0)
      expect_equal(m$maximum, fl_LL(fl_b0))
      expect_lt(nR, 10)
    }
  }
})

test_that("calcJ failing later returns the last iterate at which the gradient was computed", {
  for (k in 2:4) {
    m <- bgw_mle(fl_R, fl_b0, calcJ = fl_badJ(k), bgw_settings = fl_Q,
                 calcJ_type = "log_probability")
    expect_equal(m$code, 65)
    expect_equal(m$iterations, k - 2)
    expect_equal(m$maximum, fl_LL(m$estimate), tolerance = 1e-10)
    expect_equal(m$finalLL, m$maximum)
    expect_equal(m$gradient, unname(-colSums(fl_J(m$estimate))), tolerance = 1e-8)
    expect_identical(m$betaStop, m$estimate)
  }
})

test_that("the iteration files have no row for the iterate where the gradient failed", {
  d <- tempfile("bgwfl")
  dir.create(d)
  s <- c(fl_Q, list(writeIter = TRUE, writeIterGradient = TRUE, writeIterScaleVec = TRUE,
                    modelName = "m", outputDirectory = d))
  m <- bgw_mle(fl_R, fl_b0, calcJ = fl_badJ(4), bgw_settings = s, calcJ_type = "log_probability")
  for (suffix in c("iterations", "gradients", "scaleVecs")) {
    f <- utils::read.csv(file.path(d, paste0("m_", suffix, ".csv")), check.names = FALSE)
    expect_equal(max(f$Iteration), m$iterations)
  }
  it <- utils::read.csv(file.path(d, "m_iterations.csv"), check.names = FALSE)
  gr <- utils::read.csv(file.path(d, "m_gradients.csv"), check.names = FALSE)
  last <- nrow(it)
  expect_equal(as.numeric(it[last, c("a", "b", "c")]), unname(m$estimate), tolerance = 1e-10)
  expect_equal(as.numeric(gr[last, c("a", "b", "c")]), -m$gradient, tolerance = 1e-8)
})

test_that("failed finite-difference steps later in the search are handled the same way", {
  for (K in c(8, 14)) {
    seen <- list()
    nR <- 0
    R <- function(b) {
      nR <<- nR + 1
      b <- unname(b)
      fd <- any(vapply(seen, function(s) {
        d <- abs(b - s)
        sum(d > 0) == 1 && max(d) < 1e-5
      }, TRUE))
      seen[[length(seen) + 1]] <<- b
      if (nR > K && fd) rep(NA_real_, fl_n) else fl_R(b)
    }
    m <- bgw_mle(R, fl_b0, bgw_settings = fl_Q)
    expect_equal(m$code, 65)
    expect_equal(m$maximum, fl_LL(m$estimate), tolerance = 1e-10)
    expect_equal(m$gradient, unname(-colSums(fl_J(m$estimate))), tolerance = 1e-4)
  }
})

test_that("calcR returning Inf or a value <= 0 at betaStart is rejected", {
  for (bad in list(0, -0.1, Inf, NaN, NA)) {
    R <- function(b) {
      p <- fl_R(b)
      p[1] <- bad
      p
    }
    expect_error(bgw_mle(R, fl_b0, bgw_settings = fl_Q), "Initial call to calcR failed")
  }
})

test_that("calcR returning Inf during the search no longer produces NaN parameters", {
  R <- function(b) {
    p <- fl_R(b)
    if (unname(b)[1] > 0.5) p[1] <- Inf
    p
  }
  for (analytic in c(TRUE, FALSE)) {
    m <- if (analytic) {
      bgw_mle(R, fl_b0, calcJ = fl_J, bgw_settings = fl_Q, calcJ_type = "log_probability")
    } else {
      bgw_mle(R, fl_b0, bgw_settings = fl_Q)
    }
    expect_true(all(is.finite(m$estimate)))
    expect_lte(unname(m$estimate)[1], 0.5)
    expect_equal(m$maximum, fl_LL(m$estimate), tolerance = 1e-8)
  }
})
