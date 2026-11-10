# Tests of the model object returned at betaStart without estimation.

set.seed(2)
mo_n  <- 400
mo_x1 <- rnorm(mo_n)
mo_x2 <- rnorm(mo_n)
mo_x3 <- rnorm(mo_n)
mo_y  <- rbinom(mo_n, 1, plogis(0.8*mo_x1 - 0.5*mo_x2 + 0.2*mo_x3))
mo_P  <- function(b) {
  b <- unname(b)
  plogis(b[1]*mo_x1 + b[2]*mo_x2 + b[3]*mo_x3)
}
mo_R  <- function(b) ifelse(mo_y == 1, mo_P(b), 1 - mo_P(b))
mo_J  <- function(b) cbind(mo_x1, mo_x2, mo_x3)*(mo_y - mo_P(b))
mo_b0 <- c(a = 0.1, b = 0.2, c = -0.1)
mo_Q  <- list(silent = TRUE)

mo_fit <- function(settings = list(userRequestModelObject = TRUE), analytic = TRUE, b = mo_b0,
                   R = mo_R, J = mo_J) {
  settings <- c(mo_Q, settings)
  if (analytic) {
    bgw_mle(R, b, calcJ = J, bgw_settings = settings, calcJ_type = "log_probability")
  } else {
    bgw_mle(R, b, bgw_settings = settings)
  }
}

test_that("the fields are correct with analytic derivatives", {
  m <- mo_fit()
  G <- mo_J(mo_b0)
  V <- solve(crossprod(G))
  expect_s3_class(m, "bgw_mle")
  expect_equal(m$code, 20)
  expect_equal(m$iterations, 0)
  expect_identical(m$estimate, mo_b0)
  expect_identical(m$betaStop, mo_b0)
  expect_identical(m$betaStart, mo_b0)
  expect_equal(m$maximum, sum(log(mo_R(mo_b0))))
  expect_equal(m$finalLL, m$maximum)
  expect_equal(m$gradient, unname(-colSums(G)))
  expect_equal(unname(m$varcovBGW), unname(V))
  expect_identical(rownames(m$varcovBGW), names(mo_b0))
  expect_equal(unname(m$seBGW), unname(sqrt(diag(V))))
  expect_equal(unname(m$tstatBGW), unname(mo_b0/sqrt(diag(V))))
  expect_equal(m$hessianMethodUsedBGW, "bhhh")
  e <- eigen(crossprod(G), symmetric = TRUE, only.values = TRUE)$values
  expect_equal(m$vcHessianRecipCondBGW, min(e)/max(e))
  expect_equal(m$vcHessianConditionNumber, m$vcHessianRecipCondBGW)
})

test_that("finite differences agree with analytic derivatives", {
  a <- mo_fit()
  f <- mo_fit(analytic = FALSE)
  expect_equal(f$maximum, a$maximum)
  expect_equal(f$gradient, a$gradient, tolerance = 1e-5)
  expect_equal(f$varcovBGW, a$varcovBGW, tolerance = 1e-4)
  expect_false(f$hasAnalyticGrad)
})

test_that("at the estimates, the object matches a normal run", {
  e <- mo_fit(settings = list())
  m <- mo_fit(b = e$estimate)
  expect_equal(m$maximum, e$maximum, tolerance = 1e-12)
  expect_equal(m$gradient, e$gradient, tolerance = 1e-6)
  expect_equal(m$varcovBGW, e$varcovBGW, tolerance = 1e-10)
  expect_equal(m$seBGW, e$seBGW, tolerance = 1e-10)
  expect_equal(m$vcVecBGW, e$vcVecBGW, tolerance = 1e-10)
})

test_that("a singular BHHH matrix gives no covariance, but still gives the rest", {
  R <- function(b) mo_R(c(b[1], b[2], 0))
  J <- function(b) cbind(mo_J(c(b[1], b[2], 0))[, 1:2], 0)
  m <- mo_fit(R = R, J = J)
  expect_equal(m$code, 20)
  expect_null(m$varcovBGW)
  expect_null(m$seBGW)
  expect_null(m$tstatBGW)
  expect_null(m$vcVecBGW)
  expect_equal(m$hessianMethodUsedBGW, "N/A - BHHH Hessian is singular")
  expect_lte(m$vcHessianRecipCondBGW, .Machine$double.eps)
  expect_identical(m$estimate, mo_b0)
})

test_that("vcHessianMethod = none skips the covariance", {
  m <- mo_fit(settings = list(userRequestModelObject = TRUE, vcHessianMethod = "none"))
  expect_null(m$varcovBGW)
  expect_equal(m$hessianMethodUsedBGW, "None")
  expect_equal(m$vcHessianRecipCondBGW, -1)
  expect_length(m$gradient, 3)
})

test_that("maxIterations = 0 gives the same object, apart from its description", {
  a <- mo_fit()
  z <- mo_fit(settings = list(maxIterations = 0))
  expect_match(a$message, "userRequestModelObject = TRUE", fixed = TRUE)
  expect_match(z$message, "maxIterations = 0", fixed = TRUE)
  skip <- c("message", "bgw_settings")
  expect_identical(z[setdiff(names(z), skip)], a[setdiff(names(a), skip)])
})

test_that("silent prints nothing", {
  expect_silent(mo_fit())
  expect_silent(mo_fit(settings = list(maxIterations = 0)))
})

test_that("finite-difference Hessians are rejected with a clear message", {
  for (method in c("finiteDifferences", "fdFunction")) {
    expect_error(mo_fit(settings = list(userRequestModelObject = TRUE, vcHessianMethod = method)),
                 "requires vcHessianMethod")
    expect_error(mo_fit(settings = list(maxIterations = 0, vcHessianMethod = method)),
                 "requires vcHessianMethod")
  }
})

test_that("individual-level gradients are returned when requested", {
  m <- mo_fit(settings = list(userRequestModelObject = TRUE, userRequestIndLevelGradients = TRUE))
  expect_equal(unname(m$indLevelGradients), unname(mo_J(mo_b0)))
})

test_that("an existing iterations file is left alone", {
  d <- tempfile("bgwmo")
  dir.create(d)
  f <- file.path(d, "m_iterations.csv")
  writeLines("keep", f)
  m <- mo_fit(settings = list(userRequestModelObject = TRUE, writeIter = TRUE, modelName = "m",
                              outputDirectory = d))
  expect_identical(readLines(f), "keep")
})

test_that("a normal estimation still reports the renamed fields", {
  e <- mo_fit(settings = list())
  expect_true(e$code %in% 3:6)
  expect_equal(e$vcHessianConditionNumber, e$vcHessianRecipCondBGW)
  expect_gt(e$vcHessianRecipCondBGW, 0)
  expect_equal(e$hessianMethodUsedBGW, "bhhh")
  expect_length(e$vcVecBGW, 6)
})
