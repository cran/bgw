# Tests of the files written at each iteration: writeIter, writeIterGradient and writeIterScaleVec.

set.seed(3)
wf_n  <- 400
wf_x1 <- rnorm(wf_n)
wf_x2 <- rnorm(wf_n)
wf_x3 <- rnorm(wf_n)
wf_y  <- rbinom(wf_n, 1, plogis(0.8*wf_x1 - 0.5*wf_x2 + 0.2*wf_x3))
wf_P  <- function(b) {
  b <- unname(b)
  plogis(b[1]*wf_x1 + b[2]*wf_x2 + b[3]*wf_x3)
}
wf_R  <- function(b) ifelse(wf_y == 1, wf_P(b), 1 - wf_P(b))
wf_J  <- function(b) cbind(wf_x1, wf_x2, wf_x3)*(wf_y - wf_P(b))
wf_b0 <- c(a = 0.1, b = 0.2, c = -0.1)
wf_all <- list(writeIter = TRUE, writeIterGradient = TRUE, writeIterScaleVec = TRUE,
               modelName = "m")

wf_fit <- function(dir, settings = wf_all, analytic = TRUE, b = wf_b0) {
  settings <- c(list(silent = TRUE, outputDirectory = dir), settings)
  if (analytic) {
    bgw_mle(wf_R, b, calcJ = wf_J, bgw_settings = settings, calcJ_type = "log_probability")
  } else {
    bgw_mle(wf_R, b, bgw_settings = settings)
  }
}
wf_dir  <- function() {
  d <- tempfile("bgwwf")
  dir.create(d)
  d
}
wf_read <- function(dir, suffix) {
  utils::read.csv(file.path(dir, paste0("m_", suffix, ".csv")), check.names = FALSE)
}

test_that("all three files share one layout, and one row per iteration", {
  d <- wf_dir()
  m <- wf_fit(d)
  it <- wf_read(d, "iterations")
  gr <- wf_read(d, "gradients")
  sv <- wf_read(d, "scaleVecs")
  cols <- c("a", "b", "c", "logLike", "Iteration", "ItType")
  expect_identical(names(it), cols)
  expect_identical(names(gr), cols)
  expect_identical(names(sv), cols)
  expect_true(all(it$ItType == "estimate"))
  expect_true(all(gr$ItType == "gradient"))
  expect_true(all(sv$ItType == "scaleVec"))
  expect_identical(gr$Iteration, it$Iteration)
  expect_identical(sv$Iteration, it$Iteration)
  expect_equal(gr$logLike, it$logLike)
  expect_false(any(duplicated(it$Iteration)))
  expect_equal(it$Iteration[1], 0)
  expect_equal(max(it$Iteration), m$iterations)
  expect_silent(stacked <- rbind(it, gr, sv))
  expect_equal(nrow(stacked), 3*nrow(it))
})

test_that("each row belongs to its iterate, and the last agrees with the model", {
  d <- wf_dir()
  m <- wf_fit(d)
  it <- wf_read(d, "iterations")
  gr <- wf_read(d, "gradients")
  sv <- wf_read(d, "scaleVecs")
  for (i in seq_len(nrow(it))) {
    b <- as.numeric(it[i, c("a", "b", "c")])
    expect_equal(it$logLike[i], sum(log(wf_R(b))), tolerance = 1e-10)
    expect_equal(as.numeric(gr[i, c("a", "b", "c")]), unname(colSums(wf_J(b))), tolerance = 1e-6)
  }
  last <- nrow(it)
  expect_equal(as.numeric(it[last, c("a", "b", "c")]), unname(m$estimate), tolerance = 1e-10)
  expect_equal(as.numeric(gr[last, c("a", "b", "c")]), -m$gradient, tolerance = 1e-6)
  expect_equal(as.numeric(sv[last, c("a", "b", "c")]), m$scaleVec, tolerance = 1e-10)
})

test_that("finite-difference runs write the same layout", {
  d <- wf_dir()
  # No covariance matrix: with one, the Fortran asks for the Jacobian again at the final point, and
  # its finite-difference steps then use the updated scale vector, so the returned gradient could
  # differ slightly from the last one written.
  m <- wf_fit(d, settings = c(wf_all, vcHessianMethod = "none"), analytic = FALSE)
  gr <- wf_read(d, "gradients")
  expect_equal(max(gr$Iteration), m$iterations)
  expect_equal(as.numeric(gr[nrow(gr), c("a", "b", "c")]), -m$gradient, tolerance = 1e-6)
})

test_that("unnamed parameters get generated column names", {
  d <- wf_dir()
  m <- wf_fit(d, b = unname(wf_b0))
  expect_identical(names(wf_read(d, "iterations"))[1:3], c("beta1", "beta2", "beta3"))
  expect_identical(names(wf_read(d, "gradients"))[1:3], c("beta1", "beta2", "beta3"))
})

test_that("writeIterMode is followed by all three files", {
  d <- wf_dir()
  wf_fit(d)
  n1 <- nrow(wf_read(d, "scaleVecs"))
  wf_fit(d)
  expect_equal(nrow(wf_read(d, "scaleVecs")), n1)
  wf_fit(d, settings = c(wf_all, writeIterMode = "append"))
  expect_equal(nrow(wf_read(d, "scaleVecs")), 2*n1)
  expect_equal(nrow(wf_read(d, "gradients")), 2*n1)
  expect_equal(nrow(wf_read(d, "iterations")), 2*n1)
})

test_that("appending to an old-format file keeps its layout", {
  d <- wf_dir()
  f <- file.path(d, "m_iterations.csv")
  writeLines(c('"a","b","c","logLike"', "1,2,3,4"), f)
  wf_fit(d, settings = list(writeIter = TRUE, modelName = "m", writeIterMode = "append"))
  it <- utils::read.csv(f, check.names = FALSE)
  expect_identical(names(it), c("a", "b", "c", "logLike"))
  expect_gt(nrow(it), 1)
  expect_true(all(is.finite(it$logLike)))
})

test_that("nothing is written unless requested", {
  d <- wf_dir()
  wf_fit(d, settings = list(modelName = "m"))
  expect_length(list.files(d), 0)
})

test_that("a finite-difference Hessian does not add rows", {
  d <- wf_dir()
  m <- wf_fit(d, settings = c(wf_all, vcHessianMethod = "finiteDifferences"))
  it <- wf_read(d, "iterations")
  expect_false(any(duplicated(it$Iteration)))
  expect_equal(max(it$Iteration), m$iterations)
})

test_that("a model object at betaStart writes no files", {
  d <- wf_dir()
  wf_fit(d, settings = c(wf_all, userRequestModelObject = TRUE))
  expect_length(list.files(d), 0)
})

test_that("printing does not fail when the gradient cannot be computed", {
  R <- function(b) if (isTRUE(all.equal(unname(b), unname(wf_b0)))) wf_R(b) else rep(NA_real_, wf_n)
  expect_no_error(utils::capture.output(m <- bgw_mle(R, wf_b0)))
  expect_s3_class(m, "bgw_mle")
  expect_equal(m$code, 65)
  expect_equal(m$maximum, sum(log(wf_R(wf_b0))))
  expect_equal(m$finalLL, m$maximum)
  expect_match(m$message, "Gradient could not be computed", fixed = TRUE)
})
