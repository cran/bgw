# Constants and tables shared by the bgw R functions. Internal; none of this is exported.

# Positions of named values in the BGW Fortran work vectors iv (integer) and v (double), using the
# Fortran names. For example, iv[bgw_iv_index$mxiter] holds the iteration limit. Each name and value
# matches a PARAMETER statement in src/*.f90, except those tagged [R only]. Entries are in numeric
# order, to ease comparison with the Fortran. Some positions are only used by the optimizer
# variants that bgw does not currently call (alg = 2), and are kept because bgw_mle_setup.R
# documents the full Fortran layout.
#
# The tag on each line gives the role of the entry:
#   [set]     a limit, tolerance or flag; R assigns it and the Fortran reads it
#   [layout]  a storage position or size in iv or v, computed in bgw_mle_setup.R
#   [result]  written by the Fortran as it runs; R at most initializes or reads it
#   [both]    updated by both R and the Fortran during estimation
#   [R only]  not in src/*.f90; from the original DGLG or DITSUM, which R code replaces
#   [unused]  not referenced by the R code
# Descriptions come from the Fortran comments, from the parameter documentation published with
# TOMS 611 (SUMSL), or from how the code uses the entry.
bgw_iv_index <- list(
  toobig   = 2,    # [both]   nonzero if the step was too big (for example, it caused overflow)
  ivneed   = 3,    # [layout] extra length of iv requested
  vneed    = 4,    # [layout] extra length of v requested
  nfcall   = 6,    # [result] invocation count for the objective function
  nfgcal   = 7,    # [both]   value of nfcall at the step with the biggest function reduction
  covprt   = 14,   # [set]    odd = print the covariance; >= 2 = print regression diagnostics
  covreq   = 15,   # [set]    vcHessianMethod (3 = BHHH; 2, -2 = finite differences)
  dtype    = 16,   # [set]    0 = do not update the scale vector d
  mxfcal   = 17,   # [set]    function evaluation limit (Fortran default 200; bgw default 300)
  mxiter   = 18,   # [set]    iteration limit (Fortran default 150; bgw default 200)
  outlev   = 19,   # [set]    a summary is printed every abs(outlev) iterations
  parprt   = 20,   # [set]    1 = print nondefault v values
  prunit   = 21,   # [unused] output unit for Fortran printing; 0 = no printing
  solprt   = 22,   # [set]    1 = print the returned x, gradient and scale vector d
  statpr   = 23,   # [set]    1 = print summary statistics on return
  x0prt    = 24,   # [set]    1 = print the initial x and scale vector d
  inith    = 25,   # [set]    same position as inits
  inits    = 25,   # [set]    same position as inith
  covmat   = 26,   # [result] start of the covariance in v; 0 = none, negative = failed
  d        = 27,   # [layout] [R only] start of the scale vector d in v
  g        = 28,   # [layout] start of the gradient in v
  ngcall   = 30,   # [both]   gradient evaluation count
  niter    = 31,   # [result] iteration count
  mode     = 35,   # [result] -1 at the start; 0 = normal; > 0 = FD Hessian for the covariance
  dig      = 37,   # [layout] start of D**-1 * g in v (its norm is v(dgnorm))
  prntit   = 39,   # [R only] iterations since the last printed summary
  step     = 40,   # [layout] start of the current step in v
  stlstg   = 41,   # [layout] start of the saved last good step in v
  lmat     = 42,   # [layout] start of the Cholesky factor L of the Hessian model, in v
  x0       = 43,   # [layout] start of x as of the beginning of the iteration, in v
  lastiv   = 44,   # [layout] last position used in iv
  lastv    = 45,   # [layout] last position used in v
  nextiv   = 46,   # [layout] next free position in iv
  nextv    = 47,   # [layout] next free position in v
  y        = 48,   # [layout] start of the vector y in v (scratch; used to update S)
  parsav   = 49,   # [layout] start of the saved copy of the input parameters in v
  nvdflt   = 50,   # [set]    number of v defaults (must be 32 for alg = 1)
  algsav   = 51,   # [set]    algorithm choice (alg)
  nfcov    = 52,   # [result] function evaluations used for the covariance
  ngcov    = 53,   # [both]   gradient evaluations used for the covariance
  h        = 56,   # [layout] Hessian storage, also used for the covariance; its sign is a flag
  rdreq    = 57,   # [set]    whether to compute the covariance and regression diagnostics
  perm     = 58,   # [layout] start of the variable part of iv, after the fixed block
  jtol     = 59,   # [layout] start of dtol (p) and d0 (p) in v, used when updating d
  vsave    = 60,   # [layout] start of the saved copy of v(1:9), kept while trying another model
  r        = 61,   # [layout] [R only] start of the residuals in v
  s        = 62,   # [layout] start of S, the second-order part of the Hessian model, in v
  sused    = 64,   # [result] Hessian model behind the step; printed as MODEL in the summary
  w        = 65,   # [layout] start of work storage in v
  jcn      = 66,   # [layout] start of a p-length work vector, used for regression diagnostics
  j        = 70,   # [layout] [R only] start of the Jacobian (dr) in v
  hc       = 71,   # [layout] start of the Hessian in v when used; bgw sets 0
  fdh      = 74,   # [unused] position of the finite-difference Hessian; set by the Fortran
  ierr     = 75,   # [result] return code from the QR factorization; 0 = R has full rank
  ipivot   = 76,   # [layout] start of the QR pivot order in iv
  qtr      = 77,   # [layout] start of Q**T * r in v
  rmat     = 78,   # [layout] start of the QR factor R in v; 0 = none; its sign is a flag
  qrtyp    = 80,   # [set]    set to 1 by default; not otherwise used
  regd0    = 82,   # [layout] [R only] start of the regression diagnostics (rd) in v
  dradpr   = 101   # [unused] set only when alg > 2; not otherwise used
)

# Positions in v; the tags are the same as above.
bgw_v_index <- list(
  dstnrm   = 2,    # [result] 2-norm of D*STEP
  stppar   = 5,    # [result] Marquardt parameter; 0 = full Newton step
  nreduc   = 6,    # [result] function reduction predicted by the quadratic model, Newton step
  preduc   = 7,    # [result] function reduction predicted by the quadratic model, current step
  f        = 10,   # [result] objective function value at x
  fdif     = 11,   # [result] function reduction v(f0) - v(f)
  f0       = 13,   # [result] objective function value at the start of the iteration
  reldx    = 17,   # [result] scaled relative change in x caused by the step
  epslon   = 19,   # [set]    maximum relative error allowed for psi(step)
  phmnfc   = 20,   # [set]    with phmxfc, tolerance for accepting the step
  phmxfc   = 21,   # [set]    see phmnfc
  decfac   = 22,   # [set]    factor to decrease the radius when iv(toobig) is nonzero
  incfac   = 23,   # [set]    minimum factor to increase the radius
  rdfcmn   = 24,   # [set]    minimum value for v(radfac)
  rdfcmx   = 25,   # [set]    maximum value for v(radfac)
  tuner1   = 26,   # [set]    decides if the function reduction was much less than expected
  tuner2   = 27,   # [set]    decides if the function reduction was enough to accept the step
  tuner3   = 28,   # [set]    decides if the radius should be increased
  tuner4   = 29,   # [set]    gradient test on whether to increase the radius (with v(dgnorm))
  tuner5   = 30,   # [set]    gradient test on whether to increase the radius (with v(gtstep))
  afctol   = 31,   # [set]    absolute function convergence tolerance
  rfctol   = 32,   # [set]    relative function convergence tolerance
  xctol    = 33,   # [set]    x-convergence criterion
  xftol    = 34,   # [set]    false convergence tolerance
  lmax0    = 35,   # [set]    maximum 2-norm of D*STEP for the very first step
  lmaxs    = 36,   # [set]    maximum reasonable step size (and initial step bound)
  sctol    = 37,   # [set]    singular convergence tolerance
  dinit    = 38,   # [set]    if nonnegative, every element of d starts at this value
  dtinit   = 39,   # [set]    if positive, initial value of every element of dtol (see dfac)
  d0init   = 40,   # [set]    if positive, initial value of every element of d0 (see dfac)
  dfac     = 41,   # [set]    used when updating d: d1(i) = max(sqrt(abs(h(i,i))), dfac*d(i))
  dltfdc   = 42,   # [set]    same position as eta0
  eta0     = 42,   # [set]    same position as dltfdc
  bias     = 43,   # [set]    same position as dltfdj
  dltfdj   = 43,   # [set]    same position as bias
  delta0   = 44,   # [set]    relative step for the finite-difference Hessian
  fuzz     = 45,   # [set]    used in deciding whether to change Hessian models
  rlimit   = 46,   # [set]    given a default; not otherwise used
  cosmin   = 47,   # [set]    smallest acceptable cosine in the update of S
  huberc   = 48,   # [set]    given a default; not otherwise used
  rsptol   = 49,   # [set]    relative cutoff on rho''; below it, seminormal equations replace QR
  sigmin   = 50,   # [set]    given a default; not otherwise used
  rcond    = 53    # [result] its square bounds the reciprocal condition number of vcHessian
)

# Text for BGW return codes, indexed by iv[1]. bgw_mle stores trimws() of the annotated version in
# model$message; bgw_itsum prints the plain version.
#
# iv[3] to iv[6] are "good convergence"; iv[7] and iv[8] are "bad convergence"; iv[9] to iv[11] are
# "premature stopping of search".
#
# The following occurs on entry to ditsum if iv[1] > 62, and means a problem occurred at the START.
# If iv[1] > 61, iv1 = iv[1] - 51, so:
#   iv[1] = 63 => iv1 = 12; iv[1] = 64 => iv1 = 13; iv[1] = 65 => iv1 = 14; iv[1] = 66 => iv1 = 15
# iv1 >= 16 (i.e., iv[1] >= 67) simply writes out iv[1]. These values would typically occur when the
# user overwrites a default value with an illegal value. iv[1] = 67 => iv1 = 16 can be set two
# different ways in dparck; message 16 should cover it.
bgw_messages <- function(annotated = FALSE) {
  tbl <- list(
    "       ***** Evaluate function ***** \n",
    "       ***** Evaluate gradient ***** \n",
    "       ***** X-convergence ***** \n",
    "       ***** Relative function convergence ***** \n",
    "       ***** X- and relative function convergence ***** \n",
    "       ***** Absolute function convergence ***** \n",
    "       ***** Singular convergence ***** \n",
    "       ***** False convergence ***** \n",
    "       ***** Function evaluation limit ***** \n",
    "       ***** Iteration limit ***** \n",
    "       ***** stopx ***** \n",
    "       ***** Initial f(x) cannot be computed ***** \n",
    "       ***** Bad parameters to assess ***** \n",
    "       ***** Gradient could not be computed ***** \n",
    "       ***** Inconsistent dimensions (Could be N, ND, P, or LIV. )***** \n",
    " ***** Internal error specifying alg in divset and/or dparck. It should be 1. ***** \n"
  )
  if (annotated) {
    # Entry 6 (absolute function convergence) is deliberately not annotated: its interpretation is
    # more ambiguous than the other convergence codes, and it almost never occurs.
    note <- c("3" = "favorable", "4" = "favorable", "5" = "favorable", "7" = "unfavorable",
              "8" = "unfavorable", "9" = "unfavorable", "10" = "unfavorable")
    for (i in names(note)) {
      k <- as.integer(i)
      tbl[[k]] <- sub(" ***** \n", paste0(" (", note[[i]], ") ***** \n"), tbl[[k]], fixed = TRUE)
    }
  }
  tbl
}
