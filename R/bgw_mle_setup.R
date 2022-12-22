#' @title bgw_mle_setup
#'
#' @description Sets up R-level storage for bgw_mle.
#'
#' This function replaces multiple Fortran subroutines from the BGW
#' Fortran code due at least in part to the prohibition against using
#' Fortran write statements in R packages. The current design produces
#' two vectors (iv_r and v_r) that mirror the main vectors required by the
#' Fortran code. For the moment, the idea is to create named vectors to
#' facilitate coding in bgw_mle. It may be that these could be deleted
#' (or overwritten by an as.numeric() conversion) prior to the main Fortran
#' calls.
#'
#' @param p  Number of parameters (components of x) being estimated. (Determined in bgw_mle from size of 'start' vector.
#' @param n Number of model residuals (in vector r). (Determined in bgw_mle by the size of the output vector from CalcR.)
#' @param hasAnalyticModelDeriv Logical. TRUE if CalcRJ has been provided. FALSE means that bgw_mle must employ finite-difference gradients (which has implications for storage allocation).
#' @param control List of bgw_mle control parameters (optional). If not provided, BGW default parameters will be used. If provided, default parameters will be overwritten by those corresponding parameters provided by the caller (but these must also be checked).
#'
#' @return iv and v vectors used by BGW Fortran.
#         \itemize{
#           \item \strong{\code{"iv_r"}}: Internal R version of BGW Fortran vector iv
#           \item \strong{\code{"v_r"}}:  Internal R version of BGW Fortran vector v
#         }
#' @export
##------------------------------------------------------------------------
# October 14, 2022
bgw_mle_setup <- function(p, n, hasAnalyticModelDeriv, control=NULL) {

  # DSB remark: The first part of this setup code replicates the functionality of the
  # Fortran subroutine divset in the BGW code. We maintain detailed documentation
  # here for future reference. The call is:
  # CALL divset_f(alg, iv, liv, lv, v)

  # DSB remark: BGW Fortran code accesses a system-specific function to identify the
  # correct Fortran output unit number. However, we are precluded from using this
  # in an R package.
  # IF (prunit <= liv) iv_r[prunit] = i7mdcn(1)

  # DSB remarks: The BGW Fortran code is general, and can be used for different
  # types of problems. alg = 1 for statistical estimation, and alg = 2
  # for general optimization. divset stores alg so it can be checked later on
  # when subroutine dparck is called. However, in this package we assume alg = 1 througout.
  # IF (algsav <= liv) iv_r[algsav] = alg
  # It also performs error checking.
  # IF (alg < 1 .OR. alg > 4) GO TO 40
  # [Additional error checking on sizes of iv and v omitted.]

  # ############################# #
  #### Initial error checking ####
  # ############################# #
  if (!is.numeric(p)) stop(" Non-numeric p passed to bgw_mle_setup.")
  if (p <= 0) stop(" Value of p passed to bgw_mle_setup must be greater than zero.")
  if (!is.numeric(n)) stop(" Non-numeric n passed to bgw_mle_setup.")
  if (n <= 0) stop(" Value of n passed to bgw_mle_setup must be greater than zero.")
  if (!is.logical(hasAnalyticModelDeriv)) stop(" Parameter hasAnalyticModelDeriv must be logical.")

  alg <- 1

  # ############################# #
  #### Initialize iv and v     ####
  # ############################# #
  # DSB remark
  # The current documentation in BGW Fortran indicates that the length of
  # iv should be "at least 90 + p."  We have not yet seen a reason why
  # iv might need to be longer, but if and when it occurs we will document it.
  liv <- 90 + p
  iv_r <- rep(0,liv)

  # Similarly, a size for v is suggested.  As noted above, the size depends
  # on the availability of CalcRJ.
  # DSB remark: BGW Fortran code can accommodate parameter vectors that include
  # nuisance parameters, where nuisance parameters are treated differently than
  # non-nuisance parameters. We do not use this functionality.
  # ps = number of non-nuisance parameters = p.
  ps <- p
  # Compute lv.
  if (hasAnalyticModelDeriv) {
    lv <- 105 + p*(3*p + 16) + 2*n + 4*ps + n*(p + 1 + (p-ps+1)*(p-ps+2)/2)}
  else {
    lv <- 105 + p*(3*p + 16) + 2*n + 4*ps + n*(p + 3 + (p-ps+1)*(p-ps+2)/2)
  }
  v_r <- rep(0,lv)

  # ################################ #
  #### Supply default values to v ####
  # ################################ #

  # DSB remark. In BGW Fortran code, default v values are assigned by the
  # following call:  CALL dv7dfl(alg, lv, v).
  # As noted above, alg = 1.
  #
  # Many of these values are technical in nature, relating to, e.g.,
  # convergence criteria for the optimization algorithm, parameters to use
  # when computing finite difference derivatives, etc.
  # All of these are a function of the machine precision for the system being
  # used. The key quantity is "machine epsilon" (machep) which denotes the
  # unit roundoff on the machine. More formally,
  #  unit roundoff = the smallest positive number machep such that
  #      1 + machep > 1 AND 1 - machep < 1.
  #
  # Other quantities include
  # Fortran has built-in intrinsic functions for this (and other quantities).
  # Specifically:
  #     The smallest positive eta such that -eta exists. ("eta," or "tiny")
  #     The largest machine number big such that -big exists. ("big").
  #
  # We developed a revised version of a double precision Fortran function in BGW
  # called dr7mdc. Even with this revised setup routine, this function is still
  # needed in the remainder of the Fortran code we are using.
  # It computes 6 different constants indexed by an input k.  It
  # would be straightforward to write a wrapper for this function and use it
  # to compute the default v values. However, for now we will use built-in
  # R functions:

  machep <- .Machine$double.eps
  # print(machep)
  big    <- .Machine$double.xmax
  eta    <- .Machine$double.xmin
  #
  # Note: We can do numeric testing on these very shortly.
  # In what follows, we include information on the v
  # subscript names for purposes of clarity.
  #
  afctol <- 31
  bias <- 43
  cosmin <- 47
  decfac <- 22
  delta0 <- 44
  dfac <- 41
  dinit <- 38
  dltfdc <- 42
  dltfdj <- 43
  dtinit <- 39
  d0init <- 40
  epslon <- 19
  eta0 <- 42
  fuzz <- 45
  huberc <- 48
  incfac <- 23
  lmax0 <- 35
  lmaxs <- 36
  phmnfc <- 20
  phmxfc <- 21
  rdfcmn <- 24
  rdfcmx <- 25
  rfctol <- 32
  rlimit <- 46
  rsptol <- 49
  sctol <- 37
  sigmin <- 50
  tuner1 <- 26
  tuner2 <- 27
  tuner3 <- 28
  tuner4 <- 29
  tuner5 <- 30
  xctol <- 33
  xftol <- 34

  # machep <- dr7mdc(3)
  v_r[afctol] <- 1.e-20
  if (machep > 1e-10) v_r[afctol] <- machep^2
  v_r[decfac] <- 0.5e0
  sqteps <- sqrt(machep)
  v_r[dfac] <- 0.6e0
  v_r[dtinit] <- 1.e-6
  mepcrt <- machep**(1./3.)
  v_r[d0init] <- 1.e+0
  v_r[epslon] <- 0.1e0
  v_r[incfac] <- 2.e+0
  v_r[lmax0] <- 1.e+0
  v_r[lmaxs] <- 1.e+0
  v_r[phmnfc] <- -0.1e0
  v_r[phmxfc] <- 0.1e0
  v_r[rdfcmn] <- 0.1e0
  v_r[rdfcmx] <- 4.e+0
  v_r[rfctol] <- max(1.e-10, mepcrt**2)
  v_r[sctol] <- v_r[rfctol]
  v_r[tuner1] <- 0.1e0
  v_r[tuner2] <- 1.e-4
  v_r[tuner3] <- 0.75e0
  v_r[tuner4] <- 0.5e0
  v_r[tuner5] <- 0.75e0
  v_r[xctol] <- sqteps
  v_r[xftol] <- 1.e+2 * machep

 if (alg == 1) {

 #  ***  REGRESSION  VALUES
  v_r[cosmin] <- max(1.e-6, 1.e+2 * machep)
  v_r[dinit]  <- 0.e+0
  v_r[delta0] <- sqteps
  v_r[dltfdc] <- mepcrt
  v_r[dltfdj] <- sqteps
  v_r[fuzz]   <- 1.5e0
  v_r[huberc] <- 0.7e0
  # v_r[rlimit] <- dr7mdc(5)
  v_r[rlimit] <- sqrt(big/256.e+0)*16.e+0
  v_r[rsptol] <- 1.e-3
  v_r[sigmin] <- 1.e-4
  }

  # ***  GENERAL OPTIMIZATION VALUES
  if (alg == 2) {
  v_r[bias] <- 0.8e0
  v_r[dinit] <- -1.e0
  v_r[eta0] <- 1.0e3 * machep
  }

 # Next, we move on to iv_r values
 # iv subscripts
 ivneed <- 3
 vneed <- 4
 nvsave <- 9
 covprt <- 14
 covreq <- 15
 dtype <- 16
 mxfcal <- 17
 mxiter <- 18
 outlev <- 19
 parprt <- 20
 prunit <- 21
 solprt <- 22
 statpr <- 23
 x0prt <- 24
 inith <- 25
 inits <- 25
 lmat <- 42
 lastiv <- 44
 lastv <- 45
 parsav <- 49
 nvdflt <- 50
 algsav <- 51
 nfcov <- 52
 ngcov <- 53
 rdreq <- 57
 perm <- 58
 vsave <- 60
 hc <- 71
 ierr <- 75
 ipivot <- 76
 rmat <- 78
 qrtyp <- 80
 dradpr <- 101
 #
 if (alg == 1) {
   miv <- 82
   mv <- 98
 } else {
   miv <- 59
   mv <- 71
 }

 iv_r[algsav] = alg
 iv_r[ivneed] <- 0
 iv_r[lastiv] <- miv
 iv_r[lastv]  <- mv
 iv_r[lmat]   <- mv + 1
 iv_r[mxfcal] <- 200
 iv_r[mxiter] <- 150
 iv_r[outlev] <- 1
 iv_r[parprt] <- 1
 iv_r[perm]   <- miv + 1
 iv_r[solprt] <- 1
 iv_r[statpr] <- 1
 iv_r[vneed]  <- 0
 iv_r[x0prt]  <- 1

 alg <- 1
 if (alg == 1){
   #  ***  REGRESSION  VALUES
   # iv_r[covprt] <- 3
   iv_r[covprt] <- 3
   iv_r[covreq] <- 1
   iv_r[dtype]  <- 1
   iv_r[hc]    <- 0
   iv_r[ierr]   <- 0
   iv_r[inits]  <- 0
   iv_r[ipivot] <- 0
   iv_r[nvdflt] <- 32
   iv_r[vsave]  <- 58
   iv_r[parsav] <- iv_r[vsave] + nvsave
   iv_r[qrtyp] <- 1
   iv_r[rdreq] <- 3
   iv_r[rmat] <- 0
 }

 if (alg == 2){

   #  ***  GENERAL OPTIMIZATION VALUES
   iv_r[dtype] <- 0
   iv_r[inith] <- 1
   iv_r[nfcov] <- 0
   iv_r[ngcov] <- 0
   iv_r[nvdflt] <- 25
   iv_r[parsav] <- 47
 }
 # Never used (only here for 'future' code with contraints)
 if (alg > 2) iv_r[vsave] <- iv_r[vsave] + 3
 if (alg > 2) iv_r[parsav] <- 61

 # END of divset-based code

 # Note:  This may be the place to implement elements of control that
 # involve changing default parameters.

 # Location" is in DGLF/DGLFG with iv(1) = 12
 # iv(1) = 13 in anticipation of call to DRGLG to do "storage allocation"
 # But iv1 = 12 (still)
 ifac <- (p-ps+2)*(p-ps+1)/2
 if (hasAnalyticModelDeriv) {
   iv_r[vneed] <- iv_r[vneed] + p + n*(p + 3 + ifac)
 } else {
   iv_r[covreq] <- -abs(iv_r[covreq])
   if (iv_r[covreq] == 0 && iv_r[rdreq] > 0) iv_r[covreq] <- -1
   iv_r[vneed] <- iv_r[vneed] + p + n*(p + 1 + ifac)
 }

 # Next is pseudo-call to DRGLG with iv(1) = 13
 lh <- p*(p+1)/2

 # nd is leading dimension of dr = derivative of r wrt x (so, nd x p). Must be >= ps.
 nd <- ps
 # nn is dimension of r, and leading dimension of rd
 # Not sure yet if we really need this here.
 nn <- n
 iv_r[ivneed] <- iv_r[ivneed] + p
 iv_r[vneed]  <- iv_r[vneed] + p*(p+13)/2 + 2*n + 4*ps

 # Next is an adjustment to iv[perm]
 # After divset, iv[perm] = 83 when alg = 1
 # parameter in DRGLG:
 # xnoti <- 90
 # i = xnoti + 1 = 91
 # if (iv(perm) < i) iv(perm) = i
 iv_r[perm] <- 91

 # Next is pseudo-call to DG7LIT (from DRGLG) with iv(1) = 13
 iv_r[vneed] <- iv_r[vneed] + p*(3*p + 19)/2 + 7

 # Next is pseudo-call to DPARCK with iv(1) = 13
 miv1 <- miv
 # IF (perm <= liv) miv1 = MAX0(miv1, iv(perm) - 1)
 # Note: This is strange, but unlikely to be a mistake.
 # perm <- 58 is hardwired
 # liv = 90+9 miv is 82 for alg = 1 (but would be 59 for alg=2)
 # At this stage iv(perm) = 91 so:
 miv1 <- max(miv,iv_r[perm]-1)
 # This will be 90, no matter what.
 miv2 <- 0
 # IF (ivneed <= liv) miv2 = miv1 + MAX0(iv(ivneed), 0)
 # ivneed = 3, which is clearly less than liv
 # At this stage, iv(ivneed) = p
 miv2 <- miv1 + max(iv_r[ivneed],0)
 # Clearly, miv2 = 90 + p
 # IF (lastiv <= liv) iv(lastiv) = miv2
 # lastiv = 44 which is clearly less than liv = 90 + 9
 iv_r[lastiv] <- miv2
 # At this point iv[lastiv] = 90 + P

 # The following does not occur
 # IF (liv < miv1) GO TO 300

 # For some reason iv(ivneed) is now being reset to 0
 # iv(ivneed) = 0
 iv_r[ivneed] <- 0

 # Now we look at iv(lastv)
 # From divset iv(lastv) = mv = 98 when alg = 1
 # At this stage iv(vneed) has been accumulated, and,
 # iv(lmat) = mv + 1 = 99 for alg = 1
 # iv(lastv) = MAX0(iv(vneed), 0) + iv(lmat) - 1
 iv_r[lastv] <- max(iv_r[vneed],0) + iv_r[lmat] - 1
 # We need to numerically check on the above value.

 # We are re-setting iv(vneed) = 0
 # iv(vneed) = 0
 iv_r[vneed] <- 0

 # Recall: At this point both liv and miv2 = 90 + p
 # IF (liv < miv2) GO TO 300

 # This is the real acid test.  It appears as though
 # all of the calculations leading to iv(lastv) are those that
 # give rise to the original lv calculation. So the following should
 # not be a problem.
 # IF (lv < iv(lastv)) GO TO 320
 # 30  CONTINUE

 # IF (n >= 1) GO TO 50 [no problem]
 # 50 IF (iv1 /= 14) iv(nextiv) = iv(perm)
 #    IF (iv1 /= 14) iv(nextv) = iv(lmat)
 #    IF (iv1 == 13) GO TO 999
 # At this stage iv1 = 13. So, the above two statements are executed.
 # We haven't used this yet... Recall iv(perm) = 91
 nextiv <- 46
 iv_r[nextiv] <- iv_r[perm]
 # So, iv(nextiv) = 91.  This would be the starting point of a p-vector (if needed)
 # iv(nextv) is being moved back to 99, which is apparently the starting point for
 # storage used by BGW for problem-related data (not infrastructure).
 nextv <- 47
 iv_r[nextv] <- iv_r[lmat]
 #
 # This location is the pseudo-return from dparck back to dg7lit (with iv(1) = 13).
 # [Back inside dg7lit.]
 # This takes us to storage allocation inside dg7lit
 # !  ***  STORAGE ALLOCATION  ***
 #
 #   10   pp1o2 = p * (p + 1) / 2
 pp1o2 <- p*(p+1)/2
 # iv(s) = iv(lmat) + pp1o2
 s <- 62
 iv_r[s] <- iv_r[lmat] + pp1o2

 # iv(x0) = iv(s) + pp1o2
 x0 <- 43
 iv_r[x0] <- iv_r[s] + pp1o2

 # iv(step) = iv(x0) + p
 step <- 40
 iv_r[step] <- iv_r[x0] + p

 # iv(stlstg) = iv(step) + p
 stlstg <- 41
 iv_r[stlstg] <- iv_r[step] + p

 # iv(dig) = iv(stlstg) + p
 dig <- 37
 iv_r[dig] <- iv_r[stlstg] + p

 # iv(w) = iv(dig) + p
 w <- 65
 iv_r[w] <- iv_r[dig] + p

 # iv(h) = iv(w) + 4*p + 7
 h <- 56
 iv_r[h] <- iv_r[w] + 4*p + 7

 # iv(nextv) = iv(h) + pp1o2
 iv_r[nextv] <- iv_r[h] + pp1o2

 # The following is the return from dg7lit to drglg (with iv(1) = 14).
 # IF (iv(1) /= 13) GO TO 20
 # iv(1) = 14
 # GO TO 999
 # [Back inside drglg, just after line 20]
 #
 # Final bit of storage allocation inside drglg
 # !  ***  STORAGE ALLOCATION  ***
 # iv(ipivot) = iv(nextiv)
 # At this stage iv[nextiv] = 91
 iv_r[ipivot] <- iv_r[nextiv]

 # iv(nextiv) = iv(ipivot) + p
 # iv[nextiv] becomes 90 + p + 1
 #   = what the next one WOULD be if for some reason we needed
 #     additional storage in iv (but, we don't).
 iv_r[nextiv] = iv_r[ipivot] + p

 # iv(y) = iv(nextv)
 y <- 48
 iv_r[y] <- iv_r[nextv]

 # iv(g) = iv(y) + p + n
 g <- 28
 iv_r[g] <- iv_r[y] + p + n

 # iv(rmat) = iv(g) + p + 4*ps
 rmat <- 78
 iv_r[rmat] <- iv_r[g] + p + 4*ps

 # iv(qtr) = iv(rmat) + lh
 qtr <- 77
 iv_r[qtr] <- iv_r[rmat] + lh

 # iv(jtol) = iv(qtr) + p + n
 jtol <- 59
 iv_r[jtol] <- iv_r[qtr] + p + n

 jcn <- 66
 # iv(jcn) = iv(jtol) + 2*p
 iv_r[jcn] <- iv_r[jtol] + 2*p

 # iv(nextv) = iv(jcn) + p
 iv_r[nextv] <- iv_r[jcn] + p

 # IF (iv1 == 13) GO TO 999
 # This is the location of return to dglf/g from drglg with iv(1) = 14
 # [Back inside dglf/g, just after call to drglg with iv(1) = 13]
 # [Note that in this location iv1 = 12]
 iv_r[1] <- 14
 #
 # What is left is the final storage allocation inside dglf/g.
 #
 # !  ***  STORAGE ALLOCATION  ***
 # iv(d) = iv(nextv)
 d <- 27
 iv_r[d] <- iv_r[nextv]

 # iv(r) = iv(d) + p
 r <- 61
 iv_r[r] <- iv_r[d] + p

 # iv(regd0) = iv(r) + (p - ps + 1)*n
 regd0 <- 82
 if (hasAnalyticModelDeriv) {
   iv_r[regd0] <- iv_r[r] + (p - ps + 1)*n
 } else {
   iv_r[regd0] <- iv_r[r] + (p - ps + 3)*n
 }

 # iv(j) = iv(regd0) + ((p-ps+2)*(p-ps+1)/2)*n
 j <- 70
 iv_r[j] <- iv_r[regd0] + ((p-ps+2)*(p-ps+1)/2)*n

 # iv(nextv) = iv(j) + n*ps
 iv_r[nextv] <- iv_r[j] + n*ps

 # Note: Because iv1 = 12, we just keep going.
 # The next call is to drglg which goes ahead and starts the
 # estimation.
 # IF (iv1 == 13) GO TO 999
 #
 # 10   d1 = iv(d)
 # dr1 = iv(j)
 # r1 = iv(r)
 # rd1 = iv(regd0)
 #
 # However, here is the problem:  The first call to drglg
 # is with iv(1) = 14 at the very start of the estimation loop.
 # It is this call that goes ahead and checks parameters again.
 # And, we must make this part of the setup procedure.
 #
 # Call to drglg with iv(1) = 14
 # ps1 = ps + 1
 # iv1 = iv(1) = 14
 # IF (iv1 > 2) GO TO 10
 # 10   IF (nd < ps) GO TO 360
 # [Removed...]
 # IF (iv1 == 14) GO TO 30
 # 30   jtol1 = iv(jtol)
 # IF (v(dinit) >= zero) CALL dv7scp(p, d, v(dinit))
 # IF (v(dtinit) > zero) CALL dv7scp(p, v(jtol1), v(dtinit))
 # i = jtol1 + p
 # IF (v(d0init) > zero) CALL dv7scp(p, v(i), v(d0init))
 # NOTE: This gets initialized!
 # iv(nf0) = 0
 # iv(nf1) = 0
 #
 # 40   CONTINUE
 # g1 = iv(g)
 # y1 = iv(y)
 # CALL dg7lit_m(d, v(g1), iv, liv, lv, p, ps, v, x, v(y1))
 #
 # We now have a pseudo-call to dg7lit with iv(1) = 14
 # [Inside dg7lit...]
 # CALL dparck(alg, d, iv, liv, lv, p, v) [with iv(1) = 14]
 # iv1 = iv(1)
 # IF (iv1 /= 13 .AND. iv1 /= 12) GO TO 30
 # 30  CONTINUE
 # IF (n >= 1) GO TO 50
 # 50      IF (iv1 /= 14) iv(nextiv) = iv(perm)
 # [Removed]
 # k = iv(parsav) - epslon
 # ------------------------------------------
 # [From earlier]
 #    nvsave <- 9
 #    iv_r[vsave]  <- 58
 #    iv_r[parsav] <- iv_r[vsave] + nvsave
 #    => iv_r[parsav] = 58 + 9 = 67
 #    Also:  epslon <- 19
 # => k = 67 - 19 = 48 and k+1 = 49
 # kplus1 = k+1
 # CALL dv7dfl(alg1, lv-k, v(kplus1))
 # iv(dtype0) = 2 - alg1 = 1
 # iv(oldn) = n
 # The following loads strings into which that say:
 # 'NONDEFAULT V'
 # which(1) = dflt(1)
 # which(2) = dflt(2)
 # which(3) = dflt(3)
 # GO TO 110
 # 110  IF (iv1 == 14) iv1 = 12
 # IF (big > tiny) GO TO 120
 # tiny = dr7mdc(1)
 # machep = dr7mdc(3)
 # big = dr7mdc(6)
 # vm(12) = machep
 # vx(12) = big
 # vx(13) = big
 # vm(14) = machep
 # vm(17) = tiny
 # vx(17) = big
 # vm(18) = tiny
 # vx(18) = big
 # vx(20) = big
 # vx(21) = big
 # vx(22) = big
 # vm(24) = machep
 # vm(25) = machep
 # vm(26) = machep
 # vx(28) = dr7mdc(5)
 # vm(29) = machep
 # vx(30) = big
 # vm(33) = machep
 #
 # The code that checks v values starts at 120
 # 120  m = 0
 #
 # NOTES:
 # The code that exists at this stage does
 # a couple of things:
 # 1.  It checks v default values to make sure they are in range.
 #     If there are errors, it writes out a message for each error.
 # 2.  It checks initialization of D (multiple types).
 #     If there are errors, it writes messages.
 # 3.  It checks for an error with iv(inits) or iv(inith) [which are in the same place?].
 #     If so it writes a message [? I think.]
 # 4.  If no errors are detected, it then writes out messages summarizing
 #     which default values have been changed.
 #
 # Note that, because we have control at this stage, we do not need to
 # replicate every last bit of this code at this stage.  We can decide
 # what options we will allow the user to change using control, and then
 # implement checks for those.
 #
 setup <- list()
 setup$iv_r <- iv_r
 setup$v_r  <- v_r
return(setup)
}
