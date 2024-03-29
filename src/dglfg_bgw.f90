! Version for BGW R package. August 5, 2022
! The major issue is that R packages carefully check Fortran code
! for any WRITE or PRINT statements (regardless of whether they call them).
! So the entire package must be cleansed of any such statements.
! The main issue is output routines:  ditsum and dn3rdp.
! For now, we create dummy versions of these that do not write anything,
! so that their original place in the code is preserved.
! See DSB NOTES below for locations
! ------------------------------------------------------------------------
! Subroutine drglg removed from original dglfg.f
! Subroutines dglg and dglf removed
! Code converted using TO_F90 by Alan Miller
! Date: 2021-07-17  Time: 10:48:25

SUBROUTINE df7hes(d, g, irt, iv, liv, lv, p, v, x)

!  ***  COMPUTE FINITE-DIFFERENCE HESSIAN, STORE IT IN V STARTING
!  ***  AT V(IV(FDH)) = V(-IV(H)).

!  ***  IF IV(COVREQ) .GE. 0 THEN DF7HES USES GRADIENT DIFFERENCES,
!  ***  OTHERWISE FUNCTION DIFFERENCES.  STORAGE IN V IS AS IN DG7LIT.

! IRT VALUES...
!     1 = COMPUTE FUNCTION VALUE, I.E., V(F).
!     2 = COMPUTE G.
!     3 = DONE.


!  ***  PARAMETER DECLARATIONS  ***

INTEGER, INTENT(IN)                      :: p
INTEGER, INTENT(IN)                      :: liv
INTEGER, INTENT(IN)                      :: lv
DOUBLE PRECISION, INTENT(IN OUT)         :: d(p)
DOUBLE PRECISION, INTENT(IN OUT)         :: g(p)
INTEGER, INTENT(IN OUT)                  :: irt
INTEGER, INTENT(IN OUT)                  :: iv(liv)
DOUBLE PRECISION, INTENT(IN OUT)         :: v(lv)
DOUBLE PRECISION, INTENT(IN OUT)         :: x(p)

!  ***  LOCAL VARIABLES  ***

INTEGER :: gsave1, flo1, hes, hmi, hpi, hpm, i, k, kind, l, m, mm1, mm1o2,  &
    pp1o2, stpi, stpm, stp0
DOUBLE PRECISION :: del

!  ***  EXTERNAL SUBROUTINES  ***

EXTERNAL dv7cpy

! DV7CPY.... COPY ONE VECTOR TO ANOTHER.

!  ***  SUBSCRIPTS FOR IV AND V  ***

DOUBLE PRECISION, PARAMETER :: half=0.5D+0
DOUBLE PRECISION, PARAMETER :: negpt5=-0.5D+0
DOUBLE PRECISION, PARAMETER :: one=1.d+0
DOUBLE PRECISION, PARAMETER :: two=2.d+0
DOUBLE PRECISION, PARAMETER :: zero=0.d+0

INTEGER, PARAMETER :: covreq=15
INTEGER, PARAMETER :: delta=52
INTEGER, PARAMETER :: delta0=44
INTEGER, PARAMETER :: dltfdc=42
INTEGER, PARAMETER :: f=10
INTEGER, PARAMETER :: fdh=74
INTEGER, PARAMETER :: fx=53
INTEGER, PARAMETER :: h=56
INTEGER, PARAMETER :: kagqt=33
INTEGER, PARAMETER :: mode=35
INTEGER, PARAMETER :: nfgcal=7
INTEGER, PARAMETER :: savei=63
INTEGER, PARAMETER :: switch=12
INTEGER, PARAMETER :: toobig=2
INTEGER, PARAMETER :: w=65
INTEGER, PARAMETER :: xmsave=51

! DSB NOTE:  The following was added to stop compiler warnings
hes = 0
flo1 = 1
!+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++

irt = 4
kind = iv(covreq)
m = iv(mode)
IF (m > 0) GO TO 10
iv(h) = -IABS(iv(h))
iv(fdh) = 0
iv(kagqt) = -1
v(fx) = v(f)
10   IF (m > p) GO TO 999
IF (kind < 0) GO TO 110

!  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING BOTH FUNCTION AND
!  ***  GRADIENT VALUES.

gsave1 = iv(w) + p
IF (m > 0) GO TO 20
!        ***  FIRST CALL ON DF7HES.  SET GSAVE = G, TAKE FIRST STEP  ***
CALL dv7cpy(p, v(gsave1), g)
iv(switch) = iv(nfgcal)
GO TO 90

20   del = v(delta)
x(m) = v(xmsave)
IF (iv(toobig) == 0) GO TO 40

!     ***  HANDLE OVERSIZE V(DELTA)  ***

IF (del*x(m) > zero) GO TO 30
!             ***  WE ALREADY TRIED SHRINKING V(DELTA), SO QUIT  ***
iv(fdh) = -2
GO TO 220

!        ***  TRY SHRINKING V(DELTA)  ***
30      del = negpt5 * del
GO TO 100

40   hes = -iv(h)

!  ***  SET  G = (G - GSAVE)/DEL  ***

DO  i = 1, p
  g(i) = (g(i) - v(gsave1)) / del
  gsave1 = gsave1 + 1
END DO

!  ***  ADD G AS NEW COL. TO FINITE-DIFF. HESSIAN MATRIX  ***

k = hes + m*(m-1)/2
l = k + m - 2
IF (m == 1) GO TO 70

!  ***  SET  H(I,M) = 0.5 * (H(I,M) + G(I))  FOR I = 1 TO M-1  ***

mm1 = m - 1
DO  i = 1, mm1
  v(k) = half * (v(k) + g(i))
  k = k + 1
END DO

!  ***  ADD  H(I,M) = G(I)  FOR I = M TO P  ***

70   l = l + 1
DO  i = m, p
  v(l) = g(i)
  l = l + i
END DO

90   m = m + 1
iv(mode) = m
IF (m > p) GO TO 210

!  ***  CHOOSE NEXT FINITE-DIFFERENCE STEP, RETURN TO GET G THERE  ***

del = v(delta0) *   MAX(one/d(m),  ABS(x(m)))
IF (x(m) < zero) del = -del
v(xmsave) = x(m)
100  x(m) = x(m) + del
v(delta) = del
irt = 2
GO TO 999

!  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING FUNCTION VALUES ONLY.

110  stp0 = iv(w) + p - 1
mm1 = m - 1
mm1o2 = m*mm1/2
IF (m > 0) GO TO 120
!        ***  FIRST CALL ON DF7HES.  ***
iv(savei) = 0
GO TO 200

120  i = iv(savei)
hes = -iv(h)
IF (i > 0) GO TO 180
IF (iv(toobig) == 0) GO TO 140

!     ***  HANDLE OVERSIZE STEP  ***

stpm = stp0 + m
del = v(stpm)
IF (del*x(xmsave) > zero) GO TO 130
!             ***  WE ALREADY TRIED SHRINKING THE STEP, SO QUIT  ***
iv(fdh) = -2
GO TO 220

!        ***  TRY SHRINKING THE STEP  ***
130     del = negpt5 * del
x(m) = x(xmsave) + del
v(stpm) = del
irt = 1
GO TO 999

!  ***  SAVE F(X + STP(M)*E(M)) IN H(P,M)  ***

140  pp1o2 = p * (p-1) / 2
hpm = hes + pp1o2 + mm1
v(hpm) = v(f)

!  ***  START COMPUTING ROW M OF THE FINITE-DIFFERENCE HESSIAN H.  ***

hmi = hes + mm1o2
IF (mm1 == 0) GO TO 160
hpi = hes + pp1o2
DO  i = 1, mm1
  v(hmi) = v(fx) - (v(f) + v(hpi))
  hmi = hmi + 1
  hpi = hpi + 1
END DO
160  v(hmi) = v(f) - two*v(fx)

!  ***  COMPUTE FUNCTION VALUES NEEDED TO COMPLETE ROW M OF H.  ***

i = 1

170  iv(savei) = i
stpi = stp0 + i
v(delta) = x(i)
x(i) = x(i) + v(stpi)
IF (i == m) x(i) = v(xmsave) - v(stpi)
irt = 1
GO TO 999

180  x(i) = v(delta)
IF (iv(toobig) == 0) GO TO 190
!        ***  PUNT IN THE EVENT OF AN OVERSIZE STEP  ***
iv(fdh) = -2
GO TO 220

!  ***  FINISH COMPUTING H(M,I)  ***

190  stpi = stp0 + i
hmi = hes + mm1o2 + i - 1
stpm = stp0 + m
v(hmi) = (v(hmi) + v(f)) / (v(stpi)*v(stpm))
i = i + 1
IF (i <= m) GO TO 170
iv(savei) = 0
x(m) = v(xmsave)

200  m = m + 1
iv(mode) = m
IF (m > p) GO TO 210

!  ***  PREPARE TO COMPUTE ROW M OF THE FINITE-DIFFERENCE HESSIAN H.
!  ***  COMPUTE M-TH STEP SIZE STP(M), THEN RETURN TO OBTAIN
!  ***  F(X + STP(M)*E(M)), WHERE E(M) = M-TH STD. UNIT VECTOR.

del = v(dltfdc) *   MAX(one/d(m),  ABS(x(m)))
IF (x(m) < zero) del = -del
v(xmsave) = x(m)
x(m) = x(m) + del
stpm = stp0 + m
v(stpm) = del
irt = 1
GO TO 999

!  ***  RESTORE V(F), ETC.  ***

210  iv(fdh) = hes
220  v(f) = v(fx)
irt = 3
IF (kind < 0) GO TO 999
iv(nfgcal) = iv(switch)
gsave1 = iv(w) + p
CALL dv7cpy(p, g, v(gsave1))
GO TO 999

999  RETURN
!  ***  LAST LINE OF DF7HES FOLLOWS  ***
END SUBROUTINE df7hes

SUBROUTINE dg2lrd(dr, iv, l, lh, liv, lv, nd, n, p, ps, r, rd,  &
    rhoi, rhor, v, w, x, z)

!  ***  COMPUTE REGRESSION DIAGNOSTIC FOR  DRGLG  ***

!  ***  PARAMETERS  ***

INTEGER, INTENT(IN)                      :: lh
INTEGER, INTENT(IN)                      :: liv
INTEGER, INTENT(IN)                      :: lv
INTEGER, INTENT(IN)                      :: nd
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: p
INTEGER, INTENT(IN)                      :: ps
DOUBLE PRECISION, INTENT(IN OUT)         :: dr(nd,p)
INTEGER, INTENT(IN OUT)                  :: iv(liv)
DOUBLE PRECISION, INTENT(IN OUT)         :: l(lh)

DOUBLE PRECISION, INTENT(IN OUT)         :: r(n)
DOUBLE PRECISION, INTENT(IN OUT)         :: rd(n)
INTEGER, INTENT(IN OUT)                  :: rhoi(*)
DOUBLE PRECISION, INTENT(IN OUT)         :: rhor(*)
DOUBLE PRECISION, INTENT(IN OUT)         :: v(lv)
DOUBLE PRECISION, INTENT(IN OUT)         :: w(p)
DOUBLE PRECISION, INTENT(IN OUT)         :: x(p)
DOUBLE PRECISION, INTENT(IN OUT)         :: z(p)

!  ***  CODED BY DAVID M. GAY (SPRING 1986, SUMMER 1991)  ***

!  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***

EXTERNAL dd7tpr, dl7itv, dl7ivm,dl7srt, dl7sqr, ds7lvm, dv2axy,dv7cpy, dv7scp
DOUBLE PRECISION :: dd7tpr

! DD7TPR... COMPUTES INNER PRODUCT OF TWO VECTORS.
! DL7ITV... MULTIPLIES INVERSE TRANSPOSE OF LOWER TRIANGLE TIMES VECTOR.
! DL7IVM... APPLY INVERSE OF COMPACT LOWER TRIANG. MATRIX.
! DL7SRT.... COMPUTES CHOLESKY FACTOR OF (LOWER TRIANG. OF) SYM. MATRIX.
! DL7SQR... COMPUTES L*(L**T) FOR LOWER TRIANG. MATRIX L.
! DS7LVM... MULTIPLIES COMPACTLY STORED SYM. MATRIX TIMES VECTOR.
! DV2AXY.... ADDS A MULTIPLE OF ONE VECTOR TO ANOTHER.
! DV7CPY.... COPIES ONE VECTOR TO ANOTHER.
! DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR.

!  ***  LOCAL VARIABLES  ***

LOGICAL :: useflo
INTEGER :: bs1, bsinc, flo1, floinc, h1, hps1, i,  &
    j, j1, k, ki, ki1, kid, l1, LE, ll, loo1, n1,  &
    pmps, pp1o2, ps1, px, rdr, xni, zap1, zaplen
DOUBLE PRECISION :: frac, hi, ri, s, t, t1

!  ***  CONSTANTS  ***

!  ***  IV SUBSCRIPTS  ***

INTEGER, PARAMETER :: bs=85
INTEGER, PARAMETER :: bsstr=86
INTEGER, PARAMETER :: covreq=15
INTEGER, PARAMETER :: fdh=74
INTEGER, PARAMETER :: flo=88
INTEGER, PARAMETER :: flostr=89
INTEGER, PARAMETER :: loo=84
INTEGER, PARAMETER :: nb=87
INTEGER, PARAMETER :: nfix=83
INTEGER, PARAMETER :: rdreq=57
INTEGER, PARAMETER :: regd=67
INTEGER, PARAMETER :: xnoti=90
DOUBLE PRECISION, PARAMETER :: half=0.5D+0
DOUBLE PRECISION, PARAMETER :: negone=-1.d+0
DOUBLE PRECISION, PARAMETER :: one=1.d+0
DOUBLE PRECISION, PARAMETER :: zero=0.d+0

! DSB NOTE:  The following was added to stop compiler warnings
flo1 = 1
floinc = 1
!++++++++++++++++++++++++++++++++  BODY  +++++++++++++++++++++++++++++++

i = iv(rdreq)
rdr = MOD(i/2, 3)
IF (rdr == 0) GO TO 999
h1 = iv(fdh)
useflo = .false.
px = p
n1 = n
frac = one
xni = 0
IF (rdr == 1) GO TO 120
loo1 = iv(loo)
IF (loo1 <= 0 .OR. loo1 > 6) THEN
  iv(regd) = -1
  GO TO 999
END IF
IF (loo1 > 3) THEN
  useflo = .true.
  flo1 = iv(flo)
  floinc = iv(flostr)
  loo1 = loo1 - 3
END IF
xni = iv(xnoti)
px = p - iv(nfix)
IF (px < ps .OR. px > p) THEN
  iv(regd) = -2
  GO TO 999
END IF
IF (loo1 == 1) GO TO 120
n1 = iv(nb)
IF (n1 <= 0 .OR. n1 > n) THEN
  iv(regd) = -3
  GO TO 999
END IF
bs1 = iv(bs)
bsinc = iv(bsstr)
IF (h1 <= 0) GO TO 190
IF (IABS(iv(covreq)) >= 3) CALL dl7sqr(p, v(h1), l)
pp1o2 = px*(px+1)/2
ps1 = ps + 1
zap1 = ps*(ps1)/2 + 1
LE = 0
DO  i = 1, n1
  IF (useflo) THEN
    frac = rhor(flo1)
    flo1 = flo1 + floinc
  END IF
  l1 = LE + 1
  IF (l1 > n) EXIT
  LE = LE + rhoi(bs1)
  IF (LE > n) LE = n
  bs1 = bs1 + bsinc
  CALL dv7cpy(pp1o2, l, v(h1))
  IF (ps >= px) GO TO 50
  k = zap1
  ki = l1
  DO  j = ps1, p
    ki = ki + n
    ki1 = ki
    DO  ll = l1, LE
      CALL dv2axy(ps, l(k), -frac*rd(ki1), dr(1,ll), l(k))
      ki1 = ki1 + 1
    END DO
    k = k + ps
    DO  j1 = ps1, j
      ki = ki + n
      ki1 = ki
      t = zero
      DO  ll = l1, LE
        t = t + rd(ki1)
        ki1 = ki1 + 1
      END DO
      l(k) = l(k) - frac*t
      k = k + 1
    END DO
  END DO
  50      DO  ll = l1, LE
    t = -frac*rd(ll)
    k = 1
    DO  j = 1, ps
      CALL dv2axy(j, l(k), t*dr(j,ll), dr(1,ll), l(k))
      k = k + j
    END DO
  END DO
  CALL dl7srt(1, px, l, l, j)
  IF (j == 0) THEN
    CALL dv7scp(px, w, zero)
    DO  ll = l1, LE
      CALL dv2axy(ps, w, r(ll), dr(1,ll), w)
      IF (ps1 > px) CYCLE
      k = l1
      DO  j = ps1, p
        k = k + n
        w(j) = w(j) + r(k)
      END DO
    END DO
    CALL dl7ivm(px, w, l, w)
    CALL dl7itv(px, w, l, w)
    CALL ds7lvm(px, z, v(h1), w)
    rd(i) = half * frac * dd7tpr(px, w, z)
    IF (xni > 0) THEN
      CALL dv2axy(px, rhor(xni), frac, w, x)
      xni = xni + px
    END IF
  ELSE
    rd(i) = negone
    IF (xni > 0) THEN
      CALL dv7cpy(px, rhor(xni), x)
      xni = xni + px
    END IF
  END IF
END DO
! 110  iv(regd) = 1
iv(regd) = 1
!     *** RESTORE L ***
CALL dl7srt(1, p, l, v(h1), j)
GO TO 999

120  IF (h1 <= 0) GO TO 190
IF (IABS(iv(covreq)) >= 3) CALL dl7sqr(p, v(h1), l)
IF (ps >= px) GO TO 170
ps1 = ps + 1
pmps = px - ps
zap1 = ps*(ps1)/2
zaplen = px*(px+1)/2 - zap1
hps1 = h1 + zap1
zap1 = zap1 + 1
DO  i = 1, n
  IF (useflo) THEN
    frac = rhor(flo1)
    flo1 = flo1 + floinc
  END IF
  CALL dv7cpy(zaplen, l(zap1), v(hps1))
  CALL dv7scp(ps, w, zero)
  k = zap1
  ki = i
  kid = ki
  DO  j = ps1, px
    ki = ki + n
    CALL dv2axy(ps, l(k), -frac*rd(ki), dr(1,i), l(k))
    k = k + ps
    kid = kid + n
    w(j) = frac*r(kid)
    DO  j1 = ps1, j
      ki = ki + n
      l(k) =  l(k) - frac*rd(ki)
      k = k + 1
    END DO
  END DO
  CALL dl7srt(ps1, px, l, l, j)
  IF (j /= 0) GO TO 150
  CALL dv7cpy(ps, z, dr(1,i))
  CALL dv7scp(pmps, z(ps1), zero)
  CALL dl7ivm(px, z, l, z)
  hi = dd7tpr(px, z, z)
  CALL dl7ivm(px, w, l, w)
  ri = frac*r(i)
!        *** FIRST PS ELEMENTS OF W VANISH ***
  t = dd7tpr(pmps, w(ps1), z(ps1))
  s = frac*rd(i)
  t1 = one - s*hi
  IF (t1 <= zero) GO TO 150
  CALL dv2axy(px, w, (ri + s*t)/t1, z, w)
  CALL dl7itv(px, w, l, w)
  CALL ds7lvm(px, z, v(h1), w)
  rd(i) = half * dd7tpr(px, w, z)
  IF (xni > 0) THEN
    CALL dv2axy(px, rhor(xni), one, w, x)
    xni = xni + px
  END IF
  CYCLE
  150     rd(i) = negone
  IF (xni > 0) THEN
    CALL dv7cpy(px, rhor(xni), x)
    xni = xni + px
  END IF
END DO

!     *** RESTORE L ***

CALL dv7cpy(zaplen, l(zap1), v(hps1))
CALL dl7srt(ps1, px, l, l, j)
GO TO 200

170  DO  i = 1, n
  IF (useflo) THEN
    frac = rhor(flo1)
    flo1 = flo1 + floinc
  END IF
  CALL dl7ivm(px, z, l, dr(1,i))
  s = dd7tpr(px, z, z)
  t = one - frac*rd(i) * s
  IF (t <= zero) THEN
    rd(i) = negone
    IF (xni > 0) THEN
      CALL dv7cpy(px, rhor(xni), x)
      xni = xni + px
    END IF
  ELSE
    rd(i) = half * frac * (r(i)/t)**2 * s
    IF (xni > 0) THEN
      CALL dl7itv(px, z, l, z)
      CALL dv2axy(px, rhor(xni), frac*r(i)/t, z, x)
      xni = xni + px
    END IF
  END IF
END DO
GO TO 200

190  CALL dv7scp(n1, rd, negone)
200  iv(regd) = 1

999  RETURN
!  ***  LAST LINE OF DG2LRD FOLLOWS  ***
END SUBROUTINE dg2lrd

! DSB NOTE:  dg7lit has been removed. It is been moved to
! the drglg module.
! SUBROUTINE dg7lit(d, g, iv, liv, lv, p, ps, v, x, y)

SUBROUTINE dl7nvr(n, lin, l)

!  ***  COMPUTE  LIN = L**-1,  BOTH  N X N  LOWER TRIANG. STORED   ***
!  ***  COMPACTLY BY ROWS.  LIN AND L MAY SHARE THE SAME STORAGE.  ***

!  ***  PARAMETERS  ***


INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(OUT)            :: lin(*)
DOUBLE PRECISION, INTENT(IN)             :: l(*)


!     DIMENSION L(N*(N+1)/2), LIN(N*(N+1)/2)

!  ***  LOCAL VARIABLES  ***

INTEGER :: i, ii, im1, jj, j0, j1, k, k0, np1
DOUBLE PRECISION ::  t
DOUBLE PRECISION, PARAMETER :: one=1.d+0
DOUBLE PRECISION, PARAMETER :: zero=0.d+0

!  ***  BODY  ***

np1 = n + 1
j0 = n*(np1)/2
DO  ii = 1, n
  i = np1 - ii
  lin(j0) = one/l(j0)
  IF (i <= 1) EXIT
  j1 = j0
  im1 = i - 1
  DO  jj = 1, im1
    t = zero
    j0 = j1
    k0 = j1 - jj
    DO  k = 1, jj
      t = t - l(k0)*lin(j0)
      j0 = j0 - 1
      k0 = k0 + k - i
    END DO
    lin(j0) = t/l(k0)
  END DO
  j0 = j0 - 1
END DO
! 999  RETURN
RETURN
!  ***  LAST LINE OF DL7NVR FOLLOWS  ***
END SUBROUTINE dl7nvr

SUBROUTINE dl7tsq(n, a, l)

!  ***  SET A TO LOWER TRIANGLE OF (L**T) * L  ***

!  ***  L = N X N LOWER TRIANG. MATRIX STORED ROWWISE.  ***
!  ***  A IS ALSO STORED ROWWISE AND MAY SHARE STORAGE WITH L.  ***


INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(OUT)            :: a(*)
DOUBLE PRECISION, INTENT(IN)             :: l(*)


!     DIMENSION A(N*(N+1)/2), L(N*(N+1)/2)

INTEGER :: i, ii, iim1, i1, j, k, m
DOUBLE PRECISION :: lii, lj

ii = 0
DO  i = 1, n
  i1 = ii + 1
  ii = ii + i
  m = 1
  IF (i == 1) GO TO 30
  iim1 = ii - 1
  DO  j = i1, iim1
    lj = l(j)
    DO  k = i1, j
      a(m) = a(m) + lj*l(k)
      m = m + 1
    END DO
  END DO
  30      lii = l(ii)
  DO  j = i1, ii
    a(j) = lii * l(j)
  END DO
END DO

! 999  RETURN
RETURN
!  ***  LAST LINE OF DL7TSQ FOLLOWS  ***
END SUBROUTINE dl7tsq

! DSB NOTE:  See discussion at the beginning
! This is the routine that writes regression diagnostics.
! We do not need it at the moment, but we replace it
! with a dummy as a placeholder so we can find it later.
! SUBROUTINE dn3rdp(iv, liv, lv, n, p, rd, rhoi, rhor, v)
SUBROUTINE dn3rdp(inDummy, outDummy)
!  ***  PRINT REGRESSION DIAGNOSTICS FOR MLPSL AND NL2S1 ***

INTEGER :: inDummy, outDummy

outDummy = inDummy
RETURN
!  ***  LAST LINE OF DN3RDP FOLLOWS  ***
END SUBROUTINE dn3rdp
