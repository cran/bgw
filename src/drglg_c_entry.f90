!module drglg_mod
!    use, intrinsic :: iso_c_binding
!    use ieee_arithmetic
!
!    implicit none
!    private
!    ! public :: drglg, divset_f, dparck_m, dg7lit_m, rho
!    public :: drglg, divset_f
!
!    ! Subroutines in this file
!    ! drglg, divset_f
!    ! dg7lit_m, dparck_m, rho
!
!contains
    SUBROUTINE drglg(d, dr, iv, liv, lv, n, nd, nn, p, ps, &
       lrhoi, lrhor, r, rd, v, x, rhoi, rhor, i_itsum) bind(C, name="drglg_")

    ! SUBROUTINE drglg(d, dr, iv, liv, lv, n, nd, nn, p, ps, &
    !   lrhoi, lrhor, r, rd, v, x, rhoi, rhor) bind(C, name="drglg_")

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2021-07-14  Time: 16:45:17

    !      SUBROUTINE  DRGLG(D, DR, IV, LIV, LV, N, ND, NN, P, PS, R,
    !     1                  RD, RHO, RHOI, RHOR, V, X)

    ! *** ITERATION DRIVER FOR GENERALIZED (NON)LINEAR MODELS (ETC.)
    use, intrinsic :: iso_c_binding
    implicit none

    ! integer(kind = c_int), INTENT(IN OUT)    :: liv
    ! integer(kind = c_int), INTENT(IN OUT)    :: lv
    integer(kind = c_int), INTENT(IN)          :: liv
    integer(kind = c_int), INTENT(IN)          :: lv
    integer(kind = c_int), INTENT(IN)          :: n
    integer(kind = c_int), INTENT(IN)          :: nd
    integer(kind = c_int), INTENT(IN)          :: nn
    integer(kind = c_int), INTENT(IN)          :: p
    integer(kind = c_int), INTENT(IN)          :: ps
    ! integer(kind = c_int), INTENT(IN OUT)    :: lrhoi
    ! integer(kind = c_int), INTENT(IN OUT)    :: lrhor
    integer(kind = c_int), INTENT(IN), value   :: lrhoi
    integer(kind = c_int), INTENT(IN), value   :: lrhor
    real(kind = c_double), INTENT(IN OUT)      :: d(p)
    real(kind = c_double), INTENT(IN OUT)      :: dr(nd,n)
    integer(kind = c_int), INTENT(IN OUT)      :: iv(liv)
    real(kind = c_double), INTENT(IN OUT)      :: r(*)
    real(kind = c_double), INTENT(IN OUT)      :: rd(*)
    real(kind = c_double), INTENT(IN OUT)      :: v(lv)
    real(kind = c_double), INTENT(IN OUT)      :: x(p)
    ! integer(kind = c_int), INTENT(IN OUT)      :: rhoi(*)
    ! real(kind = c_double), INTENT(INOUT)       :: rhor(*)
    integer(kind = c_int), INTENT(IN OUT), dimension(lrhoi)  :: rhoi
    real(kind = c_double), INTENT(INOUT), dimension(lrhor)   :: rhor

    integer(kind = c_int), INTENT(IN OUT)      :: i_itsum

    !      INTEGER IV(LIV), RHOI(*)
    !      DOUBLE PRECISION D(P), DR(ND,N), R(*), RD(*), RHOR(*),
    !      DIMENSION RD(N, (P-PS)*(P-PS+1)/2 + 1)
    ! EXTERNAL rho

    !--------------------------  PARAMETER USAGE  --------------------------

    ! D....... SCALE VECTOR.
    ! DR...... DERIVATIVES OF R AT X.
    ! IV...... INTEGER VALUES ARRAY.
    ! LIV..... LENGTH OF IV... LIV MUST BE AT LEAST P + 90.
    ! LV...... LENGTH OF V...  LV  MUST BE AT LEAST
    !              105 + P*(2*P+16) + 2*N + 4*PS.
    ! N....... TOTAL NUMBER OF RESIDUALS.
    ! ND...... LEADING DIMENSION OF DR -- MUST BE AT LEAST PS.
    ! NN...... LEAD DIMENSION OF R, RD.
    ! P....... NUMBER OF PARAMETERS (COMPONENTS OF X) BEING ESTIMATED.
    ! PS...... NUMBER OF NON-NUISANCE PARAMETERS.
    ! R....... RESIDUALS (OR MEANS -- FUNCTIONS OF X) WHEN  DRGLG IS CALLED
    !          WITH IV(1) = 1.
    ! RD...... RD(I) = HALF * (G(I)**T * H(I)**-1 * G(I)) ON OUTPUT WHEN
    !          IV(RDREQ) IS 2, 3, 5, OR 6.   DRGLG SETS IV(REGD) = 1 IF RD
    !          IS SUCCESSFULLY COMPUTED, TO 0 IF NO ATTEMPT WAS MADE
    !          TO COMPUTE IT, AND TO -1 IF H (THE FINITE-DIFFERENCE HESSIAN)
    !          WAS INDEFINITE.  BEFORE CONVERGENCE, RD IS ALSO USED AS
    !          TEMPORARY STORAGE.
    ! RHO..... COMPUTES INFO ABOUT OBJECTIVE FUNCTION.
    ! RHOI.... PASSED WITHOUT CHANGE TO RHO.
    ! RHOR.... PASSED WITHOUT CHANGE TO RHO.
    ! V....... FLOATING-POINT VALUES ARRAY.
    ! X....... PARAMETER VECTOR BEING ESTIMATED (INPUT = INITIAL GUESS,
    !              OUTPUT = BEST VALUE FOUND).

    ! *** CALLING SEQUENCE FOR RHO...

    !  CALL RHO(NEED, F, N, NF, XN, R, RD, RHOI, RHOR, W)

    !  PARAMETER DECLARATIONS FOR RHO...

    ! INTEGER NEED(2), N, NF, RHOI(*)
    ! FLOATING-POINT F, XN(*), R(*), RD(N,*), RHOR(*), W(N)

    !    RHOI AND RHOR ARE FOR RHO TO USE AS IT SEES FIT.  THEY ARE PASSED
    ! TO RHO WITHOUT CHANGE.  IF IV(RDREQ) IS AT LEAST 4, I.E., IF MORE
    ! THAN THE SIMPLEST REGRESSION DIAGNOSTIC INFORMATION IS TO BE COMPUTED,
    ! THEN SOME COMPONENTS OF RHOI AND RHOR MUST CONVEY SOME EXTRA
    ! DETAILS, AS DESCRIBED BELOW.
    !    F, R, RD, AND W ARE EXPLAINED BELOW WITH NEED.
    !    XN IS THE VECTOR OF NUISANCE PARAMETERS (OF LENGTH P - PS).  IF
    ! RHO NEEDS TO KNOW THE LENGTH OF XN, THEN THIS LENGTH SHOULD BE
    ! COMMUNICATED THROUGH RHOI (OR THROUGH COMMON).  RHO SHOULD NOT CHANGE
    ! XN.
    !    NEED(1) = 1 MEANS RHO SHOULD SET F TO THE SUM OF THE LOSS FUNCTION
    ! VALUES AT THE RESIDUALS R(I).  NF IS THE CURRENT FUNCTION INVOCATION
    ! COUNT (A VALUE THAT IS INCREMENTED EACH TIME A NEW PARAMETER EXTIMATE
    ! X IS CONSIDERED).  NEED(2) IS THE VALUE NF HAD AT THE LAST R WHERE
    ! RHO MIGHT BE CALLED WITH NEED(1) = 2.  IF RHO SAVES INTERMEDIATE
    ! RESULTS FOR USE IN CALLS WITH NEED(1) = 2, THEN IT CAN USE NF TO TELL
    ! WHICH INTERMEDIATE RESULTS ARE APPROPRIATE, AND IT CAN SAVE SOME OF
    ! THESE RESULTS IN R.
    !    NEED(1) = 2 MEANS RHO SHOULD SET R(I) TO THE LOSS FUNCTION
    ! DERIVATIVE WITH RESPECT TO THE RESIDUALS THAT WERE PASSED TO RHO WHEN
    ! NF HAD THE SAME VALUE IT DOES NOW (AND NEED(1) WAS 1).  RHO SHOULD
    ! ALSO SET W(I) TO THE APPROXIMATION OF THE SECOND DERIVATIVE OF THE
    ! LOSS FUNCTION (WITH RESPECT TO THE I-TH RESIDUAL) THAT SHOULD BE USED
    ! IN THE GAUSS-NEWTON MODEL.  WHEN THERE ARE NUISANCE PARAMETERS (I.E.,
    ! WHEN PS .LT. P) RHO SHOULD ALSO SET R(I+K*N) TO THE DERIVATIVE OF THE
    ! LOSS FUNCTION WITH RESPECT TO THE I-TH RESIDUAL AND XN(K), AND IT
    ! SHOULD SET RD(I,J + K*(K+1)/2 + 1) TO THE SECOND PARTIAL DERIVATIVE
    ! OF THE I-TH RESIDUAL WITH RESPECT TO XN(J) AND XN(K), 0 .LE. J .LE. K
    ! AND 1 .LE. K .LE. P - PS, WHERE XN(0) MEANS THE I-TH RESIDUAL ITSELF.
    ! IN ANY EVENT, RHO SHOULD ALSO SET RD(I,1) TO THE (TRUE) SECOND
    ! DERIVATIVE OF THE LOSS FUNCTION WITH RESPECT TO THE I-TH RESIDUAL.
    !    NF (THE FUNCTION INVOCATION COUNT WHOSE NORMAL USE IS EXPLAINED
    ! ABOVE) SHOULD NOT BE CHANGED UNLESS RHO CANNOT CARRY OUT THE REQUESTED
    ! TASK, IN WHICH CASE RHO SHOULD SET NF TO 0.


    !  ***  REGRESSION DIAGNOSTICS  ***

    ! IV(RDREQ) INDICATES WHETHER A COVARIANCE MATRIX AND REGRESSION
    ! DIAGNOSTIC VECTOR ARE TO BE COMPUTED.  IV(RDREQ) HAS THE FORM
    ! IV(RDREQ) = CVR +2*RDR, WHERE CVR = 0 OR 1 AND RDR = 0, 1, OR 2,
    ! SO THAT

    !      CVR = MOD(IV(RDREQ), 2)
    !      RDR = MOD(IV(RDREQ)/2, 3).

    !    CVR = 0 FOR NO COVARIANCE MATRIX
    !        = 1 IF A COVARIANCE MATRIX ESTIMATE IS DESIRED

    !    RDR = 0 FOR NO LEAVE-ONE-OUT DIAGNOSTIC INFORMATION.
    !        = 1 TO HAVE ONE-STEP ESTIMATES OF F(X(I)) - F(X*) STORED IN RD,
    !            WHERE X(I) MINIMIZES F (THE OBJECTIVE FUNCTION) WITH
    !            COMPONENT I OF R REMOVED AND X* MINIMIZES THE FULL F.
    !        = 2 FOR MORE DETAILED ONE-STEP LEAVE-ONE-OUT INFORMATION, AS
    !            DICTATED BY THE IV COMPONENTS DESCRIBED BELOW.

    ! FOR RDR = 2, THE FOLLOWING COMPONENTS OF IV ARE RELEVANT...

    !  NFIX = IV(83) = NUMBER OF TRAILING NUISANCE PARAMETERS TO TREAT
    !          AS FIXED WHEN COMPUTING DIAGNOSTIC VECTORS (0 .LE. NFIX .LE.
    !          P - PS, SO X(I) IS KEPT FIXED FOR P - NFIX .LT. I .LE. P).

    !   LOO = IV(84) TELLS WHAT TO LEAVE OUT...
    !       = 1 MEANS LEAVE OUT EACH COMPONENT OF R SEPARATELY, AND
    !       = 2 MEANS LEAVE OUT CONTIGUOUS BLOCKS OF R COMPONENTS.
    !           FOR LOO = 2, IV(85) IS THE STARTING SUBSCRIPT IN RHOI
    !           OF AN ARRAY BS OF BLOCK SIZES, IV(86) IS THE STRIDE FOR BS,
    !           AND IV(87) = NB IS THE NUMBER OF BLOCKS, SO THAT
    !           BS(I) = RHOI(IV(85) + (I-1)*IV(86)), 1 .LE. I .LE. NB.
    !           NOTE THAT IF ALL BLOCKS ARE THE SAME SIZE, THEN IT SUFFICES
    !           TO SET RHOI(IV(85)) = BLOCKSIZE AND IV(86) = 0.
    !           NOTE THAT LOO = 1 IS EQUIVALENT TO LOO = 2 WITH
    !           RHOI(IV(85)) = 1, IV(86) = 0, IV(87) = N.
    !       = 3,4 ARE SIMILAR TO LOO = 1,2, RESPECTIVELY, BUT LEAVING A
    !           FRACTION OUT.  IN THIS CASE, IV(88) IS THE STARTING
    !           SUBSCRIPT IN RHOR OF AN ARRAY FLO OF FRACTIONS TO LEAVE OUT,
    !           AND IV(89) IS THE STRIDE FOR FLO...
    !           FLO(I) = RHOR(IV(88) + (I-1)*IV(89)), 1 .LE. I .LE. NB.

    ! XNOTI = IV(90) TELLS WHAT DIAGNOSTIC INFORMATION TO STORE...
    !       = 0  MEANS JUST STORE ONE-STEP ESTIMATES OF F(X(I)) - F(X*) IN
    !            RD(I), 1 .LE. I .LE. NB.
    !       .GT. 0 MEANS ALSO STORE ONE-STEP ESTIMATES OF X(I) ESTIMATES
    !            IN RHOR, STARTING AT RHOR(XNOTI)...
    !              X(I)(J) = RHOR((I-1)*(P-NFIX) + J + XNOTI-1),
    !              1 .LE. I .LE. NB, 1 .LE. J .LE. P - NFIX.

    !    SOMETIMES ONE-STEP ESTIMATES OF X(I) DO NOT EXIST, BECAUSE THE
    ! APPROXIMATE UPDATED HESSIAN MATRIX IS INDEFINITE.  IN SUCH CASES,
    ! THE CORRESPONDING RD COMPONENT IS SET TO -1, AND, IF XNOTI IS
    ! POSITIVE, THE SOLUTION X IS RETURNED AS X(I).  WHEN ONE-STEP ESTIMATES
    ! OF X(I) DO EXIST, THE CORRESPONDING COMPONENT OF RD IS POSITIVE.

    ! SUMMARY OF RHOI COMPONENTS (FOR RDR = MOD(IV(RDREQ)/2, 3) = 2)...

    ! IV(83) = NFIX
    ! IV(84) = LOO
    ! IV(85) = START IN RHOI OF BS
    ! IV(86) = STRIDE FOR BS
    ! IV(87) = NB
    ! IV(88) = START IN RHOR OF FLO
    ! IV(89) = STRIDE FOR FLO
    ! IV(90) = XNOTI (START IN RHOR OF X(I)).


    !  ***  COVARIANCE MATRIX ESTIMATE  ***

    ! IF IV(RDREQ) INDICATES THAT A COVARIANCE MATRIX IS TO BE COMPUTED,
    ! THEN IV(COVREQ) = IV(15) DETERMINES THE FORM OF THE COMPUTED
    ! COVARIANCE MATRIX ESTIMATE AND, SIMULTANEOUSLY, THE FORM OF
    ! APPROXIMATE HESSIAN MATRIX USED IN COMPUTING REGRESSION DIAGNOSTIC
    ! INFORMATION.  IN ALL CASES, SOME APPROXIMATE FINAL HESSIAN MATRIX
    ! IS OBTAINED, AND ITS INVERSE IS THE COVARIANCE MATRIX ESTIMATE
    ! (WHICH MAY HAVE TO BE SCALED APPROPRIATELY -- THAT IS UP TO YOU).
    ! IF IV(COVREQ) IS AT MOST 2 IN ABSOLUTE VALUE, THEN THE FINAL
    ! HESSIAN APPROXIMATION IS COMPUTED BY FINITE DIFFERENCES -- GRADIENT
    ! DIFFERENCES IF IV(COVREQ) IS NONNEGATIVE, FUNCTION DIFFERENCES
    ! OTHERWISE.  IF (IV(COVREQ)) IS AT LEAST 3 IN ABSOLUTE VALUE, THEN THE
    ! CURRENT GAUSS-NEWTON HESSIAN APPROXIMATION IS TAKEN AS THE FINAL
    ! HESSIAN APPROXIMATION.  FOR SOME PROBLEMS THIS SAVES TIME AND YIELDS
    ! THE SAME OR NEARLY THE SAME HESSIAN APPROXIMATION AS DO FINITE
    ! DIFFERENCES.  FOR OTHER PROBLEMS, THE TWO KINDS OF HESSIAN
    ! APPROXIMATIONS MAY GIVE DECIDEDLY DIFFERENT REGRESSION DIAGNOSTICS AND
    ! COVARIANCE MATRIX ESTIMATES.


    !  ***  GENERAL  ***

    !     CODED BY DAVID M. GAY.

    !+++++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++

    !  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***

    ! DSB NOTE: variables for dummy ditsum and dn3rdp
    INTEGER :: inDummy, outDummy
    EXTERNAL dd7up5,divset, dg2lrd, dn3rdp, dd7tpr, dq7adr, dvsum,  &
    ! dg7lit_m,ditsum, dl7nvr, dl7itv, dl7ivm,dl7srt, dl7sqr,  &
    ditsum, dl7nvr, dl7itv, dl7ivm,dl7srt, dl7sqr,  &
    dl7svx, dl7svn, dl7tsq,dl7vml,do7prd,dv2axy,dv7cpy, dv7scl, dv7scp
    real(kind = c_double) :: dd7tpr, dl7svx, dl7svn,dvsum

    ! DD7UP5...  UPDATES SCALE VECTOR D.
    ! DIVSET.... PROVIDES DEFAULT IV AND V INPUT COMPONENTS.
    ! DG2LRD.... COMPUTES REGRESSION DIAGNOSTIC.
    ! DN3RDP... PRINTS REGRESSION DIAGNOSTIC.
    ! DD7TPR... COMPUTES INNER PRODUCT OF TWO VECTORS.
    ! DQ7ADR.... ADDS ROWS TO QR FACTORIZATION.
    ! DVSUM..... RETURNS SUM OF ELEMENTS OF A VECTOR.
    ! dg7lit_m.... PERFORMS BASIC MINIMIZATION ALGORITHM.
    ! DITSUM.... PRINTS ITERATION SUMMARY, INFO ABOUT INITIAL AND FINAL X.
    ! DL7NVR... INVERTS COMPACTLY STORED TRIANGULAR MATRIX.
    ! DL7ITV... MULTIPLIES INVERSE TRANSPOSE OF LOWER TRIANGLE TIMES VECTOR.
    ! DL7IVM... APPLY INVERSE OF COMPACT LOWER TRIANG. MATRIX.
    ! DL7SRT.... COMPUTES CHOLESKY FACTOR OF (LOWER TRIANG. OF) SYM. MATRIX.
    ! DL7SQR... COMPUTES L*(L**T) FOR LOWER TRIANG. MATRIX L.
    ! DL7SVX... UNDERESTIMATES LARGEST SINGULAR VALUE OF TRIANG. MATRIX.
    ! DL7SVN... OVERESTIMATES SMALLEST SINGULAR VALUE OF TRIANG. MATRIX.
    ! DL7TSQ... COMPUTES (L**T)*L FOR LOWER TRIANG. MATRIX L.
    ! DL7VML.... COMPUTES L * V, V = VECTOR, L = LOWER TRIANGULAR MATRIX.
    ! DO7PRD.... ADDS OUTER PRODUCT OF VECTORS TO A MATRIX.
    ! DV2AXY.... ADDS A MULTIPLE OF ONE VECTOR TO ANOTHER.
    ! DV7CPY.... COPIES ONE VECTOR TO ANOTHER.
    ! DV7SCL... MULTIPLIES A VECTOR BY A SCALAR.
    ! DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR.

    !  ***  LOCAL VARIABLES  ***

    LOGICAL :: justg, updatd, zerog
    integer(kind = c_int) :: g1, hn1, i, ii, iv1, j, j1, &
       jtol1, k, l, lh,  &
       need1(2), need2(2),  pmps, ps1, pslen, qtr1,  &
       rmat1, step1, temp1, temp2, temp3, temp4, w, wi, y1
    real(kind = c_double) :: rhmax, rhtol, rho1, rho2, t

    ! DSB NOTE: These were temporarily included for debugging
    ! integer(kind = c_int) :: itest
    ! CHARACTER (LEN=132) :: output_string

    !  ***  SUBSCRIPTS FOR IV AND V  ***

    !  ***  IV SUBSCRIPT VALUES  ***

    integer(kind = c_int), PARAMETER :: cnvcod=55
    integer(kind = c_int), PARAMETER :: covmat=26
    integer(kind = c_int), PARAMETER :: dtype=16
    integer(kind = c_int), PARAMETER :: f0=13
    integer(kind = c_int), PARAMETER :: fdh=74
    integer(kind = c_int), PARAMETER :: g=28
    integer(kind = c_int), PARAMETER :: h=56
    integer(kind = c_int), PARAMETER :: hc=71
    integer(kind = c_int), PARAMETER :: ipivot=76
    integer(kind = c_int), PARAMETER :: ivneed=3
    integer(kind = c_int), PARAMETER :: jcn=66
    integer(kind = c_int), PARAMETER :: jtol=59
    integer(kind = c_int), PARAMETER :: lmat=42
    integer(kind = c_int), PARAMETER :: mode=35
    integer(kind = c_int), PARAMETER :: nextiv=46
    integer(kind = c_int), PARAMETER :: nextv=47
    integer(kind = c_int), PARAMETER :: nfcall=6
    integer(kind = c_int), PARAMETER :: nfcov=52
    integer(kind = c_int), PARAMETER :: nf0=68
    integer(kind = c_int), PARAMETER :: nf1=69
    integer(kind = c_int), PARAMETER :: nfgcal=7
    integer(kind = c_int), PARAMETER :: ngcall=30
    integer(kind = c_int), PARAMETER :: ngcov=53
    integer(kind = c_int), PARAMETER :: perm=58
    integer(kind = c_int), PARAMETER :: qtr=77
    integer(kind = c_int), PARAMETER :: restor=9
    integer(kind = c_int), PARAMETER :: rmat=78
    integer(kind = c_int), PARAMETER :: rdreq=57
    integer(kind = c_int), PARAMETER :: regd=67
    integer(kind = c_int), PARAMETER :: step=40
    integer(kind = c_int), PARAMETER :: toobig=2
    integer(kind = c_int), PARAMETER :: vneed=4
    integer(kind = c_int), PARAMETER :: xnoti=90
    integer(kind = c_int), PARAMETER :: y=48

    !  ***  V SUBSCRIPT VALUES  ***

    integer(kind = c_int), PARAMETER :: dinit=38
    integer(kind = c_int), PARAMETER :: dtinit=39
    integer(kind = c_int), PARAMETER :: d0init=40
    integer(kind = c_int), PARAMETER :: f=10
    integer(kind = c_int), PARAMETER :: rsptol=49
    DOUBLE PRECISION, PARAMETER :: one=1.d+0
    DOUBLE PRECISION, PARAMETER :: zero=0.d+0
    SAVE need1, need2
    DATA need1(1)/1/, need1(2)/0/, need2(1)/2/, need2(2)/0/

    ! DSB NOTE:  The following lines are to stop compiler warnings
    w = 1
    g1 = 1
    temp2 = 1
    !+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
    ! ! output_string = ' Successful entry into DRGLG...'
    ! ! CALL labelpr(output_string,132)
    ! ! output_string = ' Made it this far...'
    ! itest = 0
    ! ! CALL intpr1(output_string,20,itest)
    !      ! output_string = ""
    !      WRITE(output_string,1) "\n", LIV, "\n", LV, "\n", N, "\n", ND,
    !     1  "\n", NN, "\n", P, "\n", PS,"\n", LRHOI, "\n", LRHOR
    !  1     FORMAT(9(A2,1X,I8))
    !        CALL mexPrintf(! output_string)
    ! output_string = ' shit v(20)= '
    ! g1 = 20
    ! rhmax = 200.d0
    ! x(1) = 200.d0
    ! v(g1) = rhmax
    ! iv(1) = 200
    ! CALL dblepr1(output_string,15,v(20))
    i = i_itsum
    i_itsum = i
    lh = p * (p+1) / 2
    !  FIX ME  The following needs to be looked at.
    IF (iv(1) == 0) CALL divset(1, iv, liv, lv, v)
    ! output_string = ' dv7dfl v(20)= '
    ! CALL dblepr1(output_string,15,v(20))
    ps1 = ps + 1
    iv1 = iv(1)
    ! itest = itest + 1
    ! ! CALL intpr1(output_string,20,itest)
    IF (iv1 > 2) GO TO 10
    w = iv(y) + p
    iv(restor) = 0
    i = iv1 + 2
    IF (iv(toobig) == 0) THEN
    SELECT CASE (i)
    CASE (    1)
    GO TO 120
    CASE (    2)
    GO TO  110
    CASE (    3)
    GO TO  110
    CASE (    4)
    GO TO  130
    END SELECT
    END IF
    v(f) = v(f0)
    IF (i /= 3) iv(1) = 2
    GO TO 40

    !  ***  FRESH START OR RESTART -- CHECK INPUT INTEGERS  ***

    10   IF (nd < ps) GO TO 360
    IF (ps > p) GO TO 360
    IF (ps <= 0) GO TO 360
    IF (n <= 0) GO TO 360
    IF (iv1 == 14) GO TO 30
    IF (iv1 > 16) GO TO 420
    IF (iv1 < 12) GO TO 40
    IF (iv1 == 12) iv(1) = 13
    IF (iv(1) /= 13) GO TO 20
    iv(ivneed) = iv(ivneed) + p
    iv(vneed) = iv(vneed) + p*(p+13)/2 + 2*n + 4*ps

    !     *** ADJUST IV(PERM) TO MAKE ROOM FOR IV INPUT COMPONENTS
    !     *** NEEDED WHEN IV(RDREQ) IS 4 OR 5...
    i = xnoti + 1
    IF (iv(perm) < i) iv(perm) = i

    20   CALL dg7lit_m(d, x, iv, liv, lv, p, ps, v, x, x, i_itsum)
    ! 20   CALL dg7lit_m(d, x, iv, liv, lv, p, ps, v, x, x)
    IF (iv(1) /= 14) GO TO 999

    !  ***  STORAGE ALLOCATION  ***
    ! ! output_string = ' DRGLG Storage Allocation'
    ! ! CALL labelpr(output_string,132)
    iv(ipivot) = iv(nextiv)
    iv(nextiv) = iv(ipivot) + p
    iv(y) = iv(nextv)
    iv(g) = iv(y) + p + n
    iv(rmat) = iv(g) + p + 4*ps
    iv(qtr) = iv(rmat) + lh
    iv(jtol) = iv(qtr) + p + n
    iv(jcn) = iv(jtol) + 2*p
    iv(nextv) = iv(jcn) + p
    IF (iv1 == 13) GO TO 999

    30   jtol1 = iv(jtol)
    ! ! output_string = ' Reached line 30 in DRGLG...'
    ! ! CALL labelpr(output_string,132)
    IF (v(dinit) >= zero) CALL dv7scp(p, d, v(dinit))
    IF (v(dtinit) > zero) CALL dv7scp(p, v(jtol1), v(dtinit))
    i = jtol1 + p
    IF (v(d0init) > zero) CALL dv7scp(p, v(i), v(d0init))
    iv(nf0) = 0
    iv(nf1) = 0

    40   CONTINUE
    ! ! output_string = ' DRGLG; have arrived at line 40.'
    ! ! CALL labelpr(output_string,132)
    g1 = iv(g)
    y1 = iv(y)
    ! ! CALL labelpr(output_string,132)
    ! ! output_string = ' Value of iv(1) is: '
    ! CALL intpr1(output_string,20,iv(1))
    ! CALL intpr1(output_string,20,g1)
    ! CALL intpr1(output_string,20,y1)
    ! CALL dg7lit_m(d, v(g1), iv, liv, lv, p, ps, v, x, v(y1))
    CALL dg7lit_m(d, v(g1), iv, liv, lv, p, ps, v, x, v(y1), i_itsum)
    ! output_string = ' Return from dg7lit at line 40 in DRGLG...'
    ! CALL labelpr(output_string,132)
    ! output_string = ' Value of iv(1) is: '
    ! CALL intpr1(output_string,20,iv(1))
    IF (iv(1) - 2 < 0) THEN
    GO TO    50
    ELSE IF (iv(1) - 2 == 0) THEN
    GO TO    60
    ELSE
    GO TO   380
    END IF

    50   v(f) = zero
    IF (iv(nf1) == 0) GO TO 999
    IF (iv(restor) /= 2) GO TO 999
    iv(nf0) = iv(nf1)
    CALL dv7cpy(n, rd, r)
    iv(regd) = 0
    GO TO 999

    60   IF (iv(mode) > 0) GO TO 370
    CALL dv7scp(p, v(g1), zero)
    rmat1 = IABS(iv(rmat))
    qtr1 = IABS(iv(qtr))
    CALL dv7scp(ps, v(qtr1), zero)
    iv(regd) = 0
    CALL dv7scp(ps, v(y1), zero)
    CALL dv7scp(lh, v(rmat1), zero)
    IF (iv(restor) /= 3) GO TO 70
    CALL dv7cpy(n, r, rd)
    iv(nf1) = iv(nf0)
    70   CONTINUE
    ! output_string = ' DRGLG; have arrived at line 70.'
    ! CALL labelpr(output_string,132)
    ! output_string = ' r(1) before = '
    ! CALL dblepr1(output_string,15,r(1))
    CALL rho(need2, t, n, iv(nfgcal), x(ps1), r, rd, rhoi, rhor, v(w))
    ! output_string = ' r(1) after  = '
    ! CALL dblepr1(output_string,15,r(1))
    ! output_string = ' iv(nfgcal)  = '
    ! CALL intpr1(output_string,15,iv(nfgcal))
    ! output_string = ' iv(mode)    = '
    ! CALL intpr1(output_string,15,iv(mode))
    IF (iv(nfgcal) > 0) GO TO 90
    80      iv(toobig) = 1
    GO TO 40
    90   IF (iv(mode) < 0) GO TO 999
    DO  i = 1, n
    CALL dv2axy(ps, v(y1), r(i), dr(1,i), v(y1))
    END DO
    GO TO 999

    !  ***  COMPUTE F(X)  ***

    110  i = iv(nfcall)
    ! output_string = ' Attempting to compute F(x) at 110 in DRGLG...'
    ! CALL labelpr(output_string,132)
    ! output_string = ' First values of r: '
    ! CALL dblepr1(output_string,20,r(19))
    need1(2) = iv(nfgcal)
    CALL rho(need1, v(f), n, i, x(ps1), r, rd, rhoi, rhor, v(w))
    !      ! output_string = ""
    !      WRITE(output_string,2) "\n", V(F)
    !   2  FORMAT(A2,1X,D13.6)
    !        CALL mexPrintf(! output_string)
    iv(nf1) = i
    ! output_string = ' Returned from call to rho with i = '
    ! CALL intpr1(output_string,36,i)
    ! output_string = ' Value of F = '
    ! CALL dblepr1(output_string,14,v(f))
    IF (i <= 0) GO TO 80
    GO TO 40

    !  ***  COMPUTE GRADIENT INFORMATION FOR FINITE-DIFFERENCE HESSIAN  ***

    120  CONTINUE
    iv(1) = 2
    justg = .true.
    i = iv(nfcall)
    CALL rho(need1, t, n, i, x(ps1), r, rd, rhoi, rhor, v(w))
    IF (i <= 0) GO TO 80
    CALL rho(need2, t, n, i, x(ps1), r, rd, rhoi, rhor, v(w))
    IF (i <= 0) GO TO 80
    GO TO 250

    !  ***  PREPARE TO COMPUTE GRADIENT INFORMATION WHILE ITERATING  ***

    130  CONTINUE
       ! output_string = ' Entry into DRGLG with iv(1) = 2; at line 130.'
       ! CALL labelpr(output_string,132)
    justg = .false.
    g1 = iv(g)

    !  ***  DECIDE WHETHER TO UPDATE D BELOW  ***
    i = iv(dtype)
    updatd = .false.
      ! output_string = ' Confirming iv(dytpe) = '
      ! CALL intpr1(output_string,24,i)
    IF (i <= 0) GO TO 140
    IF (i == 1 .OR. iv(mode) < 0) updatd = .true.

    !  ***  COMPUTE RMAT AND QTR  ***

    140  CONTINUE
    ! output_string = ' Entry into DRGLG with iv(1) = 2; at line 140.'
       ! CALL labelpr(output_string,132)
    qtr1 = IABS(iv(qtr))
      ! output_string = ' Confirming iv(qtr) =   '
      ! CALL intpr1(output_string,24,qtr1)
    rmat1 = IABS(iv(rmat))
    iv(rmat) = rmat1
      ! output_string = ' Confirming iv(rmat) =  '
      ! CALL intpr1(output_string,24,rmat1)
    iv(hc) = 0
    iv(nf0) = 0
    iv(nf1) = 0
      ! output_string = ' Confirming iv(mode) =  '
      ! CALL intpr1(output_string,24,iv(mode))
    IF (iv(mode) < 0) GO TO 160

    !  ***  ADJUST Y  ***
       ! output_string = ' Entry into DRGLG with iv(1) = 2; adjust y.'
       ! CALL labelpr(output_string,132)
    y1 = iv(y)
    wi = w
    step1 = iv(step)
    DO  i = 1, n
    t = v(wi) - rd(i)
    wi = wi + 1
    IF (t /= zero) CALL dv2axy(ps, v(y1),  &
    t*dd7tpr(ps,v(step1),dr(1,i)), dr(1,i), v(y1))
    END DO


    !  ***  CHECK FOR NEGATIVE W COMPONENTS  ***

    160  CONTINUE
    ! output_string = ' Entry into DRGLG with iv(1) = 2; line 160.'
    ! CALL labelpr(output_string,132)
    j1 = w + n - 1
    DO  wi = w, j1
    IF (v(wi) < zero) GO TO 240
    END DO

    !  ***  W IS NONNEGATIVE.  COMPUTE QR FACTORIZATION  ***
    !  ***  AND, IF NECESSARY, USE SEMINORMAL EQUATIONS  ***
    ! output_string = ' Entry into DRGLG with iv(1) = 2; W is nonnegative.'
    ! CALL labelpr(output_string,132)

    rhmax = zero
    rhtol = v(rsptol)
    temp1 = g1 + p
    zerog = .true.
    wi = w
    DO  i = 1, n
    rho1 = r(i)
    rho2 = v(wi)
      IF (i .EQ. 0) THEN
        ! output_string = ' i = '
        ! CALL intpr1(output_string,5,i)
        ! output_string = ' r(i) = '
        ! CALL dblepr1(output_string,8,rho1)
        ! output_string = ' V(wi) = '
        ! CALL dblepr1(output_string,9,rho2)
      END IF
    wi = wi + 1
    t =  SQRT(rho2)
    IF (rhmax < rho2) rhmax = rho2
    IF (rho2 > rhtol*rhmax) GO TO 180
    !           *** SEMINORMAL EQUATIONS ***
    CALL dv2axy(ps, v(g1), rho1, dr(1,i), v(g1))
    rho1 = zero
    zerog = .false.
    GO TO 190
    180     rho1 =  rho1 / t
    !        *** QR ACCUMULATION ***
    190     CALL dv7scl(ps, v(temp1), t, dr(1,i))
    CALL dq7adr(ps, v(qtr1), v(rmat1), v(temp1), rho1)
    END DO


    !  ***  COMPUTE G FROM RMAT AND QTR  ***
    ! output_string = ' Entry into DRGLG with iv(1) = 2; COMPUTE G from RMAT and QTR.'
    ! CALL labelpr(output_string,132)
    temp2 = temp1 + ps
    CALL dl7vml(ps, v(temp1), v(rmat1), v(qtr1))
    IF (zerog) GO TO 220
    iv(qtr) = -qtr1
    IF (dl7svx(ps, v(rmat1), v(temp2), v(temp2)) * rhtol >=  &
    dl7svn(ps, v(rmat1), v(temp2), v(temp2))) GO TO 230
    CALL dl7ivm(ps, v(temp2), v(rmat1), v(g1))

    !        *** SEMINORMAL EQUATIONS CORRECTION OF BJOERCK --
    !        *** ONE CYCLE OF ITERATIVE REFINEMENT...
    ! output_string = ' Entry into DRGLG with iv(1) = 2; SEMINORMAL.'
    ! CALL labelpr(output_string,132)

    temp3 = temp2 + ps
    temp4 = temp3 + ps
    CALL dl7itv(ps, v(temp3), v(rmat1), v(temp2))
    CALL dv7scp(ps, v(temp4), zero)
    rhmax = zero
    wi = w
    DO  i = 1, n
    rho2 = v(wi)
    wi = wi + 1
    IF (rhmax < rho2) rhmax = rho2
    rho1 = zero
    IF (rho2 <= rhtol*rhmax) rho1 = r(i)
    t = rho1 - rho2*dd7tpr(ps, v(temp3), dr(1,i))
    CALL dv2axy(ps, v(temp4), t, dr(1,i), v(temp4))
    END DO
    CALL dl7ivm(ps, v(temp3), v(rmat1), v(temp4))
    CALL dv2axy(ps, v(temp2), one, v(temp3), v(temp2))
    CALL dv2axy(ps, v(qtr1), one, v(temp2), v(qtr1))

    220     CONTINUE
    ! output_string = ' Entry into DRGLG with iv(1) = 2; line 220.'
    ! CALL labelpr(output_string,132)
    iv(qtr) = qtr1
    230  CALL dv2axy(ps, v(g1), one, v(temp1), v(g1))
    ! output_string = ' g(1) = '
    ! CALL dblepr1(output_string,8,v(g1))
    IF (ps >= p) GO TO 350
    GO TO 270

    !  ***  INDEFINITE GN HESSIAN...  ***

    240  iv(rmat) = -rmat1
    iv(hc) = rmat1
    CALL do7prd(n, lh, ps, v(rmat1), v(w), dr, dr)

    !  ***  COMPUTE GRADIENT  ***

    250  g1 = iv(g)
    CALL dv7scp(p, v(g1), zero)
    DO  i = 1, n
    CALL dv2axy(ps, v(g1), r(i), dr(1,i), v(g1))
    END DO
    IF (ps >= p) GO TO 350

    !  ***  COMPUTE GRADIENT COMPONENTS OF NUISANCE PARAMETERS ***

    270  k = p - ps
    j1 = 1
    g1 = g1 + ps
    DO  j = 1, k
    j1 = j1 + nn
    v(g1) =dvsum(n, r(j1))
    g1 = g1 + 1
    END DO
    IF (justg) GO TO 390

    !  ***  COMPUTE HESSIAN COMPONENTS OF NUISANCE PARAMETERS  ***

    i = ps*ps1/2
    pslen = p*(p+1)/2 - i
    hn1 = rmat1 + i
    CALL dv7scp(pslen, v(hn1), zero)
    pmps = p - ps
    k = hn1
    j1 = 1
    DO  ii = 1, pmps
    j1 = j1 + nn
    j = j1
    DO  i = 1, n
    CALL dv2axy(ps, v(k), rd(j), dr(1,i), v(k))
    j = j + 1
    END DO
    k = k + ps
    DO  i = 1, ii
    j1 = j1 + nn
    v(k) =dvsum(n, rd(j1))
    k = k + 1
    END DO
    END DO
    IF (iv(rmat) <= 0) GO TO 350
    j = iv(lmat)
    CALL dv7cpy(pslen, v(j), v(hn1))
    IF (dl7svn(ps, v(rmat1), v(temp2), v(temp2)) <= zero) GO TO 320
    CALL dl7srt(ps1, p, v(rmat1), v(rmat1), i)
    IF (i <= 0) GO TO 330

    !  *** HESSIAN IS NOT POSITIVE DEFINITE ***

    320  CALL dl7sqr(ps, v(rmat1), v(rmat1))
    CALL dv7cpy(pslen, v(hn1), v(j))
    iv(hc) = rmat1
    iv(rmat) = -rmat1
    GO TO 350

    !  *** NUISANCE PARS LEAVE HESSIAN POS. DEF.  GET REST OF QTR ***

    330  j = qtr1 + ps
    g1 = iv(g) + ps
    DO  i = ps1, p
    t = dd7tpr(i-1, v(hn1), v(qtr1))
    hn1 = hn1 + i
    v(j) = (v(g1) - t) / v(hn1-1)
    j = j + 1
    g1 = g1 + 1
    END DO

    350  CONTINUE
    ! output_string = ' Entry into DRGLG with iv(1) = 2; line 350.'
    ! CALL labelpr(output_string,132)
    IF (justg) GO TO 390
    IF (updatd) CALL dd7up5(d, iv, liv, lv, p, ps, v)
    GO TO 40

    !  ***  MISC. DETAILS  ***

    !     ***  BAD N, ND, OR P  ***


    360  iv(1) = 66
    GO TO 420

    !  ***  COVARIANCE OR INITIAL S COMPUTATION  ***

    370  iv(nfcov) = iv(nfcov) + 1
    iv(nfcall) = iv(nfcall) + 1
    iv(nfgcal) = iv(nfcall)
    iv(1) = -1
    GO TO 999

    !  ***  CONVERGENCE OBTAINED -- SEE WHETHER TO COMPUTE COVARIANCE  ***

    380  IF (iv(covmat) /= 0) GO TO 410
    IF (iv(regd) /= 0) GO TO 410

    !     ***  SEE IF CHOLESKY FACTOR OF HESSIAN IS AVAILABLE  ***

    k = iv(fdh)
    IF (k <= 0) GO TO 400
    IF (iv(rdreq) <= 0) GO TO 410

    !     ***  COMPUTE REGRESSION DIAGNOSTICS AND DEFAULT COVARIANCE IF
    !          DESIRED  ***

    iv(mode) = p + 1
    iv(ngcall) = iv(ngcall) + 1
    iv(ngcov) = iv(ngcov) + 1
    iv(cnvcod) = iv(1)
    iv(nfcov) = iv(nfcov) + 1
    iv(nfcall) = iv(nfcall) + 1
    iv(nfgcal) = iv(nfcall)
    iv(1) = -1
    GO TO 999

    390  CONTINUE
    ! output_string = ' Entry into DRGLG with iv(1) = 2; line 390.'
    ! CALL labelpr(output_string,132)
    IF (iv(mode) <= p) GO TO 40
    !     *** SAVE RD IN W FOR POSSIBLE USE IN OTHER DIAGNOSTICS ***
    CALL dv7cpy(n, v(w), rd)
    !     *** OVERWRITE RD WITH REGRESSION DIAGNOSTICS ***
    l = iv(lmat)
    i = iv(jcn)
    step1 = iv(step)
    CALL dg2lrd(dr, iv, v(l), lh, liv, lv, nd, n, p, ps, r, rd,  &
    rhoi, rhor, v, v(step1), x, v(i))
    iv(1) = iv(cnvcod)
    iv(cnvcod) = 0
    IF (MOD(iv(rdreq),2) == 0) GO TO 410

    !        *** FINISH COVARIANCE COMPUTATION ***

    i = IABS(iv(h))
    iv(fdh) = 0
    CALL dl7nvr(p, v(i), v(l))
    CALL dl7tsq(p, v(i), v(i))
    iv(covmat) = i
    GO TO 410

    !  ***  COME HERE FOR INDEFINITE FINITE-DIFFERENCE HESSIAN  ***

    400  iv(covmat) = k
    iv(regd) = k

    !  ***  PRINT SUMMARY OF FINAL ITERATION AND OTHER REQUESTED ITEMS  ***

    410  g1 = iv(g)
    ! 420  CALL ditsum(d, v(g1), iv, liv, lv, p, v, x)
    i = i_itsum
    i_itsum = i
    420 CALL ditsum(inDummy,outDummy)
    IF (iv(1) <= 6 .AND. iv(rdreq) > 0)  &
      CALL dn3rdp(inDummy,outDummy)
    ! CALL dn3rdp(iv, liv, lv, n, p, rd, rhoi, rhor, v)

    999  RETURN
    !  ***  LAST LINE OF  DRGLG FOLLOWS  ***
    END SUBROUTINE  drglg

 SUBROUTINE divset_f(alg, iv, liv, lv, v) bind(C, name="divset_f_")
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2009-11-02  Time: 10:37:38
    !  ***  SUPPLY ***SOL (VERSION 2.3) DEFAULT VALUES TO IV AND V  ***
    !  ***  ALG = 1 MEANS REGRESSION CONSTANTS.
    !  ***  ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS.

    use, intrinsic :: iso_c_binding
    ! integer(kind = c_int), INTENT(IN)                  :: alg
    ! integer(kind = c_int), INTENT(OUT)                 :: iv(liv)
    ! integer(kind = c_int), INTENT(IN OUT)              :: liv
    ! integer(kind = c_int), INTENT(IN OUT)              :: lv
    ! real(kind = c_double), INTENT(IN OUT)              :: v(lv)

    integer(kind = c_int), intent(in)                    :: alg
    integer(kind = c_int), intent(in), value             :: liv
    integer(kind = c_int), intent(in), value             :: lv
    integer(kind = c_int), intent(inout), dimension(liv) :: iv
    real(kind = c_double), intent(inout), dimension(lv)  :: v

    integer(kind = c_int) :: i7mdcn
    EXTERNAL i7mdcn, dv7dfl
    ! I7MDCN... RETURNS MACHINE-DEPENDENT INTEGER CONSTANTS.
    ! DV7DFL.... PROVIDES DEFAULT VALUES TO V.

    integer(kind = c_int) :: alg1, miv, mv, i_one
    integer(kind = c_int) :: miniv(4), minv(4)

    !  ***  SUBSCRIPTS FOR IV  ***
    !  ***  IV SUBSCRIPT VALUES  ***

    integer(kind = c_int), PARAMETER :: algsav=51
    integer(kind = c_int), PARAMETER :: covprt=14
    integer(kind = c_int), PARAMETER :: covreq=15
    integer(kind = c_int), PARAMETER :: dradpr=101
    integer(kind = c_int), PARAMETER :: dtype=16
    integer(kind = c_int), PARAMETER :: hc=71
    integer(kind = c_int), PARAMETER :: ierr=75
    integer(kind = c_int), PARAMETER :: inith=25
    integer(kind = c_int), PARAMETER :: inits=25
    integer(kind = c_int), PARAMETER :: ipivot=76
    integer(kind = c_int), PARAMETER :: ivneed=3
    integer(kind = c_int), PARAMETER :: lastiv=44
    integer(kind = c_int), PARAMETER :: lastv=45
    integer(kind = c_int), PARAMETER :: lmat=42
    integer(kind = c_int), PARAMETER :: mxfcal=17
    integer(kind = c_int), PARAMETER :: mxiter=18
    integer(kind = c_int), PARAMETER :: nfcov=52
    integer(kind = c_int), PARAMETER :: ngcov=53
    integer(kind = c_int), PARAMETER :: nvdflt=50
    integer(kind = c_int), PARAMETER :: nvsave=9
    integer(kind = c_int), PARAMETER :: outlev=19
    integer(kind = c_int), PARAMETER :: parprt=20
    integer(kind = c_int), PARAMETER :: parsav=49
    integer(kind = c_int), PARAMETER :: perm=58
    integer(kind = c_int), PARAMETER :: prunit=21
    integer(kind = c_int), PARAMETER :: qrtyp=80
    integer(kind = c_int), PARAMETER :: rdreq=57
    integer(kind = c_int), PARAMETER :: rmat=78
    integer(kind = c_int), PARAMETER :: solprt=22
    integer(kind = c_int), PARAMETER :: statpr=23
    integer(kind = c_int), PARAMETER :: vneed=4
    integer(kind = c_int), PARAMETER :: vsave=60
    integer(kind = c_int), PARAMETER :: x0prt=24

    DATA miniv(1)/82/, miniv(2)/59/, miniv(3)/103/, miniv(4)/103/,  &
        minv(1)/98/, minv(2)/71/, minv(3)/101/, minv(4)/85/

    !-------------------------------  BODY  --------------------------------
    i_one = 1
    ! IF (prunit <= liv) iv(prunit) = i7mdcn(1)
    IF (prunit <= liv) iv(prunit) = i7mdcn(i_one)
    IF (algsav <= liv) iv(algsav) = alg
    IF (alg < 1 .OR. alg > 4) GO TO 40
    miv = miniv(alg)
    IF (liv < miv) GO TO 20
    mv = minv(alg)
    IF (lv < mv) GO TO 30
    alg1 = MOD(alg-1,2) + 1
    CALL dv7dfl(alg1, lv, v)
    iv(1) = 12
    IF (alg > 2) iv(dradpr) = 1
    iv(ivneed) = 0
    iv(lastiv) = miv
    iv(lastv) = mv
    iv(lmat) = mv + 1
    iv(mxfcal) = 200
    iv(mxiter) = 150
    iv(outlev) = 1
    iv(parprt) = 1
    iv(perm) = miv + 1
    iv(solprt) = 1
    iv(statpr) = 1
    iv(vneed) = 0
    iv(x0prt) = 1

    IF (alg1 >= 2) GO TO 10

    !  ***  REGRESSION  VALUES

    iv(covprt) = 3
    iv(covreq) = 1
    iv(dtype) = 1
    iv(hc) = 0
    iv(ierr) = 0
    iv(inits) = 0
    iv(ipivot) = 0
    iv(nvdflt) = 32
    iv(vsave) = 58
    IF (alg > 2) iv(vsave) = iv(vsave) + 3
    iv(parsav) = iv(vsave) + nvsave
    iv(qrtyp) = 1
    iv(rdreq) = 3
    iv(rmat) = 0
    GO TO 999

    !  ***  GENERAL OPTIMIZATION VALUES

    10   iv(dtype) = 0
    iv(inith) = 1
    iv(nfcov) = 0
    iv(ngcov) = 0
    iv(nvdflt) = 25
    iv(parsav) = 47
    IF (alg > 2) iv(parsav) = 61
    GO TO 999

    20   iv(1) = 15
    GO TO 999

    30   iv(1) = 16
    GO TO 999

    40   iv(1) = 67

    999  RETURN
    !  ***  LAST LINE OF DIVSET FOLLOWS  ***
    END SUBROUTINE divset_f

    SUBROUTINE dg7lit_m(d, g, iv, liv, lv, p, ps, v, x, y, i_itsum)

    ! SUBROUTINE dg7lit_m(d, g, iv, liv, lv, p, ps, v, x, y) bind(C, name="dg7lit_m_")

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2021-07-16  Time: 17:05:36

    !  ***  CARRY OUT NL2SOL-LIKE ITERATIONS FOR GENERALIZED LINEAR   ***
    !  ***  REGRESSION PROBLEMS (AND OTHERS OF SIMILAR STRUCTURE)     ***

    !  ***  PARAMETER DECLARATIONS  ***

    use, intrinsic :: iso_c_binding
    implicit none

    integer(kind = c_int), INTENT(IN)             :: liv
    integer(kind = c_int), INTENT(IN)             :: lv
    integer(kind = c_int), INTENT(IN)             :: p
    integer(kind = c_int), INTENT(IN)             :: ps
    real(kind = c_double), INTENT(IN OUT)         :: d(p)
    real(kind = c_double), INTENT(IN OUT)         :: g(p)
    integer(kind = c_int), INTENT(IN OUT)         :: iv(liv)
    real(kind = c_double), INTENT(IN OUT)         :: v(lv)
    real(kind = c_double), INTENT(IN OUT)         :: x(p)
    real(kind = c_double), INTENT(IN OUT)         :: y(p)

    integer(kind = c_int), INTENT(IN OUT)         :: i_itsum

    !--------------------------  PARAMETER USAGE  --------------------------

    ! D.... SCALE VECTOR.
    ! IV... INTEGER VALUE ARRAY.
    ! LIV.. LENGTH OF IV.  MUST BE AT LEAST 82.
    ! LH... LENGTH OF H = P*(P+1)/2.
    ! LV... LENGTH OF V.  MUST BE AT LEAST P*(3*P + 19)/2 + 7.
    ! G.... GRADIENT AT X (WHEN IV(1) = 2).
    ! P.... NUMBER OF PARAMETERS (COMPONENTS IN X).
    ! PS... NUMBER OF NONZERO ROWS AND COLUMNS IN S.
    ! V.... FLOATING-POINT VALUE ARRAY.
    ! X.... PARAMETER VECTOR.
    ! Y.... PART OF YIELD VECTOR (WHEN IV(1)= 2, SCRATCH OTHERWISE).

    !  ***  DISCUSSION  ***

    !       DG7LIT PERFORMS NL2SOL-LIKE ITERATIONS FOR A VARIETY OF
    !     REGRESSION PROBLEMS THAT ARE SIMILAR TO NONLINEAR LEAST-SQUARES
    !     IN THAT THE HESSIAN IS THE SUM OF TWO TERMS, A READILY-COMPUTED
    !     FIRST-ORDER TERM AND A SECOND-ORDER TERM.  THE CALLER SUPPLIES
    !     THE FIRST-ORDER TERM OF THE HESSIAN IN HC (LOWER TRIANGLE, STORED
    !     COMPACTLY BY ROWS IN V, STARTING AT IV(HC)), AND DG7LIT BUILDS AN
    !     APPROXIMATION, S, TO THE SECOND-ORDER TERM.  THE CALLER ALSO
    !     PROVIDES THE FUNCTION VALUE, GRADIENT, AND PART OF THE YIELD
    !     VECTOR USED IN UPDATING S. DG7LIT DECIDES DYNAMICALLY WHETHER OR
    !     NOT TO USE S WHEN CHOOSING THE NEXT STEP TO TRY...  THE HESSIAN
    !     APPROXIMATION USED IS EITHER HC ALONE (GAUSS-NEWTON MODEL) OR
    !     HC + S (AUGMENTED MODEL).

    !        IF PS .LT. P, THEN ROWS AND COLUMNS PS+1...P OF S ARE KEPT
    !     CONSTANT.  THEY WILL BE ZERO UNLESS THE CALLER SETS IV(INITS) TO
    !     1 OR 2 AND SUPPLIES NONZERO VALUES FOR THEM, OR THE CALLER SETS
    !     IV(INITS) TO 3 OR 4 AND THE FINITE-DIFFERENCE INITIAL S THEN
    !     COMPUTED HAS NONZERO VALUES IN THESE ROWS.

    !        IF IV(INITS) IS 3 OR 4, THEN THE INITIAL S IS COMPUTED BY
    !     FINITE DIFFERENCES.  3 MEANS USE FUNCTION DIFFERENCES, 4 MEANS
    !     USE GRADIENT DIFFERENCES.  FINITE DIFFERENCING IS DONE THE SAME
    !     WAY AS IN COMPUTING A COVARIANCE MATRIX (WITH IV(COVREQ) = -1, -2,
    !     1, OR 2).

    !        FOR UPDATING S,DG7LIT ASSUMES THAT THE GRADIENT HAS THE FORM
    !     OF A SUM OVER I OF RHO(I,X)*GRAD(R(I,X)), WHERE GRAD DENOTES THE
    !     GRADIENT WITH RESPECT TO X.  THE TRUE SECOND-ORDER TERM THEN IS
    !     THE SUM OVER I OF RHO(I,X)*HESSIAN(R(I,X)).  IF X = X0 + STEP,
    !     THEN WE WISH TO UPDATE S SO THAT S*STEP IS THE SUM OVER I OF
    !     RHO(I,X)*(GRAD(R(I,X)) - GRAD(R(I,X0))).  THE CALLER MUST SUPPLY
    !     PART OF THIS IN Y, NAMELY THE SUM OVER I OF
    !     RHO(I,X)*GRAD(R(I,X0)), WHEN CALLING DG7LIT WITH IV(1) = 2 AND
    !     IV(MODE) = 0 (WHERE MODE = 38).  G THEN CONTANS THE OTHER PART,
    !     SO THAT THE DESIRED YIELD VECTOR IS G - Y.  IF PS .LT. P, THEN
    !     THE ABOVE DISCUSSION APPLIES ONLY TO THE FIRST PS COMPONENTS OF
    !     GRAD(R(I,X)), STEP, AND Y.

    !        PARAMETERS IV, P, V, AND X ARE THE SAME AS THE CORRESPONDING
    !     ONES TO NL2SOL (WHICH SEE), EXCEPT THAT V CAN BE SHORTER
    !     (SINCE THE PART OF V THAT NL2SOL USES FOR STORING D, J, AND R IS
    !     NOT NEEDED).  MOREOVER, COMPARED WITH NL2SOL, IV(1) MAY HAVE THE
    !     TWO ADDITIONAL OUTPUT VALUES 1 AND 2, WHICH ARE EXPLAINED BELOW,
    !     AS IS THE USE OF IV(TOOBIG) AND IV(NFGCAL).  THE VALUES IV(D),
    !     IV(J), AND IV(R), WHICH ARE OUTPUT VALUES FROM NL2SOL (AND
    !     NL2SNO), ARE NOT REFERENCED BY DG7LIT OR THE SUBROUTINES IT CALLS.

    !        WHEN DG7LIT IS FIRST CALLED, I.E., WHEN DG7LIT IS CALLED WITH
    !     IV(1) = 0 OR 12, V(F), G, AND HC NEED NOT BE INITIALIZED.  TO
    !     OBTAIN THESE STARTING VALUES,DG7LIT RETURNS FIRST WITH IV(1) = 1,
    !     THEN WITH IV(1) = 2, WITH IV(MODE) = -1 IN BOTH CASES.  ON
    !     SUBSEQUENT RETURNS WITH IV(1) = 2, IV(MODE) = 0 IMPLIES THAT
    !     Y MUST ALSO BE SUPPLIED.  (NOTE THAT Y IS USED FOR SCRATCH -- ITS
    !     INPUT CONTENTS ARE LOST.  BY CONTRAST, HC IS NEVER CHANGED.)
    !     ONCE CONVERGENCE HAS BEEN OBTAINED, IV(RDREQ) AND IV(COVREQ) MAY
    !     IMPLY THAT A FINITE-DIFFERENCE HESSIAN SHOULD BE COMPUTED FOR USE
    !     IN COMPUTING A COVARIANCE MATRIX.  IN THIS CASE DG7LIT WILL MAKE A
    !     NUMBER OF RETURNS WITH IV(1) = 1 OR 2 AND IV(MODE) POSITIVE.
    !     WHEN IV(MODE) IS POSITIVE, Y SHOULD NOT BE CHANGED.

    ! IV(1) = 1 MEANS THE CALLER SHOULD SET V(F) (I.E., V(10)) TO F(X), THE
    !             FUNCTION VALUE AT X, AND CALL DG7LIT AGAIN, HAVING CHANGED
    !             NONE OF THE OTHER PARAMETERS.  AN EXCEPTION OCCURS IF F(X)
    !             CANNOT BE EVALUATED (E.G. IF OVERFLOW WOULD OCCUR), WHICH
    !             MAY HAPPEN BECAUSE OF AN OVERSIZED STEP.  IN THIS CASE
    !             THE CALLER SHOULD SET IV(TOOBIG) = IV(2) TO 1, WHICH WILL
    !             CAUSE DG7LIT TO IGNORE V(F) AND TRY A SMALLER STEP.  NOTE
    !             THAT THE CURRENT FUNCTION EVALUATION COUNT IS AVAILABLE
    !             IN IV(NFCALL) = IV(6).  THIS MAY BE USED TO IDENTIFY
    !             WHICH COPY OF SAVED INFORMATION SHOULD BE USED IN COM-
    !             PUTING G, HC, AND Y THE NEXT TIME DG7LIT RETURNS WITH
    !             IV(1) = 2.  SEE MLPIT FOR AN EXAMPLE OF THIS.
    ! IV(1) = 2 MEANS THE CALLER SHOULD SET G TO G(X), THE GRADIENT OF F AT
    !             X.  THE CALLER SHOULD ALSO SET HC TO THE GAUSS-NEWTON
    !             HESSIAN AT X.  IF IV(MODE) = 0, THEN THE CALLER SHOULD
    !             ALSO COMPUTE THE PART OF THE YIELD VECTOR DESCRIBED ABOVE.
    !             THE CALLER SHOULD THEN CALL DG7LIT AGAIN (WITH IV(1) = 2).
    !             THE CALLER MAY ALSO CHANGE D AT THIS TIME, BUT SHOULD NOT
    !             CHANGE X.  NOTE THAT IV(NFGCAL) = IV(7) CONTAINS THE
    !             VALUE THAT IV(NFCALL) HAD DURING THE RETURN WITH
    !             IV(1) = 1 IN WHICH X HAD THE SAME VALUE AS IT NOW HAS.
    !             IV(NFGCAL) IS EITHER IV(NFCALL) OR IV(NFCALL) - 1.  MLPIT
    !             IS AN EXAMPLE WHERE THIS INFORMATION IS USED.  IF G OR HC
    !             CANNOT BE EVALUATED AT X, THEN THE CALLER MAY SET
    !             IV(TOOBIG) TO 1, IN WHICH CASE DG7LIT WILL RETURN WITH
    !             IV(1) = 15.

    !  ***  GENERAL  ***

    !     CODED BY DAVID M. GAY.
    !     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
    !     SUPPORTED IN PART BY D.O.E. GRANT EX-76-A-01-2295 TO MIT/CCREMS.

    !        (SEE NL2SOL FOR REFERENCES.)

    !+++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++

    !  ***  LOCAL VARIABLES  ***

    integer(kind = c_int) :: dummy, dig1, g01, h1, hc1, i, ipiv1, j, k, l, lmat1,  &
        lstgst, pp1o2, qtr1, rmat1, rstrst, step1, stpmod, s1, temp1, temp2, w1, x01
    real(kind = c_double) :: e, sttsst, t, t1

    integer(kind = c_int) :: alg, liv_1, lv_1

    ! DSB NOTE: Variables for dummy ditsum
    INTEGER :: inDummy, outDummy

    ! DSB NOTE: These were temporarily included for debugging
    ! integer(kind = c_int) :: itest
    ! CHARACTER (LEN=132) :: output_string

    !     ***  CONSTANTS  ***

    !  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***

    LOGICAL :: stopx
    real(kind = c_double) :: dd7tpr, dl7svx, dl7svn, drldst, dr7mdc, dv2nrm
    EXTERNAL da7sst, dd7tpr,df7hes,dg7qts,ditsum, dl7mst,dl7srt,  &
        ! dl7sqr, dl7svx, dl7svn, dl7tvm,dl7vml,dparck_m, drldst,  &
        dl7sqr, dl7svx, dl7svn, dl7tvm,dl7vml, drldst,  &
        dr7mdc, ds7lup, ds7lvm, stopx,dv2axy,dv7cpy, dv7scp, dv2nrm

    ! DA7SST.... ASSESSES CANDIDATE STEP.
    ! DD7TPR... RETURNS INNER PRODUCT OF TWO VECTORS.
    ! DF7HES.... COMPUTE FINITE-DIFFERENCE HESSIAN (FOR COVARIANCE).
    ! DG7QTS.... COMPUTES GOLDFELD-QUANDT-TROTTER STEP (AUGMENTED MODEL).
    ! DITSUM.... PRINTS ITERATION SUMMARY AND INFO ON INITIAL AND FINAL X.
    ! DL7MST... COMPUTES LEVENBERG-MARQUARDT STEP (GAUSS-NEWTON MODEL).
    ! DL7SRT.... COMPUTES CHOLESKY FACTOR OF (LOWER TRIANG. OF) SYM. MATRIX.
    ! DL7SQR... COMPUTES L * L**T FROM LOWER TRIANGULAR MATRIX L.
    ! DL7TVM... COMPUTES L**T * V, V = VECTOR, L = LOWER TRIANGULAR MATRIX.
    ! DL7SVX... ESTIMATES LARGEST SING. VALUE OF LOWER TRIANG. MATRIX.
    ! DL7SVN... ESTIMATES SMALLEST SING. VALUE OF LOWER TRIANG. MATRIX.
    ! DL7VML.... COMPUTES L * V, V = VECTOR, L = LOWER TRIANGULAR MATRIX.
    ! dparck_m.... CHECK VALIDITY OF IV AND V INPUT COMPONENTS.
    ! DRLDST... COMPUTES V(RELDX) = RELATIVE STEP SIZE.
    ! DR7MDC... RETURNS MACHINE-DEPENDENT CONSTANTS.
    ! DS7LUP... PERFORMS QUASI-NEWTON UPDATE ON COMPACTLY STORED LOWER TRI-
    !             ANGLE OF A SYMMETRIC MATRIX.
    ! STOPX.... RETURNS .TRUE. IF THE BREAK KEY HAS BEEN PRESSED.
    ! DV2AXY.... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER.
    ! DV7CPY.... COPIES ONE VECTOR TO ANOTHER.
    ! DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR.
    ! DV2NRM... RETURNS THE 2-NORM OF A VECTOR.

    !  ***  SUBSCRIPTS FOR IV AND V  ***

    !  ***  IV SUBSCRIPT VALUES  ***

    integer(kind = c_int), PARAMETER :: cnvcod=55
    integer(kind = c_int), PARAMETER :: covmat=26
    integer(kind = c_int), PARAMETER :: covreq=15
    integer(kind = c_int), PARAMETER :: dig=37
    integer(kind = c_int), PARAMETER :: fdh=74
    integer(kind = c_int), PARAMETER :: h=56
    integer(kind = c_int), PARAMETER :: hc=71
    integer(kind = c_int), PARAMETER :: ierr=75
    integer(kind = c_int), PARAMETER :: inits=25
    integer(kind = c_int), PARAMETER :: ipivot=76
    integer(kind = c_int), PARAMETER :: irc=29
    integer(kind = c_int), PARAMETER :: kagqt=33
    integer(kind = c_int), PARAMETER :: kalm=34
    integer(kind = c_int), PARAMETER :: lmat=42
    integer(kind = c_int), PARAMETER :: mode=35
    integer(kind = c_int), PARAMETER :: model=5
    integer(kind = c_int), PARAMETER :: mxfcal=17
    integer(kind = c_int), PARAMETER :: mxiter=18
    integer(kind = c_int), PARAMETER :: nextv=47
    integer(kind = c_int), PARAMETER :: nfcall=6
    integer(kind = c_int), PARAMETER :: nfgcal=7
    integer(kind = c_int), PARAMETER :: nfcov=52
    integer(kind = c_int), PARAMETER :: ngcov=53
    integer(kind = c_int), PARAMETER :: ngcall=30
    integer(kind = c_int), PARAMETER :: niter=31
    integer(kind = c_int), PARAMETER :: qtr=77
    integer(kind = c_int), PARAMETER :: radinc=8
    integer(kind = c_int), PARAMETER :: rdreq=57
    integer(kind = c_int), PARAMETER :: regd=67
    integer(kind = c_int), PARAMETER :: restor=9
    integer(kind = c_int), PARAMETER :: rmat=78
    integer(kind = c_int), PARAMETER :: s=62
    integer(kind = c_int), PARAMETER :: step=40
    integer(kind = c_int), PARAMETER :: stglim=11
    integer(kind = c_int), PARAMETER :: stlstg=41
    integer(kind = c_int), PARAMETER :: sused=64
    integer(kind = c_int), PARAMETER :: switch=12
    integer(kind = c_int), PARAMETER :: toobig=2
    integer(kind = c_int), PARAMETER :: vneed=4
    integer(kind = c_int), PARAMETER :: vsave=60
    integer(kind = c_int), PARAMETER :: w=65
    integer(kind = c_int), PARAMETER :: xirc=13
    integer(kind = c_int), PARAMETER :: x0=43

    !  ***  V SUBSCRIPT VALUES  ***

    integer(kind = c_int), PARAMETER :: cosmin=47
    integer(kind = c_int), PARAMETER :: dgnorm=1
    integer(kind = c_int), PARAMETER :: dstnrm=2
    integer(kind = c_int), PARAMETER :: f=10
    integer(kind = c_int), PARAMETER :: fdif=11
    integer(kind = c_int), PARAMETER :: fuzz=45
    integer(kind = c_int), PARAMETER :: f0=13
    integer(kind = c_int), PARAMETER :: gtstep=4
    integer(kind = c_int), PARAMETER :: incfac=23
    integer(kind = c_int), PARAMETER :: lmax0=35
    integer(kind = c_int), PARAMETER :: lmaxs=36
    integer(kind = c_int), PARAMETER :: nvsave=9
    integer(kind = c_int), PARAMETER :: phmxfc=21
    integer(kind = c_int), PARAMETER :: preduc=7
    integer(kind = c_int), PARAMETER :: radfac=16
    integer(kind = c_int), PARAMETER :: radius=8
    integer(kind = c_int), PARAMETER :: rad0=9
    integer(kind = c_int), PARAMETER :: rcond=53
    integer(kind = c_int), PARAMETER :: reldx=17
    integer(kind = c_int), PARAMETER :: size=55
    integer(kind = c_int), PARAMETER :: stppar=5
    integer(kind = c_int), PARAMETER :: tuner4=29
    integer(kind = c_int), PARAMETER :: tuner5=30
    integer(kind = c_int), PARAMETER :: wscale=56

    DOUBLE PRECISION, PARAMETER :: half=0.5D+0
    DOUBLE PRECISION, PARAMETER :: negone=-1.d+0
    DOUBLE PRECISION, PARAMETER :: one=1.d+0
    DOUBLE PRECISION, PARAMETER :: onep2=1.2D+0
    DOUBLE PRECISION, PARAMETER :: zero=0.d+0

    !+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
    alg = 1
    liv_1 = liv
    lv_1  = lv
    i = iv(1)
    IF (i == 1) GO TO 40
    IF (i == 2) GO TO 50
    IF (i == 12 .OR. i == 13) iv(vneed) = iv(vneed) + p*(3*p + 19)/2 + 7

    ! DSB NOTE:  dparck is a major source of WRITE statements

    i = iv(1) - 2
    IF (i > 12) GO TO 999
    SELECT CASE ( i )
      CASE (    1)
        GO TO 290
      CASE (    2)
        GO TO  290
      CASE (    3)
        GO TO  290
      CASE (    4)
        GO TO  290
      CASE (    5)
        GO TO  290
      CASE (    6)
        GO TO  290
      CASE (    7)
        GO TO  170
      CASE (    8)
        GO TO  120
      CASE (    9)
        GO TO  170
      CASE (   10)
        GO TO  10
      CASE (   11)
        GO TO  10
      CASE (   12)
        GO TO  20
    END SELECT

    !  ***  STORAGE ALLOCATION  ***

    10   pp1o2 = p * (p + 1) / 2
    iv(s) = iv(lmat) + pp1o2
    iv(x0) = iv(s) + pp1o2
    iv(step) = iv(x0) + p
    iv(stlstg) = iv(step) + p
    iv(dig) = iv(stlstg) + p
    iv(w) = iv(dig) + p
    iv(h) = iv(w) + 4*p + 7
    iv(nextv) = iv(h) + pp1o2
    IF (iv(1) /= 13) GO TO 20
    iv(1) = 14
    GO TO 999

    !  ***  INITIALIZATION  ***

    20  CONTINUE
    iv(niter) = 0
    iv(nfcall) = 1
    iv(ngcall) = 1
    iv(nfgcal) = 1
    iv(mode) = -1
    iv(stglim) = 2
    iv(toobig) = 0
    iv(cnvcod) = 0
    iv(covmat) = 0
    iv(nfcov) = 0
    iv(ngcov) = 0
    iv(radinc) = 0
    iv(restor) = 0
    iv(fdh) = 0
    v(rad0) = zero
    v(stppar) = zero
    v(radius) = v(lmax0) / (one + v(phmxfc))

    !  ***  SET INITIAL MODEL AND S MATRIX  ***

    iv(model) = 1
    IF (iv(s) < 0) GO TO 999
    IF (iv(inits) > 1) iv(model) = 2
    s1 = iv(s)
    IF (iv(inits) == 0 .OR. iv(inits) > 2) CALL dv7scp(p*(p+1)/2, v(s1), zero)
    iv(1) = 1
    j = iv(ipivot)
    IF (j <= 0) GO TO 999
    DO  i = 1, p
      iv(j) = i
      j = j + 1
    END DO
    GO TO 999

    !  ***  NEW FUNCTION VALUE  ***

    40   CONTINUE

    IF (iv(mode) == 0) GO TO 290
    IF (iv(mode) > 0) GO TO 520

    iv(1) = 2
    IF (iv(toobig) == 0) GO TO 999
    iv(1) = 63
    GO TO 999

    !  ***  NEW GRADIENT  ***

    50   iv(kalm) = -1
    iv(kagqt) = -1
    iv(fdh) = 0
    IF (iv(mode) > 0) GO TO 520

    !  ***  MAKE SURE GRADIENT COULD BE COMPUTED  ***

    IF (iv(toobig) == 0) GO TO 60
    iv(1) = 65
    GO TO 999
    60   IF (iv(hc) <= 0 .AND. iv(rmat) <= 0) GO TO 610

    !  ***  COMPUTE  D**-1 * GRADIENT  ***

    dig1 = iv(dig)
    k = dig1
    DO  i = 1, p
      v(k) = g(i) / d(i)
      k = k + 1
    END DO
    v(dgnorm) = dv2nrm(p, v(dig1))

    IF (iv(cnvcod) /= 0) GO TO 510
    IF (iv(mode) == 0) GO TO 440
    iv(mode) = 0
    v(f0) = v(f)
    IF (iv(inits) <= 2) GO TO 100

    !  ***  ARRANGE FOR FINITE-DIFFERENCE INITIAL S  ***

    iv(xirc) = iv(covreq)
    iv(covreq) = -1
    IF (iv(inits) > 3) iv(covreq) = 1
    iv(cnvcod) = 70
    GO TO 530

    !  ***  COME TO NEXT STMT AFTER COMPUTING F.D. HESSIAN FOR INIT. S  ***

    80   iv(cnvcod) = 0
    iv(mode) = 0
    iv(nfcov) = 0
    iv(ngcov) = 0
    iv(covreq) = iv(xirc)
    s1 = iv(s)
    pp1o2 = ps * (ps + 1) / 2
    hc1 = iv(hc)
    IF (hc1 <= 0) GO TO 90
    CALL dv2axy(pp1o2, v(s1), negone, v(hc1), v(h1))
    GO TO 100
    90   rmat1 = iv(rmat)
    CALL dl7sqr(ps, v(s1), v(rmat1))
    CALL dv2axy(pp1o2, v(s1), negone, v(s1), v(h1))
    100  iv(1) = 2


    !-----------------------------  MAIN LOOP  -----------------------------


    !  ***  PRINT ITERATION SUMMARY, CHECK ITERATION LIMIT  ***

    ! 110  CALL ditsum(d, g, iv, liv, lv, p, v, x)
    110  CALL ditsum(inDummy,outDummy)
    ! 110  CONTINUE
    i = i_itsum
    i_itsum = i
    120  k = iv(niter)
    IF (k < iv(mxiter)) GO TO 130
    iv(1) = 10
    GO TO 999
    130  iv(niter) = k + 1

    !  ***  UPDATE RADIUS  ***

    IF (k == 0) GO TO 150
    step1 = iv(step)
    DO  i = 1, p
      v(step1) = d(i) * v(step1)
      step1 = step1 + 1
    END DO
    step1 = iv(step)
    t = v(radfac) * dv2nrm(p, v(step1))
    IF (v(radfac) < one .OR. t > v(radius)) v(radius) = t

    !  ***  INITIALIZE FOR START OF NEXT ITERATION  ***

    150  x01 = iv(x0)
    v(f0) = v(f)
    iv(irc) = 4
    iv(h) = -IABS(iv(h))
    iv(sused) = iv(model)

    !     ***  COPY X TO X0  ***

    CALL dv7cpy(p, v(x01), x)

    !  ***  CHECK STOPX AND FUNCTION EVALUATION LIMIT  ***

    160  IF (.NOT. stopx(dummy)) GO TO 180
    iv(1) = 11
    GO TO 190

    !     ***  COME HERE WHEN RESTARTING AFTER FUNC. EVAL. LIMIT OR STOPX.

    170  IF (v(f) >= v(f0)) GO TO 180
    v(radfac) = one
    k = iv(niter)
    GO TO 130

    180  IF (iv(nfcall) < iv(mxfcal) + iv(nfcov)) GO TO 200
    iv(1) = 9
    190     IF (v(f) >= v(f0)) GO TO 999

    !        ***  IN CASE OF STOPX OR FUNCTION EVALUATION LIMIT WITH
    !        ***  IMPROVED V(F), EVALUATE THE GRADIENT AT X.

    iv(cnvcod) = iv(1)
    GO TO 430

    !. . . . . . . . . . . . .  COMPUTE CANDIDATE STEP  . . . . . . . . . .

    200  step1 = iv(step)
    w1 = iv(w)
    h1 = iv(h)
    t1 = one
    IF (iv(model) == 2) GO TO 210
    t1 = zero

    !        ***  COMPUTE LEVENBERG-MARQUARDT STEP IF POSSIBLE...

    rmat1 = iv(rmat)
    IF (rmat1 <= 0) GO TO 210
    qtr1 = iv(qtr)
    IF (qtr1 <= 0) GO TO 210
    ipiv1 = iv(ipivot)
    CALL dl7mst(d, g, iv(ierr), iv(ipiv1), iv(kalm), p, v(qtr1),  &
        v(rmat1), v(step1), v, v(w1))
    !        *** H IS STORED IN THE END OF W AND HAS JUST BEEN OVERWRITTEN,
    !        *** SO WE MARK IT INVALID...
    iv(h) = -IABS(h1)
    !        *** EVEN IF H WERE STORED ELSEWHERE, IT WOULD BE NECESSARY TO
    !        *** MARK INVALID THE INFORMATION DG7QTS MAY HAVE STORED IN V...
    iv(kagqt) = -1
    GO TO 260

    210  IF (h1 > 0) GO TO 250

    !     ***  SET H TO  D**-1 * (HC + T1*S) * D**-1.  ***

    h1 = -h1
    iv(h) = h1
    iv(fdh) = 0
    j = iv(hc)
    IF (j > 0) GO TO 220
    j = h1
    rmat1 = iv(rmat)
    CALL dl7sqr(p, v(h1), v(rmat1))
    220     s1 = iv(s)
    DO  i = 1, p
      t = one / d(i)
      DO  k = 1, i
        v(h1) = t * (v(j) + t1*v(s1)) / d(k)
        j = j + 1
        h1 = h1 + 1
        s1 = s1 + 1
      END DO
    END DO
    h1 = iv(h)
    iv(kagqt) = -1

    !  ***  COMPUTE ACTUAL GOLDFELD-QUANDT-TROTTER STEP  ***

    250  dig1 = iv(dig)
    lmat1 = iv(lmat)
    CALL dg7qts(d, v(dig1), v(h1), iv(kagqt), v(lmat1), p, v(step1), v, v(w1))
    IF (iv(kalm) > 0) iv(kalm) = 0

    260  IF (iv(irc) /= 6) GO TO 270
    IF (iv(restor) /= 2) GO TO 290
    rstrst = 2
    GO TO 300

    !  ***  CHECK WHETHER EVALUATING F(X0 + STEP) LOOKS WORTHWHILE  ***

    270  iv(toobig) = 0
    IF (v(dstnrm) <= zero) GO TO 290
    IF (iv(irc) /= 5) GO TO 280
    IF (v(radfac) <= one) GO TO 280
    IF (v(preduc) > onep2 * v(fdif)) GO TO 280
    IF (iv(restor) /= 2) GO TO 290
    rstrst = 0
    GO TO 300

    !  ***  COMPUTE F(X0 + STEP)  ***

    280  x01 = iv(x0)
    step1 = iv(step)
    CALL dv2axy(p, x, one, v(step1), v(x01))
    iv(nfcall) = iv(nfcall) + 1
    iv(1) = 1
    GO TO 999

    !. . . . . . . . . . . . .  ASSESS CANDIDATE STEP  . . . . . . . . . . .

    290  rstrst = 3
    300  x01 = iv(x0)
    v(reldx) = drldst(p, d, x, v(x01))
    CALL da7sst(iv, liv, lv, v)
    iv(switch) = 0
    step1 = iv(step)
    lstgst = iv(stlstg)
    i = iv(restor) + 1
    SELECT CASE ( i )
      CASE (    1)
        GO TO 340
      CASE (    2)
        GO TO  310
      CASE (    3)
        GO TO  320
      CASE (    4)
        GO TO  330
    END SELECT
    310  CALL dv7cpy(p, x, v(x01))
    GO TO 340
    320   CALL dv7cpy(p, v(lstgst), v(step1))
    GO TO 340
    330     CALL dv7cpy(p, v(step1), v(lstgst))
    CALL dv2axy(p, x, one, v(step1), v(x01))
    v(reldx) = drldst(p, d, x, v(x01))
    iv(restor) = rstrst

    !  ***  IF NECESSARY, SWITCH MODELS  ***

    340  IF (iv(switch) == 0) GO TO 350
    iv(h) = -IABS(iv(h))
    iv(sused) = iv(sused) + 2
    l = iv(vsave)
    CALL dv7cpy(nvsave, v, v(l))
    350  l = iv(irc) - 4
    stpmod = iv(model)
    IF (l > 0) THEN
       SELECT CASE ( l )
        CASE (    1)
          GO TO 370
        CASE (    2)
          GO TO 380
        CASE (    3)
          GO TO 390
        CASE (    4)
          GO TO 390
        CASE (    5)
          GO TO 390
        CASE (    6)
          GO TO 390
        CASE (    7)
          GO TO 390
        CASE (    8)
          GO TO 390
        CASE (    9)
          GO TO 500
        CASE (   10)
          GO TO 440
      END SELECT
    END IF

    !  ***  DECIDE WHETHER TO CHANGE MODELS  ***

    e = v(preduc) - v(fdif)
    s1 = iv(s)
    CALL ds7lvm(ps, y, v(s1), v(step1))
    sttsst = half * dd7tpr(ps, v(step1), y)
    IF (iv(model) == 1) sttsst = -sttsst
    IF ( ABS(e + sttsst) * v(fuzz) >=  ABS(e)) GO TO 360

    !     ***  SWITCH MODELS  ***

    iv(model) = 3 - iv(model)
    IF (-2 < l) GO TO 400
    iv(h) = -IABS(iv(h))
    iv(sused) = iv(sused) + 2
    l = iv(vsave)
    CALL dv7cpy(nvsave, v(l), v)
    GO TO 160

    360  IF (-3 < l) GO TO 400

    !  ***  RECOMPUTE STEP WITH NEW RADIUS  ***

    370  v(radius) = v(radfac) * v(dstnrm)
    GO TO 160

    !  ***  COMPUTE STEP OF LENGTH V(LMAXS) FOR SINGULAR CONVERGENCE TEST

    380  v(radius) = v(lmaxs)
    GO TO 200

    !  ***  CONVERGENCE OR FALSE CONVERGENCE  ***

    390  iv(cnvcod) = l
    IF (v(f) >= v(f0)) GO TO 510
    IF (iv(xirc) == 14) GO TO 510
    iv(xirc) = 14

    !. . . . . . . . . . . .  PROCESS ACCEPTABLE STEP  . . . . . . . . . . .

    400  iv(covmat) = 0
    iv(regd) = 0

    !  ***  SEE WHETHER TO SET V(RADFAC) BY GRADIENT TESTS  ***

    IF (iv(irc) /= 3) GO TO 430
    step1 = iv(step)
    temp1 = iv(stlstg)
    temp2 = iv(w)

    !     ***  SET  TEMP1 = HESSIAN * STEP  FOR USE IN GRADIENT TESTS  ***

    hc1 = iv(hc)
    IF (hc1 <= 0) GO TO 410
    CALL ds7lvm(p, v(temp1), v(hc1), v(step1))
    GO TO 420
    410     rmat1 = iv(rmat)
    CALL dl7tvm(p, v(temp1), v(rmat1), v(step1))
    CALL dl7vml(p, v(temp1), v(rmat1), v(temp1))

    420     IF (stpmod == 1) GO TO 430
    s1 = iv(s)
    CALL ds7lvm(ps, v(temp2), v(s1), v(step1))
    CALL dv2axy(ps, v(temp1), one, v(temp2), v(temp1))

    !  ***  SAVE OLD GRADIENT AND COMPUTE NEW ONE  ***

    430  iv(ngcall) = iv(ngcall) + 1
    g01 = iv(w)
    CALL dv7cpy(p, v(g01), g)
    iv(1) = 2
    iv(toobig) = 0
    GO TO 999

    !  ***  INITIALIZATIONS -- G0 = G - G0, ETC.  ***

    440  g01 = iv(w)
    CALL dv2axy(p, v(g01), negone, v(g01), g)
    step1 = iv(step)
    temp1 = iv(stlstg)
    temp2 = iv(w)
    IF (iv(irc) /= 3) GO TO 470

    !  ***  SET V(RADFAC) BY GRADIENT TESTS  ***

    !     ***  SET  TEMP1 = D**-1 * (HESSIAN * STEP  +  (G(X0) - G(X)))  ***

    k = temp1
    l = g01
    DO  i = 1, p
      v(k) = (v(k) - v(l)) / d(i)
      k = k + 1
      l = l + 1
    END DO

    !        ***  DO GRADIENT TESTS  ***

    IF (dv2nrm(p, v(temp1)) <= v(dgnorm) * v(tuner4))  GO TO 460
    IF (dd7tpr(p, g, v(step1))  >= v(gtstep) * v(tuner5))  GO TO 470
    460               v(radfac) = v(incfac)

    !  ***  COMPUTE Y VECTOR NEEDED FOR UPDATING S  ***

    470  CALL dv2axy(ps, y, negone, y, g)

    !  ***  DETERMINE SIZING FACTOR V(SIZE)  ***

    !     ***  SET TEMP1 = S * STEP  ***
    s1 = iv(s)
    CALL ds7lvm(ps, v(temp1), v(s1), v(step1))

    t1 =  ABS(dd7tpr(ps, v(step1), v(temp1)))
    t =  ABS(dd7tpr(ps, v(step1), y))
    v(size) = one
    IF (t < t1) v(size) = t / t1

    !  ***  SET G0 TO WCHMTD CHOICE OF FLETCHER AND AL-BAALI  ***

    hc1 = iv(hc)
    IF (hc1 <= 0) GO TO 480
    CALL ds7lvm(ps, v(g01), v(hc1), v(step1))
    GO TO 490

    480  rmat1 = iv(rmat)
    CALL dl7tvm(ps, v(g01), v(rmat1), v(step1))
    CALL dl7vml(ps, v(g01), v(rmat1), v(g01))

    490  CALL dv2axy(ps, v(g01), one, y, v(g01))

    !  ***  UPDATE S  ***

    CALL ds7lup(v(s1), v(cosmin), ps, v(size), v(step1), v(temp1),  &
        v(temp2), v(g01), v(wscale), y)
    iv(1) = 2
    GO TO 110

    !. . . . . . . . . . . . . .  MISC. DETAILS  . . . . . . . . . . . . . .

    !  ***  BAD PARAMETERS TO ASSESS  ***

    500  iv(1) = 64
    GO TO 999


    !  ***  CONVERGENCE OBTAINED -- SEE WHETHER TO COMPUTE COVARIANCE  ***

    510  IF (iv(rdreq) == 0) GO TO 600
    IF (iv(fdh) /= 0) GO TO 600
    IF (iv(cnvcod) >= 7) GO TO 600
    IF (iv(regd) > 0) GO TO 600
    IF (iv(covmat) > 0) GO TO 600
    IF (IABS(iv(covreq)) >= 3) GO TO 560
    IF (iv(restor) == 0) iv(restor) = 2
    GO TO 530

    !  ***  COMPUTE FINITE-DIFFERENCE HESSIAN FOR COMPUTING COVARIANCE  ***

    520  iv(restor) = 0
    530  CALL df7hes(d, g, i, iv, liv, lv, p, v, x)
    SELECT CASE ( i )
      CASE (    1)
        GO TO 540
      CASE (    2)
        GO TO  550
      CASE (    3)
        GO TO  580
    END SELECT
    540  iv(nfcov) = iv(nfcov) + 1
    iv(nfcall) = iv(nfcall) + 1
    iv(1) = 1
    GO TO 999

    550  iv(ngcov) = iv(ngcov) + 1
    iv(ngcall) = iv(ngcall) + 1
    iv(nfgcal) = iv(nfcall) + iv(ngcov)
    iv(1) = 2
    GO TO 999

    560  h1 = IABS(iv(h))
    iv(h) = -h1
    pp1o2 = p * (p + 1) / 2
    rmat1 = iv(rmat)
    IF (rmat1 <= 0) GO TO 570
    lmat1 = iv(lmat)
    CALL dv7cpy(pp1o2, v(lmat1), v(rmat1))
    v(rcond) = zero
    GO TO 590
    570  hc1 = iv(hc)
    iv(fdh) = h1
    CALL dv7cpy(p*(p+1)/2, v(h1), v(hc1))

    !  ***  COMPUTE CHOLESKY FACTOR OF FINITE-DIFFERENCE HESSIAN
    !  ***  FOR USE IN CALLER*S COVARIANCE CALCULATION...

    580  lmat1 = iv(lmat)
    h1 = iv(fdh)
    IF (h1 <= 0) GO TO 600
    IF (iv(cnvcod) == 70) GO TO 80
    CALL dl7srt(1, p, v(lmat1), v(h1), i)
    iv(fdh) = -1
    v(rcond) = zero
    IF (i /= 0) GO TO 600

    590  iv(fdh) = -1
    step1 = iv(step)
    t = dl7svn(p, v(lmat1), v(step1), v(step1))
    IF (t <= zero) GO TO 600
    t = t / dl7svx(p, v(lmat1), v(step1), v(step1))
    IF (t > dr7mdc(4)) iv(fdh) = h1
    v(rcond) = t

    600  iv(mode) = 0
    iv(1) = iv(cnvcod)
    iv(cnvcod) = 0
    GO TO 999

    !  ***  SPECIAL RETURN FOR MISSING HESSIAN INFORMATION -- BOTH
    !  ***  IV(HC) .LE. 0 AND IV(RMAT) .LE. 0

    610  iv(1) = 1400

    999  RETURN

    !  ***  LAST LINE OF DG7LIT FOLLOWS  ***
    END SUBROUTINE dg7lit_m

    SUBROUTINE dparck_m(alg, d, iv, liv, lv, n, v)

    ! bind(C,name="dparck_m_")
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2021-07-16  Time: 12:07:13
    ! Special code created by David Bunch to replace wr to a Fortran
    ! unit to writing to a character string.

    !  ***  CHECK ***SOL (VERSION 2.3) PARAMETERS, PRINT CHANGED VALUES  ***

    !  ***  ALG = 1 FOR REGRESSION, ALG = 2 FOR GENERAL UNCONSTRAINED OPT.

    use, intrinsic :: iso_c_binding
    implicit none

    integer(kind = c_int), INTENT(IN)             :: alg
    integer(kind = c_int), INTENT(IN)             :: liv
    integer(kind = c_int), INTENT(IN)             :: lv
    integer(kind = c_int), INTENT(IN)             :: n       ! This is actually p
    integer(kind = c_int), INTENT(IN OUT)         :: iv(liv)
    real(kind = c_double), INTENT(IN OUT)         :: d(n)
    real(kind = c_double), INTENT(IN OUT)         :: v(lv)
    ! integer(kind = c_int), INTENT(IN)             :: liv_1
    ! integer(kind = c_int), INTENT(IN)             :: lv_1

    real(kind = c_double) :: dr7mdc
    EXTERNAL divset, dr7mdc,dv7cpy,dv7dfl

    ! DIVSET  -- SUPPLIES DEFAULT VALUES TO BOTH IV AND V.
    ! DR7MDC -- RETURNS MACHINE-DEPENDENT CONSTANTS.
    ! DV7CPY  -- COPIES ONE VECTOR TO ANOTHER.
    ! DV7DFL  -- SUPPLIES DEFAULT PARAMETER VALUES TO V ALONE.

    !  ***  LOCAL VARIABLES  ***

    integer(kind = c_int) :: alg1, i, ii, iv1, j, k, &
    l, m, miv1, miv2, ndfalt, parsv1, pu

    integer(kind = c_int) :: ijmp, jlim(4), miniv(4), ndflt(4)
    integer(kind = c_int) :: kplus1
    CHARACTER (LEN=1)     :: varnm(2), sh(2)
    CHARACTER (LEN=4)     :: cngd(3), dflt(3), vn(2,34), which(3)
    real(kind = c_double) :: big, machep, tiny, vk, &
    vm(34), vx(34), zero

    ! DSB NOTE: Temporarily created for debugging
    ! CHARACTER (LEN=132) :: output_string
    CHARACTER (LEN=132) :: OutputString

    !  ***  IV AND V SUBSCRIPTS  ***

    integer(kind = c_int), PARAMETER :: algsav=51
    integer(kind = c_int), PARAMETER :: dinit=38
    integer(kind = c_int), PARAMETER :: dtype=16
    integer(kind = c_int), PARAMETER :: dtype0=54
    integer(kind = c_int), PARAMETER :: epslon=19
    integer(kind = c_int), PARAMETER :: inits=25
    integer(kind = c_int), PARAMETER :: ivneed=3
    integer(kind = c_int), PARAMETER :: lastiv=44
    integer(kind = c_int), PARAMETER :: lastv=45
    integer(kind = c_int), PARAMETER :: lmat=42
    integer(kind = c_int), PARAMETER :: nextiv=46
    integer(kind = c_int), PARAMETER :: nextv=47
    integer(kind = c_int), PARAMETER :: nvdflt=50
    integer(kind = c_int), PARAMETER :: oldn=38
    integer(kind = c_int), PARAMETER :: parprt=20
    integer(kind = c_int), PARAMETER :: parsav=49
    integer(kind = c_int), PARAMETER :: perm=58
    integer(kind = c_int), PARAMETER :: prunit=21
    integer(kind = c_int), PARAMETER :: vneed=4

    SAVE big, machep, tiny

    DATA big/0.d+0/, machep/-1.d+0/, tiny/1.d+0/, zero/0.d+0/
    DATA vn(1,1),vn(2,1)/'EPSL','ON..'/
    DATA vn(1,2),vn(2,2)/'PHMN','FC..'/
    DATA vn(1,3),vn(2,3)/'PHMX','FC..'/
    DATA vn(1,4),vn(2,4)/'DECF','AC..'/
    DATA vn(1,5),vn(2,5)/'INCF','AC..'/
    DATA vn(1,6),vn(2,6)/'RDFC','MN..'/
    DATA vn(1,7),vn(2,7)/'RDFC','MX..'/
    DATA vn(1,8),vn(2,8)/'TUNE','R1..'/
    DATA vn(1,9),vn(2,9)/'TUNE','R2..'/
    DATA vn(1,10),vn(2,10)/'TUNE','R3..'/
    DATA vn(1,11),vn(2,11)/'TUNE','R4..'/
    DATA vn(1,12),vn(2,12)/'TUNE','R5..'/
    DATA vn(1,13),vn(2,13)/'AFCT','OL..'/
    DATA vn(1,14),vn(2,14)/'RFCT','OL..'/
    DATA vn(1,15),vn(2,15)/'XCTO','L...'/
    DATA vn(1,16),vn(2,16)/'XFTO','L...'/
    DATA vn(1,17),vn(2,17)/'LMAX','0...'/
    DATA vn(1,18),vn(2,18)/'LMAX','S...'/
    DATA vn(1,19),vn(2,19)/'SCTO','L...'/
    DATA vn(1,20),vn(2,20)/'DINI','T...'/
    DATA vn(1,21),vn(2,21)/'DTIN','IT..'/
    DATA vn(1,22),vn(2,22)/'D0IN','IT..'/
    DATA vn(1,23),vn(2,23)/'DFAC','....'/
    DATA vn(1,24),vn(2,24)/'DLTF','DC..'/
    DATA vn(1,25),vn(2,25)/'DLTF','DJ..'/
    DATA vn(1,26),vn(2,26)/'DELT','A0..'/
    DATA vn(1,27),vn(2,27)/'FUZZ','....'/
    DATA vn(1,28),vn(2,28)/'RLIM','IT..'/
    DATA vn(1,29),vn(2,29)/'COSM','IN..'/
    DATA vn(1,30),vn(2,30)/'HUBE','RC..'/
    DATA vn(1,31),vn(2,31)/'RSPT','OL..'/
    DATA vn(1,32),vn(2,32)/'SIGM','IN..'/
    DATA vn(1,33),vn(2,33)/'ETA0','....'/
    DATA vn(1,34),vn(2,34)/'BIAS','....'/

    DATA vm(1)/1.0D-3/, vm(2)/-0.99D+0/, vm(3)/1.0D-3/, vm(4)/1.0D-2/,  &
    vm(5)/1.2D+0/, vm(6)/1.d-2/, vm(7)/1.2D+0/, vm(8)/0.d+0/,  &
    vm(9)/0.d+0/, vm(10)/1.d-3/, vm(11)/-1.d+0/, vm(13)/0.d+0/,  &
    vm(15)/0.d+0/, vm(16)/0.d+0/, vm(19)/0.d+0/, vm(20)/-10.d+0/,  &
    vm(21)/0.d+0/, vm(22)/0.d+0/, vm(23)/0.d+0/, vm(27)/1.01D+0/,  &
    vm(28)/1.d+10/, vm(30)/0.d+0/, vm(31)/0.d+0/, vm(32)/0.d+0/, &
    vm(34)/0.d+0/
    DATA vx(1)/0.9D+0/, vx(2)/-1.d-3/, vx(3)/1.d+1/, vx(4)/0.8D+0/,  &
    vx(5)/1.d+2/, vx(6)/0.8D+0/, vx(7)/1.d+2/, vx(8)/0.5D+0/,  &
    vx(9)/0.5D+0/, vx(10)/1.d+0/, vx(11)/1.d+0/, vx(14)/0.1D+0/,  &
    vx(15)/1.d+0/, vx(16)/1.d+0/, vx(19)/1.d+0/, vx(23)/1.d+0/,  &
    vx(24)/1.d+0/, vx(25)/1.d+0/, vx(26)/1.d+0/, vx(27)/1.d+10/,  &
    vx(29)/1.d+0/, vx(31)/1.d+0/, vx(32)/1.d+0/, vx(33)/1.d+0/, &
    vx(34)/1.d+0/

    DATA varnm(1)/'P'/, varnm(2)/'P'/, sh(1)/'S'/, sh(2)/'H'/
    DATA cngd(1),cngd(2),cngd(3)/'---C','HANG','ED V'/,  &
    dflt(1),dflt(2),dflt(3)/'NOND','EFAU','LT V'/
    DATA ijmp/33/, jlim(1)/0/, jlim(2)/24/, jlim(3)/0/, jlim(4)/24/,  &
    ndflt(1)/32/, ndflt(2)/25/, ndflt(3)/32/, ndflt(4)/25/
    DATA miniv(1)/82/, miniv(2)/59/, miniv(3)/103/, miniv(4)/103/

    !...............................  BODY  ................................
    pu = 6
    IF (prunit <= liv) pu = iv(prunit)
    IF (algsav > liv) GO TO 20
    IF (alg .EQ. iv(algsav)) GO TO 20
    ! IF (pu /= 0) WRITE(pu,10) alg, iv(algsav)
    ! 10      FORMAT(/40H the first PARAMETER TO divset should be,i3,  &
    !     12H rather than,i3)
    OutputString = ' The first argument to divset should be 1'
    iv(1) = 67
    GO TO 999
    20  CONTINUE
    IF (alg < 1 .OR. alg > 4) GO TO 340
    miv1 = miniv(alg)
    IF (iv(1) == 15) GO TO 360
    alg1 = MOD(alg-1,2) + 1
    IF (iv(1) == 0) CALL divset(alg, iv, liv, lv, v)
    iv1 = iv(1)
    IF (iv1 /= 13 .AND. iv1 /= 12) GO TO 30
    IF (perm <= liv) miv1 = MAX0(miv1, iv(perm) - 1)
    miv2 = 0
    IF (ivneed <= liv) miv2 = miv1 + MAX0(iv(ivneed), 0)
    IF (lastiv <= liv) iv(lastiv) = miv2
    IF (liv < miv1) GO TO 300
    iv(ivneed) = 0
    iv(lastv) = MAX0(iv(vneed), 0) + iv(lmat) - 1
    iv(vneed) = 0
    IF (liv < miv2) GO TO 300
    IF (lv < iv(lastv)) GO TO 320
    30  CONTINUE
    IF (iv1 < 12 .OR. iv1 > 14) GO TO 60
    IF (n >= 1) GO TO 50
    iv(1) = 81
    IF (pu == 0) GO TO 999
    ! WRITE(pu,40) varnm(alg1), n
    ! 40           FORMAT(/8H /// bad,a1,2H =,i5)
    OutputString = ' /// bad '//varnm(alg1)//' value...'
    GO TO 999
    50      IF (iv1 /= 14) iv(nextiv) = iv(perm)
    IF (iv1 /= 14) iv(nextv) = iv(lmat)
    IF (iv1 == 13) GO TO 999
    k = iv(parsav) - epslon
    kplus1 = k+1
    CALL dv7dfl(alg1, lv-k, v(kplus1))
    iv(dtype0) = 2 - alg1
    iv(oldn) = n
    which(1) = dflt(1)
    which(2) = dflt(2)
    which(3) = dflt(3)
    GO TO 110
    60   IF (n == iv(oldn)) GO TO 80
    iv(1) = 17
    IF (pu == 0) GO TO 999
    !WRITE(pu,70) varnm(alg1), iv(oldn), n
    ! 70      FORMAT(/5H /// ,1A1,14H changed from ,i5,4H TO ,i5)
    OutputString = ' /// '//varnm(alg1)//' value changed inappropriately...'
    GO TO 999

    80   IF (iv1 <= 11 .AND. iv1 >= 1) GO TO 100
    iv(1) = 80
    ! IF (pu /= 0) WRITE(pu,90) iv1
    ! 90      FORMAT(/13H ///  iv(1) =,i5,28H should be between 0 AND 14.)
    OutputString = ' /// iv(1) value is bad. It should be between 0 and 14.'
    GO TO 999

    100  which(1) = cngd(1)
    which(2) = cngd(2)
    which(3) = cngd(3)

    110  IF (iv1 == 14) iv1 = 12
    IF (big > tiny) GO TO 120
    tiny = dr7mdc(1)
    machep = dr7mdc(3)
    big = dr7mdc(6)
    vm(12) = machep
    vx(12) = big
    vx(13) = big
    vm(14) = machep
    vm(17) = tiny
    vx(17) = big
    vm(18) = tiny
    vx(18) = big
    vx(20) = big
    vx(21) = big
    vx(22) = big
    vm(24) = machep
    vm(25) = machep
    vm(26) = machep
    vx(28) = dr7mdc(5)
    vm(29) = machep
    vx(30) = big
    vm(33) = machep
    120  m = 0
    i = 1
    j = jlim(alg1)
    k = epslon
    ndfalt = ndflt(alg1)
    DO  l = 1, ndfalt
    vk = v(k)
    IF (vk >= vm(i) .AND. vk <= vx(i)) GO TO 140
    m = k
    ! IF (pu /= 0) WRITE(pu,130) vn(1,i), vn(2,i), k, vk, vm(i), vx(i)
    ! 130          FORMAT(/6H ///  ,2A4,5H.. v(,i2,3H) =,e11.3,7H should,  &
    !    11H be between,e11.3,4H AND,e11.3)
    OutputString = vn(1,i)//vn(2,i)//' has a bad value. Should be between values stored in vm(i) and //vx(i)'
    140     k = k + 1
    i = i + 1
    IF (i == j) i = ijmp
    END DO
    IF (iv(nvdflt) == ndfalt) GO TO 170
    iv(1) = 51
    IF (pu == 0) GO TO 999
    ! WRITE(pu,160) iv(nvdflt), ndfalt
    ! 160     FORMAT(/13H iv(nvdflt) =,i5,13H rather than ,i5)
    OutputString = ' iv(nvdflt) has a bad value. Should be 32 for alg=1.'
    GO TO 999
    170  IF ((iv(dtype) > 0 .OR. v(dinit) > zero) .AND. iv1 == 12) GO TO 200
    DO  i = 1, n
    IF (d(i) > zero) CYCLE
    m = 18
    ! IF (pu /= 0) WRITE(pu,180) i, d(i)
    ! 180     FORMAT(/8H ///  d(,i3,3H) =,e11.3,19H should be positive)
    OutputString = ' /// Bad value in scaling vector d(). It has a non-positive value.'
    END DO
    200  IF (m == 0) GO TO 210
    iv(1) = m
    GO TO 999

    210  IF (pu == 0 .OR. iv(parprt) == 0) GO TO 999
    IF (iv1 /= 12 .OR. iv(inits) == alg1-1) GO TO 230
    m = 1
    ! WRITE(pu,220) sh(alg1), iv(inits)
    ! 220     FORMAT(/22H nondefault values..../5H init,a1,14H..... iv(25) =, i3)
    OutputString = ' Non-default value is being used in iv(inits) = iv(25).'
    230  IF (iv(dtype) == iv(dtype0)) GO TO 250
    ! IF (m == 0) WRITE(pu,260) which
    IF (m==0) OutputString = which(1)//which(2)//which(3)//' being used. '
    m = 1
    ! WRITE(pu,240) iv(dtype)
    ! 240     FORMAT(20H dtype..... iv(16) =,i3)
    OutputString = trim(OutputString)//' See dtype in iv(16).'
    250  i = 1
    j = jlim(alg1)
    k = epslon
    l = iv(parsav)
    ndfalt = ndflt(alg1)
    DO  ii = 1, ndfalt
    IF (v(k) == v(l)) GO TO 280
    ! IF (m == 0) WRITE(pu,260) which
    ! 260          FORMAT(/1H ,3A4,9HALUES..../)
    IF (m==0) OutputString = which(1)//which(2)//which(3)//'alues being used  '
    m = 1
    ! WRITE(pu,270) vn(1,i), vn(2,i), k, v(k)
    ! 270          FORMAT(1X,2A4,5H.. v(,i2,3H) =,e15.7)
    OutputString = trim(OutputString)//' for '//vn(1,i)//vn(2,i)
    280     k = k + 1
    l = l + 1
    i = i + 1
    IF (i == j) i = ijmp
    END DO

    iv(dtype0) = iv(dtype)
    parsv1 = iv(parsav)
    CALL dv7cpy(iv(nvdflt), v(parsv1), v(epslon))
    GO TO 999

    300  iv(1) = 15
    IF (pu == 0) GO TO 999
    ! WRITE(pu,310) liv, miv2
    ! 310  FORMAT(/10H /// liv =,i5,17H must be at least,i5)
    OutputString = ' /// Bad value for liv. It is too small.'
    IF (liv < miv1) GO TO 999
    IF (lv < iv(lastv)) GO TO 320
    GO TO 999

    320  iv(1) = 16
    ! IF (pu /= 0) WRITE(pu,330) lv, iv(lastv)
    ! 330  FORMAT(/9H /// lv =,i5,17H must be at least,i5)
    OutputString = ' /// Bad value for lv. It is too small.'
    GO TO 999

    340  iv(1) = 67
    ! IF (pu /= 0) WRITE(pu,350) alg
    ! 350  FORMAT(/10H /// alg =,i5,21H must be 1 2, 3, OR 4)
    OutputString = ' Bad value for alg. It must be 1, 2, 3, or 4.'
    GO TO 999
    360 CONTINUE
    ! 360  IF (pu /= 0) WRITE(pu,370) liv, miv1
    ! 370  FORMAT(/10H /// liv =,i5,17H must be at least,i5,  &
    !    37H TO compute true MIN. liv AND MIN. lv)
    OutputString = ' /// Bad value for liv. It is too small to compute true MIN liv and min lv.'
    IF (lastiv <= liv) iv(lastiv) = miv1
    IF (lastv <= liv) iv(lastv) = 0

    999  RETURN
    !  ***  LAST LINE OF dparck_m FOLLOWS  ***
    END SUBROUTINE dparck_m

     SUBROUTINE RHO(NEED, F, NOBS, &
	  &  NF, XN, R, RD, RHOI, RHOR, W)
      ! bind(C, name="rho_")
	  ! FIX ME I need to review what the differences are between
	  ! this one and the previous one.
! Version for use with dglfg_pcm.
! Standard BGW subroutine for Maximum Likelihood for Probabilistic
! Choice Models, i.e., rho_sub_i = -log(), and r_sub_i = Prob.
!
! FIX ME We need to deal with any global variables.
!  USE PCM_DATABASE_GLOBAL_VARS
! [See MODULE PCM_DATABASE_GLOBAL_VARS for variable definitions.]
!
    use, intrinsic :: iso_c_binding
    use ieee_arithmetic
    implicit none

    integer(kind = c_int), INTENT (IN)         :: NEED(2)
    integer(kind = c_int), INTENT (IN)         :: NOBS
    integer(kind = c_int), INTENT (IN OUT)     :: NF
    integer(kind = c_int), INTENT (IN OUT)     :: RHOI(*)

    real(kind = c_double), INTENT(OUT)         :: F
    real(kind = c_double), INTENT(IN)          :: XN(*)
    real(kind = c_double), INTENT (IN OUT)     :: RHOR(*)
    real(kind = c_double), INTENT (IN OUT)     :: R(*)
    real(kind = c_double), INTENT (IN OUT)     :: RD(NOBS,*)
    real(kind = c_double), INTENT (IN OUT)     :: W(NOBS)
!
    ! integer(kind = c_int) :: ICP, IOBS, IRW, KS
    integer(kind = c_int) :: ICP, IOBS, IRW
    integer(kind = c_int) :: Ind_Weight
    real(kind = c_double) :: OOR, VT
!
!    DOUBLE PRECISION NEGONE, ZERO
    real(kind = c_double), PARAMETER :: NEGONE=-1.D0
    real(kind = c_double), PARAMETER :: ZERO=0.D0

    ! DSB NOTE: Temporarily added for debugging
    ! CHARACTER (len=132) output_string

! DSB NOTE:  The following has been added to stop a compiler error.
! This version of RHO does not happen to use XN, so the compiler complains
! about it!
  VT = XN(1)
!
! *** BODY ***
!
! *** Note:  It seems CLEAR that I should re-write the
! *** following using f90 features!
!  WEIGHT = UI(14)
    Ind_Weight = 0
!      ! output_string = ""
!      WRITE(output_string,5) "\n", "Arrived in RHO"
!      CALL mexPrintf(! output_string)
!   5  FORMAT(A2,A14)
    IF (RHOI(1) .EQ. 1) THEN
        Ind_Weight = 1
    ENDIF
!      ! output_string = ""
!      WRITE(output_string,4) "\n", RHOI(1), Ind_Weight
!      CALL mexPrintf(! output_string)
!   4  FORMAT(A2,1X,I4,1X,I4)
!      ! output_string = ""
!      WRITE(output_string,4) "\n", NEED(1)
!      CALL mexPrintf(! output_string)
    IF (NEED(1).EQ.1) THEN
       VT = ZERO
       IF (Ind_Weight.EQ.0) THEN
          DO IOBS = 1, NOBS
              VT = VT - DLOG(R(IOBS))
          END DO
       ELSE
          IRW = 2*NOBS
          IRW = 1
          DO IOBS = 1, NOBS
!              IF (IOBS < 6) THEN
!                ! output_string = ""
!                WRITE(output_string,3) "\n", RHOR(IRW), R(IOBS)
!                CALL mexPrintf(! output_string)
!   3            FORMAT(A2,1X,D13.6,1X,D13.6)
!              ENDIF
             VT = VT - RHOR(IRW) * DLOG(R(IOBS))
             IRW = IRW + 1
          END DO
       ENDIF
       F = VT
!      ! output_string = ""
!      WRITE(output_string,2) "\n", F
!   2  FORMAT(A2,1X,D13.6)
!        CALL mexPrintf(! output_string)
	   ! IF (ISNAN(F)) THEN
	   !   PRINT *, NF, ' NaN detected...'
	   !  NF = 0
	!	 RETURN
	   ! ENDIF
	   ! IF ( .NOT. (F .GE.
       ! NOTE: At times in the past it was necessary to
       ! include an internal version of ISNAN_DP. We
       ! have eliminated it since we believe compilers now
       ! have generic, intrinsic ISNAN functions.
	   ! PRINT *, NF, F, ISNAN_DP(F)
       ! July 10, 2023.
       ! CRAN disallows use of ISNAN (which it says is an GNU extension):
       ! "isnan is a GNU extension.  There are standard ways to do this as from F2003,
       ! or you can use if(my_var /= my_var)."
       ! (This latter approach is what we were originally doing with ISNAN_DP(F).
       ! The F2003 standard includes ieee_is_nan(x), which we will use (and
       ! hope this satisfies CRAN requirements).
       ! IF (ISNAN(F)) THEN
       IF (ieee_is_nan(F)) THEN
         NF = 0
         RETURN
       ENDIF
    ELSE
!       KS = 1
!    IF (UI(11).NE.NF) KS = 2
!    IF (UI(10+KS).NE.NF) THEN
!     IF (RHOI(1).NE.NF) KS = 2
!     IF (RHOI(KS).NE.NF) THEN
!        ! IOUNIT = UI(13)
!        WRITE(IOUNIT,*) ' PROBLEM WITH INITIAL POINT...'
!        NF = 0
!        RETURN
!     ENDIF
!    ICP = UI(8)
     ICP = 1
!     IF (KS.EQ.2) ICP = ICP + NOBS
       IF (Ind_Weight.EQ.0) THEN
          DO IOBS = 1, NOBS
!            OOR = NEGONE/RHOR(ICP)
              OOR = NEGONE/R(IOBS)
              R(IOBS) = OOR
              W(IOBS) = R(IOBS) * OOR
              RD(IOBS,1) = W(IOBS)
!            ICP = ICP + 1
          END DO
       ELSE
          ! IRW = 1 + 2*NOBS
          IRW = 1
          ICP = 1
          DO IOBS = 1, NOBS
              ! OOR = NEGONE/RHOR(ICP)
              OOR = NEGONE/R(ICP)
              R(IOBS) = RHOR(IRW) * OOR
              W(IOBS) = R(IOBS) * OOR
              RD(IOBS,1) = W(IOBS)
              ICP = ICP + 1
              IRW = IRW + 1
          END DO
          !DO IOBS = 1, NOBS
          !    OOR = NEGONE/RHOR(IOBS)
          !    R(IOBS) = RHOR(IOBS) * OOR
          !    W(IOBS) = R(IOBS) * OOR
          !    RD(IOBS,1) = W(IOBS)
          !END DO
       ENDIF
    ENDIF
    RETURN
! *** LAST LINE OF PCMRHO FOLLOWS ***
  ! END SUBROUTINE CM_MLE_RHO_F
  END SUBROUTINE RHO
! end module
