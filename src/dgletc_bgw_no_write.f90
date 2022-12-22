    ! dgletc_bgw.f90
    ! Version for BGW R package. August 5, 2022
    ! The major issue is that R packages carefully check Fortran code
    ! for any WRITE or PRINT statements (regardless of whether they call them).
    ! So the entire package must be cleansed of any such statements.
    ! The main issue is output routines:  ditsum and
    ! For now, we create dummy versions of these that do not write anything,
    ! so that their original place in the code is preserved.
    ! In addition, SUBROUTINE dparck is completely removed.
    ! We can do that in this file because this is not called in this version.
    ! [There is a modified version currently in drglg_mod, but we may need to
    ! figure out how to remove that one as well.] 
    SUBROUTINE da7sst(iv, liv, lv, v)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2021-07-13  Time: 16:55:52

    !  ***  ASSESS CANDIDATE STEP (***SOL VERSION 2.3)  ***

    INTEGER, INTENT(IN)                  :: liv
    INTEGER, INTENT(IN)                  :: lv
    INTEGER, INTENT(IN OUT)              :: iv(liv)
    DOUBLE PRECISION, INTENT(IN OUT)     :: v(lv)

    !  ***  PURPOSE  ***

    !        THIS SUBROUTINE IS CALLED BY AN UNCONSTRAINED MINIMIZATION
    !     ROUTINE TO ASSESS THE NEXT CANDIDATE STEP.  IT MAY RECOMMEND ONE
    !     OF SEVERAL COURSES OF ACTION, SUCH AS ACCEPTING THE STEP, RECOM-
    !     PUTING IT USING THE SAME OR A NEW QUADRATIC MODEL, OR HALTING DUE
    !     TO CONVERGENCE OR FALSE CONVERGENCE.  SEE THE RETURN CODE LISTING
    !     BELOW.

    !--------------------------  PARAMETER USAGE  --------------------------

    !  IV (I/O) INTEGER PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION
    !             BELOW OF IV VALUES REFERENCED.
    ! LIV (IN)  LENGTH OF IV ARRAY.
    !  LV (IN)  LENGTH OF V ARRAY.
    !   V (I/O) REAL PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION
    !             BELOW OF V VALUES REFERENCED.

    !  ***  IV VALUES REFERENCED  ***

    !    IV(IRC) (I/O) ON INPUT FOR THE FIRST STEP TRIED IN A NEW ITERATION,
    !             IV(IRC) SHOULD BE SET TO 3 OR 4 (THE VALUE TO WHICH IT IS
    !             SET WHEN STEP IS DEFINITELY TO BE ACCEPTED).  ON INPUT
    !             AFTER STEP HAS BEEN RECOMPUTED, IV(IRC) SHOULD BE
    !             UNCHANGED SINCE THE PREVIOUS RETURN OF DA7SST.
    !                ON OUTPUT, IV(IRC) IS A RETURN CODE HAVING ONE OF THE
    !             FOLLOWING VALUES...
    !                  1 = SWITCH MODELS OR TRY SMALLER STEP.
    !                  2 = SWITCH MODELS OR ACCEPT STEP.
    !                  3 = ACCEPT STEP AND DETERMINE V(RADFAC) BY GRADIENT
    !                       TESTS.
    !                  4 = ACCEPT STEP, V(RADFAC) HAS BEEN DETERMINED.
    !                  5 = RECOMPUTE STEP (USING THE SAME MODEL).
    !                  6 = RECOMPUTE STEP WITH RADIUS = V(LMAXS) BUT DO NOT
    !                       EVAULATE THE OBJECTIVE FUNCTION.
    !                  7 = X-CONVERGENCE (SEE V(XCTOL)).
    !                  8 = RELATIVE FUNCTION CONVERGENCE (SEE V(RFCTOL)).
    !                  9 = BOTH X- AND RELATIVE FUNCTION CONVERGENCE.
    !                 10 = ABSOLUTE FUNCTION CONVERGENCE (SEE V(AFCTOL)).
    !                 11 = SINGULAR CONVERGENCE (SEE V(LMAXS)).
    !                 12 = FALSE CONVERGENCE (SEE V(XFTOL)).
    !                 13 = IV(IRC) WAS OUT OF RANGE ON INPUT.
    !             RETURN CODE I HAS PRECDENCE OVER I+1 FOR I = 9, 10, 11.
    ! IV(MLSTGD) (I/O) SAVED VALUE OF IV(MODEL).
    !  IV(MODEL) (I/O) ON INPUT, IV(MODEL) SHOULD BE AN INTEGER IDENTIFYING
    !             THE CURRENT QUADRATIC MODEL OF THE OBJECTIVE FUNCTION.
    !             IF A PREVIOUS STEP YIELDED A BETTER FUNCTION REDUCTION,
    !             THEN IV(MODEL) WILL BE SET TO IV(MLSTGD) ON OUTPUT.
    ! IV(NFCALL) (IN)  INVOCATION COUNT FOR THE OBJECTIVE FUNCTION.
    ! IV(NFGCAL) (I/O) VALUE OF IV(NFCALL) AT STEP THAT GAVE THE BIGGEST
    !             FUNCTION REDUCTION THIS ITERATION.  IV(NFGCAL) REMAINS
    !             UNCHANGED UNTIL A FUNCTION REDUCTION IS OBTAINED.
    ! IV(RADINC) (I/O) THE NUMBER OF RADIUS INCREASES (OR MINUS THE NUMBER
    !             OF DECREASES) SO FAR THIS ITERATION.
    ! IV(RESTOR) (OUT) SET TO 1 IF V(F) HAS BEEN RESTORED AND X SHOULD BE
    !             RESTORED TO ITS INITIAL VALUE, TO 2 IF X SHOULD BE SAVED,
    !             TO 3 IF X SHOULD BE RESTORED FROM THE SAVED VALUE, AND TO
    !             0 OTHERWISE.
    !  IV(STAGE) (I/O) COUNT OF THE NUMBER OF MODELS TRIED SO FAR IN THE
    !             CURRENT ITERATION.
    ! IV(STGLIM) (IN)  MAXIMUM NUMBER OF MODELS TO CONSIDER.
    ! IV(SWITCH) (OUT) SET TO 0 UNLESS A NEW MODEL IS BEING TRIED AND IT
    !             GIVES A SMALLER FUNCTION VALUE THAN THE PREVIOUS MODEL,
    !             IN WHICH CASE DA7SST SETS IV(SWITCH) = 1.
    ! IV(TOOBIG) (IN)  IS NONZERO IF STEP WAS TOO BIG (E.G. IF IT CAUSED
    !             OVERFLOW).
    !   IV(XIRC) (I/O) VALUE THAT IV(IRC) WOULD HAVE IN THE ABSENCE OF
    !             CONVERGENCE, FALSE CONVERGENCE, AND OVERSIZED STEPS.

    !  ***  V VALUES REFERENCED  ***

    ! V(AFCTOL) (IN)  ABSOLUTE FUNCTION CONVERGENCE TOLERANCE.  IF THE
    !             ABSOLUTE VALUE OF THE CURRENT FUNCTION VALUE V(F) IS LESS
    !             THAN V(AFCTOL) AND DA7SST DOES NOT RETURN WITH
    !             IV(IRC) = 11, THEN DA7SST RETURNS WITH IV(IRC) = 10.
    ! V(DECFAC) (IN)  FACTOR BY WHICH TO DECREASE RADIUS WHEN IV(TOOBIG) IS
    !             NONZERO.
    ! V(DSTNRM) (IN)  THE 2-NORM OF D*STEP.
    ! V(DSTSAV) (I/O) VALUE OF V(DSTNRM) ON SAVED STEP.
    !   V(DST0) (IN)  THE 2-NORM OF D TIMES THE NEWTON STEP (WHEN DEFINED,
    !             I.E., FOR V(NREDUC) .GE. 0).
    !      V(F) (I/O) ON BOTH INPUT AND OUTPUT, V(F) IS THE OBJECTIVE FUNC-
    !             TION VALUE AT X.  IF X IS RESTORED TO A PREVIOUS VALUE,
    !             THEN V(F) IS RESTORED TO THE CORRESPONDING VALUE.
    !   V(FDIF) (OUT) THE FUNCTION REDUCTION V(F0) - V(F) (FOR THE OUTPUT
    !             VALUE OF V(F) IF AN EARLIER STEP GAVE A BIGGER FUNCTION
    !             DECREASE, AND FOR THE INPUT VALUE OF V(F) OTHERWISE).
    ! V(FLSTGD) (I/O) SAVED VALUE OF V(F).
    !     V(F0) (IN)  OBJECTIVE FUNCTION VALUE AT START OF ITERATION.
    ! V(GTSLST) (I/O) VALUE OF V(GTSTEP) ON SAVED STEP.
    ! V(GTSTEP) (IN)  INNER PRODUCT BETWEEN STEP AND GRADIENT.
    ! V(INCFAC) (IN)  MINIMUM FACTOR BY WHICH TO INCREASE RADIUS.
    !  V(LMAXS) (IN)  MAXIMUM REASONABLE STEP SIZE (AND INITIAL STEP BOUND).
    !             IF THE ACTUAL FUNCTION DECREASE IS NO MORE THAN TWICE
    !             WHAT WAS PREDICTED, IF A RETURN WITH IV(IRC) = 7, 8, OR 9
    !             DOES NOT OCCUR, IF V(DSTNRM) .GT. V(LMAXS) OR THE CURRENT
    !             STEP IS A NEWTON STEP, AND IF
    !             V(PREDUC) .LE. V(SCTOL) * ABS(V(F0)), THEN DA7SST RETURNS
    !             WITH IV(IRC) = 11.  IF SO DOING APPEARS WORTHWHILE, THEN
    !            DA7SST REPEATS THIS TEST (DISALLOWING A FULL NEWTON STEP)
    !             WITH V(PREDUC) COMPUTED FOR A STEP OF LENGTH V(LMAXS)
    !             (BY A RETURN WITH IV(IRC) = 6).
    ! V(NREDUC) (I/O)  FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR
    !             NEWTON STEP.  IF DA7SST IS CALLED WITH IV(IRC) = 6, I.E.,
    !             IF V(PREDUC) HAS BEEN COMPUTED WITH RADIUS = V(LMAXS) FOR
    !             USE IN THE SINGULAR CONVERVENCE TEST, THEN V(NREDUC) IS
    !             SET TO -V(PREDUC) BEFORE THE LATTER IS RESTORED.
    ! V(PLSTGD) (I/O) VALUE OF V(PREDUC) ON SAVED STEP.
    ! V(PREDUC) (I/O) FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR
    !             CURRENT STEP.
    ! V(RADFAC) (OUT) FACTOR TO BE USED IN DETERMINING THE NEW RADIUS,
    !             WHICH SHOULD BE V(RADFAC)*DST, WHERE  DST  IS EITHER THE
    !             OUTPUT VALUE OF V(DSTNRM) OR THE 2-NORM OF
    !             DIAG(NEWD)*STEP  FOR THE OUTPUT VALUE OF STEP AND THE
    !             UPDATED VERSION, NEWD, OF THE SCALE VECTOR D.  FOR
    !             IV(IRC) = 3, V(RADFAC) = 1.0 IS RETURNED.
    ! V(RDFCMN) (IN)  MINIMUM VALUE FOR V(RADFAC) IN TERMS OF THE INPUT
    !             VALUE OF V(DSTNRM) -- SUGGESTED VALUE = 0.1.
    ! V(RDFCMX) (IN)  MAXIMUM VALUE FOR V(RADFAC) -- SUGGESTED VALUE = 4.0.
    !  V(RELDX) (IN) SCALED RELATIVE CHANGE IN X CAUSED BY STEP, COMPUTED
    !             (E.G.) BY FUNCTION  DRLDST  AS
    !                 MAX (D(I)*ABS(X(I)-X0(I)), 1 .LE. I .LE. P) /
    !                    MAX (D(I)*(ABS(X(I))+ABS(X0(I))), 1 .LE. I .LE. P).
    ! V(RFCTOL) (IN)  RELATIVE FUNCTION CONVERGENCE TOLERANCE.  IF THE
    !             ACTUAL FUNCTION REDUCTION IS AT MOST TWICE WHAT WAS PRE-
    !             DICTED AND  V(NREDUC) .LE. V(RFCTOL)*ABS(V(F0)),  THEN
    !            DA7SST RETURNS WITH IV(IRC) = 8 OR 9.
    !  V(SCTOL) (IN)  SINGULAR CONVERGENCE TOLERANCE -- SEE V(LMAXS).
    ! V(STPPAR) (IN)  MARQUARDT PARAMETER -- 0 MEANS FULL NEWTON STEP.
    ! V(TUNER1) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION
    !             REDUCTION WAS MUCH LESS THAN EXPECTED.  SUGGESTED
    !             VALUE = 0.1.
    ! V(TUNER2) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION
    !             REDUCTION WAS LARGE ENOUGH TO ACCEPT STEP.  SUGGESTED
    !             VALUE = 10**-4.
    ! V(TUNER3) (IN)  TUNING CONSTANT USED TO DECIDE IF THE RADIUS
    !             SHOULD BE INCREASED.  SUGGESTED VALUE = 0.75.
    !  V(XCTOL) (IN)  X-CONVERGENCE CRITERION.  IF STEP IS A NEWTON STEP
    !             (V(STPPAR) = 0) HAVING V(RELDX) .LE. V(XCTOL) AND GIVING
    !             AT MOST TWICE THE PREDICTED FUNCTION DECREASE, THEN
    !            DA7SST RETURNS IV(IRC) = 7 OR 9.
    !  V(XFTOL) (IN)  FALSE CONVERGENCE TOLERANCE.  IF STEP GAVE NO OR ONLY
    !             A SMALL FUNCTION DECREASE AND V(RELDX) .LE. V(XFTOL),
    !             THEN DA7SST RETURNS WITH IV(IRC) = 12.

    !-------------------------------  NOTES  -------------------------------

    !  ***  APPLICATION AND USAGE RESTRICTIONS  ***

    !        THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR
    !     LEAST-SQUARES) PACKAGE.  IT MAY BE USED IN ANY UNCONSTRAINED
    !     MINIMIZATION SOLVER THAT USES DOGLEG, GOLDFELD-QUANDT-TROTTER,
    !     OR LEVENBERG-MARQUARDT STEPS.

    !  ***  ALGORITHM NOTES  ***

    !        SEE (1) FOR FURTHER DISCUSSION OF THE ASSESSING AND MODEL
    !     SWITCHING STRATEGIES.  WHILE NL2SOL CONSIDERS ONLY TWO MODELS,
    !    DA7SST IS DESIGNED TO HANDLE ANY NUMBER OF MODELS.

    !  ***  USAGE NOTES  ***

    !        ON THE FIRST CALL OF AN ITERATION, ONLY THE I/O VARIABLES
    !     STEP, X, IV(IRC), IV(MODEL), V(F), V(DSTNRM), V(GTSTEP), AND
    !     V(PREDUC) NEED HAVE BEEN INITIALIZED.  BETWEEN CALLS, NO I/O
    !     VALUES EXECPT STEP, X, IV(MODEL), V(F) AND THE STOPPING TOLER-
    !     ANCES SHOULD BE CHANGED.
    !        AFTER A RETURN FOR CONVERGENCE OR FALSE CONVERGENCE, ONE CAN
    !     CHANGE THE STOPPING TOLERANCES AND CALL DA7SST AGAIN, IN WHICH
    !     CASE THE STOPPING TESTS WILL BE REPEATED.

    !  ***  REFERENCES  ***

    !     (1) DENNIS, J.E., JR., GAY, D.M., AND WELSCH, R.E. (1981),
    !        AN ADAPTIVE NONLINEAR LEAST-SQUARES ALGORITHM,
    !        ACM TRANS. MATH. SOFTWARE, VOL. 7, NO. 3.

    !     (2) POWELL, M.J.D. (1970)  A FORTRAN SUBROUTINE FOR SOLVING
    !        SYSTEMS OF NONLINEAR ALGEBRAIC EQUATIONS, IN NUMERICAL
    !        METHODS FOR NONLINEAR ALGEBRAIC EQUATIONS, EDITED BY
    !        P. RABINOWITZ, GORDON AND BREACH, LONDON.

    !  ***  HISTORY  ***

    !        JOHN DENNIS DESIGNED MUCH OF THIS ROUTINE, STARTING WITH
    !     IDEAS IN (2). ROY WELSCH SUGGESTED THE MODEL SWITCHING STRATEGY.
    !        DAVID GAY AND STEPHEN PETERS CAST THIS SUBROUTINE INTO A MORE
    !     PORTABLE FORM (WINTER 1977), AND DAVID GAY CAST IT INTO ITS
    !     PRESENT FORM (FALL 1978), WITH MINOR CHANGES TO THE SINGULAR
    !     CONVERGENCE TEST IN MAY, 1984 (TO DEAL WITH FULL NEWTON STEPS).

    !  ***  GENERAL  ***

    !     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
    !     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
    !     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
    !     MCS-7906671.

    !------------------------  EXTERNAL QUANTITIES  ------------------------

    !  ***  NO EXTERNAL FUNCTIONS AND SUBROUTINES  ***

    !--------------------------  LOCAL VARIABLES  --------------------------

    LOGICAL :: goodx
    INTEGER :: i, nfc
    DOUBLE PRECISION :: emax, emaxs, gts, rfac1, xmax

    !  ***  SUBSCRIPTS FOR IV AND V  ***

    !  ***  DATA INITIALIZATIONS  ***

    DOUBLE PRECISION, PARAMETER :: half=0.5D+0
    DOUBLE PRECISION, PARAMETER :: one=1.d+0
    DOUBLE PRECISION, PARAMETER :: onep2=1.2D+0
    DOUBLE PRECISION, PARAMETER :: two=2.d+0
    DOUBLE PRECISION, PARAMETER :: zero=0.d+0

    INTEGER, PARAMETER :: irc=29
    INTEGER, PARAMETER :: mlstgd=32
    INTEGER, PARAMETER :: model=5
    INTEGER, PARAMETER :: nfcall=6
    INTEGER, PARAMETER :: nfgcal=7
    INTEGER, PARAMETER :: radinc=8
    INTEGER, PARAMETER :: restor=9
    INTEGER, PARAMETER :: stage=10
    INTEGER, PARAMETER :: stglim=11
    INTEGER, PARAMETER :: switch=12
    INTEGER, PARAMETER :: toobig=2
    INTEGER, PARAMETER :: xirc=13
    INTEGER, PARAMETER :: afctol=31
    INTEGER, PARAMETER :: decfac=22
    INTEGER, PARAMETER :: dstnrm=2
    INTEGER, PARAMETER :: dst0=3
    INTEGER, PARAMETER :: dstsav=18
    INTEGER, PARAMETER :: f=10
    INTEGER, PARAMETER :: fdif=11
    INTEGER, PARAMETER :: flstgd=12
    INTEGER, PARAMETER :: f0=13
    INTEGER, PARAMETER :: gtslst=14
    INTEGER, PARAMETER :: gtstep=4
    INTEGER, PARAMETER :: incfac=23
    INTEGER, PARAMETER :: lmaxs=36
    INTEGER, PARAMETER :: nreduc=6
    INTEGER, PARAMETER :: plstgd=15
    INTEGER, PARAMETER :: preduc=7
    INTEGER, PARAMETER :: radfac=16
    INTEGER, PARAMETER :: rdfcmn=24
    INTEGER, PARAMETER :: rdfcmx=25
    INTEGER, PARAMETER :: reldx=17
    INTEGER, PARAMETER :: rfctol=32
    INTEGER, PARAMETER :: sctol=37
    INTEGER, PARAMETER :: stppar=5
    INTEGER, PARAMETER :: tuner1=26
    INTEGER, PARAMETER :: tuner2=27
    INTEGER, PARAMETER :: tuner3=28
    INTEGER, PARAMETER :: xctol=33
    INTEGER, PARAMETER :: xftol=34

    !+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++

    nfc = iv(nfcall)
    iv(switch) = 0
    iv(restor) = 0
    rfac1 = one
    goodx = .true.
    i = iv(irc)
    IF (i >= 1 .AND. i <= 12) THEN
      SELECT CASE (i)
       CASE (1)
        GO TO 20
       CASE (2)
        GO TO 30
       CASE (3)
        GO TO 10
       CASE (4)
        GO TO 10
       CASE (5)
        GO TO 40
       CASE (6)
        GO TO 280
       CASE (7)
        GO TO 220
       CASE (8)
        GO TO 220
       CASE (9)
        GO TO 220
       CASE (10)
        GO TO 220
       CASE (11)
        GO TO 220
       CASE (12)
        GO TO 170
      END SELECT
    END IF
    iv(irc) = 13
    GO TO 999

    !  ***  INITIALIZE FOR NEW ITERATION  ***

    10   iv(stage) = 1
    iv(radinc) = 0
    v(flstgd) = v(f0)
    IF (iv(toobig) == 0) GO TO 110
    iv(stage) = -1
    iv(xirc) = i
    GO TO 60

    !  ***  STEP WAS RECOMPUTED WITH NEW MODEL OR SMALLER RADIUS  ***
    !  ***  FIRST DECIDE WHICH  ***

    20   IF (iv(model) /= iv(mlstgd)) GO TO 30
    !        ***  OLD MODEL RETAINED, SMALLER RADIUS TRIED  ***
    !        ***  DO NOT CONSIDER ANY MORE NEW MODELS THIS ITERATION  ***
    iv(stage) = iv(stglim)
    iv(radinc) = -1
    GO TO 110

    !  ***  A NEW MODEL IS BEING TRIED.  DECIDE WHETHER TO KEEP IT.  ***

    30   iv(stage) = iv(stage) + 1

    !     ***  NOW WE ADD THE POSSIBILTIY THAT STEP WAS RECOMPUTED WITH  ***
    !     ***  THE SAME MODEL, PERHAPS BECAUSE OF AN OVERSIZED STEP.     ***

    40   IF (iv(stage) > 0) GO TO 50

    !        ***  STEP WAS RECOMPUTED BECAUSE IT WAS TOO BIG.  ***

    IF (iv(toobig) /= 0) GO TO 60

    !        ***  RESTORE IV(STAGE) AND PICK UP WHERE WE LEFT OFF.  ***

    iv(stage) = -iv(stage)
    i = iv(xirc)
    SELECT CASE ( i )
      CASE (    1)
        GO TO 20
      CASE (    2)
        GO TO  30
      CASE (    3)
        GO TO  110
      CASE (    4)
        GO TO  110
      CASE (    5)
        GO TO  70
    END SELECT

    50   IF (iv(toobig) == 0) GO TO 70

    !  ***  HANDLE OVERSIZE STEP  ***

    IF (iv(radinc) > 0) GO TO 80
    iv(stage) = -iv(stage)
    iv(xirc) = iv(irc)

    60      v(radfac) = v(decfac)
    iv(radinc) = iv(radinc) - 1
    iv(irc) = 5
    iv(restor) = 1
    GO TO 999

    70   IF (v(f) < v(flstgd)) GO TO 110

    !     *** THE NEW STEP IS A LOSER.  RESTORE OLD MODEL.  ***

    IF (iv(model) == iv(mlstgd)) GO TO 80
    iv(model) = iv(mlstgd)
    iv(switch) = 1

    !     ***  RESTORE STEP, ETC. ONLY IF A PREVIOUS STEP DECREASED V(F).

    80   IF (v(flstgd) >= v(f0)) GO TO 110
    iv(restor) = 1
    v(f) = v(flstgd)
    v(preduc) = v(plstgd)
    v(gtstep) = v(gtslst)
    IF (iv(switch) == 0) rfac1 = v(dstnrm) / v(dstsav)
    v(dstnrm) = v(dstsav)
    nfc = iv(nfgcal)
    goodx = .false.

    110  v(fdif) = v(f0) - v(f)
    IF (v(fdif) > v(tuner2) * v(preduc)) GO TO 140
    IF (iv(radinc) > 0) GO TO 140

    !        ***  NO (OR ONLY A TRIVIAL) FUNCTION DECREASE
    !        ***  -- SO TRY NEW MODEL OR SMALLER RADIUS

    IF (v(f) < v(f0)) GO TO 120
    iv(mlstgd) = iv(model)
    v(flstgd) = v(f)
    v(f) = v(f0)
    iv(restor) = 1
    GO TO 130
    120     iv(nfgcal) = nfc
    130     iv(irc) = 1
    IF (iv(stage) < iv(stglim)) GO TO 160
    iv(irc) = 5
    iv(radinc) = iv(radinc) - 1
    GO TO 160

    !  ***  NONTRIVIAL FUNCTION DECREASE ACHIEVED  ***

    140  iv(nfgcal) = nfc
    rfac1 = one
    v(dstsav) = v(dstnrm)
    IF (v(fdif) > v(preduc)*v(tuner1)) GO TO 190

    !  ***  DECREASE WAS MUCH LESS THAN PREDICTED -- EITHER CHANGE MODELS
    !  ***  OR ACCEPT STEP WITH DECREASED RADIUS.

    IF (iv(stage) >= iv(stglim)) GO TO 150
    !        ***  CONSIDER SWITCHING MODELS  ***
    iv(irc) = 2
    GO TO 160

    !     ***  ACCEPT STEP WITH DECREASED RADIUS  ***

    150  iv(irc) = 4

    !  ***  SET V(RADFAC) TO FLETCHER*S DECREASE FACTOR  ***

    160  iv(xirc) = iv(irc)
    emax = v(gtstep) + v(fdif)
    v(radfac) = half * rfac1
    IF (emax < v(gtstep)) v(radfac) = rfac1 *   MAX(v(rdfcmn),  &
        half * v(gtstep)/emax)

    !  ***  DO FALSE CONVERGENCE TEST  ***

    170  IF (v(reldx) <= v(xftol)) GO TO 180
    iv(irc) = iv(xirc)
    IF (v(f) < v(f0)) GO TO 200
    GO TO 230

    180  iv(irc) = 12
    GO TO 240

    !  ***  HANDLE GOOD FUNCTION DECREASE  ***

    190  IF (v(fdif) < (-v(tuner3) * v(gtstep))) GO TO 210

    !     ***  INCREASING RADIUS LOOKS WORTHWHILE.  SEE IF WE JUST
    !     ***  RECOMPUTED STEP WITH A DECREASED RADIUS OR RESTORED STEP
    !     ***  AFTER RECOMPUTING IT WITH A LARGER RADIUS.

    IF (iv(radinc) < 0) GO TO 210
    IF (iv(restor) == 1) GO TO 210

    !        ***  WE DID NOT.  TRY A LONGER STEP UNLESS THIS WAS A NEWTON
    !        ***  STEP.

    v(radfac) = v(rdfcmx)
    gts = v(gtstep)
    IF (v(fdif) < (half/v(radfac) - one) * gts)  &
        v(radfac) =   MAX(v(incfac), half*gts/(gts + v(fdif)))
    iv(irc) = 4
    IF (v(stppar) == zero) GO TO 230
    IF (v(dst0) >= zero .AND. (v(dst0) < two*v(dstnrm)  &
        .OR. v(nreduc) < onep2*v(fdif)))  GO TO 230
    !             ***  STEP WAS NOT A NEWTON STEP.  RECOMPUTE IT WITH
    !             ***  A LARGER RADIUS.
    iv(irc) = 5
    iv(radinc) = iv(radinc) + 1

    !  ***  SAVE VALUES CORRESPONDING TO GOOD STEP  ***

    200  v(flstgd) = v(f)
    iv(mlstgd) = iv(model)
    IF (iv(restor) /= 1) iv(restor) = 2
    v(dstsav) = v(dstnrm)
    iv(nfgcal) = nfc
    v(plstgd) = v(preduc)
    v(gtslst) = v(gtstep)
    GO TO 230

    !  ***  ACCEPT STEP WITH RADIUS UNCHANGED  ***

    210  v(radfac) = one
    iv(irc) = 3
    GO TO 230

    !  ***  COME HERE FOR A RESTART AFTER CONVERGENCE  ***

    220  iv(irc) = iv(xirc)
    IF (v(dstsav) >= zero) GO TO 240
    iv(irc) = 12
    GO TO 240

    !  ***  PERFORM CONVERGENCE TESTS  ***

    230  iv(xirc) = iv(irc)
    240  IF (iv(restor) == 1 .AND. v(flstgd) < v(f0)) THEN
            iv(restor) = 3
         END IF
    IF ( ABS(v(f)) < v(afctol)) iv(irc) = 10
    IF (half * v(fdif) > v(preduc)) GO TO 999
    emax = v(rfctol) *  ABS(v(f0))
    emaxs = v(sctol) *  ABS(v(f0))
    IF (v(preduc) <= emaxs .AND. ( v(dstnrm) > v(lmaxs) .OR.  &
        v(stppar) == zero)) THEN
           iv(irc) = 11
    END IF
    IF (v(dst0) < zero) GO TO 250
    i = 0
    IF ((v(nreduc) > zero .AND. v(nreduc) <= emax) .OR.  &
        (v(nreduc) == zero .AND. v(preduc) == zero)) THEN
       i = 2
    END IF
    IF (v(stppar) == zero .AND. v(reldx) <= v(xctol)  &
        .AND. goodx) THEN
      i = i + 1
    END IF
    IF (i > 0) iv(irc) = i + 6

    !  ***  CONSIDER RECOMPUTING STEP OF LENGTH V(LMAXS) FOR SINGULAR
    !  ***  CONVERGENCE TEST.

    250  IF (iv(irc) > 5 .AND. iv(irc) /= 12) GO TO 999
    IF (v(stppar) == zero) GO TO 999
    IF (v(dstnrm) > v(lmaxs)) GO TO 260
    IF (v(preduc) >= emaxs) GO TO 999
    IF (v(dst0) <= zero) GO TO 270
    IF (half * v(dst0) <= v(lmaxs)) GO TO 999
    GO TO 270
    260  IF (half * v(dstnrm) <= v(lmaxs)) GO TO 999
    xmax = v(lmaxs) / v(dstnrm)
    IF (xmax * (two - xmax) * v(preduc) >= emaxs) GO TO 999
    270  IF (v(nreduc) < zero) GO TO 290

    !  ***  RECOMPUTE V(PREDUC) FOR USE IN SINGULAR CONVERGENCE TEST  ***

    v(gtslst) = v(gtstep)
    v(dstsav) = v(dstnrm)
    IF (iv(irc) == 12) v(dstsav) = -v(dstsav)
    v(plstgd) = v(preduc)
    i = iv(restor)
    iv(restor) = 2
    IF (i == 3) iv(restor) = 0
    iv(irc) = 6
    GO TO 999

    !  ***  PERFORM SINGULAR CONVERGENCE TEST WITH RECOMPUTED V(PREDUC)  ***

    280  v(gtstep) = v(gtslst)
    v(dstnrm) =  ABS(v(dstsav))
    iv(irc) = iv(xirc)
    IF (v(dstsav) <= zero) iv(irc) = 12
    v(nreduc) = -v(preduc)
    v(preduc) = v(plstgd)
    iv(restor) = 3
    290  IF (-v(nreduc) <= v(sctol) *  ABS(v(f0))) iv(irc) = 11

    999  RETURN

    !  ***  LAST LINE OF DA7SST FOLLOWS  ***
    END SUBROUTINE da7sst

    DOUBLE PRECISION FUNCTION dd7tpr(p, x, y)

    !  ***  RETURN THE INNER PRODUCT OF THE P-VECTORS X AND Y.  ***

    INTEGER, INTENT(IN)                      :: p
    DOUBLE PRECISION, INTENT(IN)             :: x(p)
    DOUBLE PRECISION, INTENT(IN)             :: y(p)

    INTEGER :: i
    DOUBLE PRECISION :: dr7mdc
    EXTERNAL dr7mdc
    !       DSB FIX ME?  Note these comments, and the commented lines below...
    !  ***  ACTIVATE THE *'ED COMMENT LINES BELOW IF UNDERFLOW IS A PROBLEM.
    !  ***  DR7MDC(2) RETURNS A MACHINE-DEPENDENT CONSTANT, SQTETA, WHICH
    !  ***  IS SLIGHTLY LARGER THAN THE SMALLEST POSITIVE NUMBER THAT
    !  ***  CAN BE SQUARED WITHOUT UNDERFLOWING.

    DOUBLE PRECISION, PARAMETER :: one=1.d+0
    DOUBLE PRECISION, PARAMETER :: zero=0.d+0
    !     DOUBLE PRECISION SQTETA, T
    !      DATA SQTETA/0.D+0/

    dd7tpr = zero
    !      IF (P .LE. 0) GO TO 999
    !      IF (SQTETA .EQ. ZERO) SQTETA = DR7MDC(2)
    DO  i = 1, p
    !         T = DMAX1(DABS(X(I)), DABS(Y(I)))
    !         IF (T .GT. ONE) GO TO 10
    !         IF (T .LT. SQTETA) GO TO 20
    !         T = (X(I)/SQTETA)*Y(I)
    !         IF (DABS(T) .LT. SQTETA) GO TO 20
    !  10      dd7tpr = dd7tpr + x(i)*y(i)
              dd7tpr = dd7tpr + x(i)*y(i)
    END DO

    ! 999  RETURN
    RETURN
    !  ***  LAST LINE OF DD7TPR FOLLOWS  ***
    END FUNCTION dd7tpr

    SUBROUTINE dd7up5(d, iv, liv, lv, p, ps, v)

    !  ***  UPDATE SCALE VECTOR D FOR DG7LIT  ***

    !  ***  PARAMETER DECLARATIONS  ***

    INTEGER, INTENT(IN)                  :: liv
    INTEGER, INTENT(IN)                  :: lv
    INTEGER, INTENT(IN)                  :: p
    INTEGER, INTENT(IN)                  :: ps
    INTEGER, INTENT(IN)                  :: iv(liv)
    DOUBLE PRECISION, INTENT(IN OUT)     :: v(lv)
    DOUBLE PRECISION, INTENT(IN OUT)     :: d(p)

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: d0, hii, i, jtoli, jtol0, r1i, s1
    DOUBLE PRECISION :: t, vdfac

    !     ***  CONSTANTS  ***

    !     ***  EXTERNAL FUNCTIONS  ***

    EXTERNAL dd7tpr
    DOUBLE PRECISION :: dd7tpr

    !  ***  SUBSCRIPTS FOR IV AND V  ***

    INTEGER, PARAMETER :: dfac=41
    INTEGER, PARAMETER :: dtype=16
    INTEGER, PARAMETER :: hc=71
    INTEGER, PARAMETER :: jtol=59
    INTEGER, PARAMETER :: niter=31
    INTEGER, PARAMETER :: rmat=78
    INTEGER, PARAMETER :: s=62

    DOUBLE PRECISION, PARAMETER :: zero=0.d+0

    !  ***  BODY  ***

    IF (iv(dtype) /= 1 .AND. iv(niter) > 0) GO TO 999
    r1i = iv(rmat)
    hii = iv(hc) - 1
    vdfac = v(dfac)
    jtol0 = iv(jtol) - 1
    d0 = jtol0 + p
    s1 = iv(s) - 1
    DO  i = 1, p
      IF (r1i <= 0) GO TO 10
      t = dd7tpr(i, v(r1i), v(r1i))
      r1i = r1i + i
      GO TO 20
      10      hii = hii + i
      t =  ABS(v(hii))
      20      s1 = s1 + i
      IF (i <= ps) t = t +   MAX(v(s1), zero)
      t =  SQRT(t)
      jtoli = jtol0 + i
      d0 = d0 + 1
      IF (t < v(jtoli)) t =   MAX(v(d0), v(jtoli))
      d(i) =   MAX(vdfac*d(i), t)
    END DO

    999  RETURN
    !  ***  LAST LINE OF DD7UP5 FOLLOWS  ***
    END SUBROUTINE dd7up5

    SUBROUTINE dg7qts(d, dig, dihdi, ka, l, p, step, v, w)

    !  *** COMPUTE GOLDFELD-QUANDT-TROTTER STEP BY MORE-HEBDEN TECHNIQUE ***
    !  ***  (NL2SOL VERSION 2.2), MODIFIED A LA MORE AND SORENSEN  ***

    !  ***  PARAMETER DECLARATIONS  ***

    INTEGER, INTENT(IN)                      :: p
    DOUBLE PRECISION, INTENT(IN)             :: d(p)
    DOUBLE PRECISION, INTENT(IN)             :: dig(p)
    DOUBLE PRECISION, INTENT(IN OUT)         :: dihdi(1)
    INTEGER, INTENT(IN OUT)                  :: ka
    DOUBLE PRECISION, INTENT(IN OUT)         :: l(1)
    DOUBLE PRECISION, INTENT(IN OUT)         :: step(p)
    DOUBLE PRECISION, INTENT(IN OUT)         :: v(21)
    DOUBLE PRECISION, INTENT(IN OUT)         :: w(*)


    !     DIMENSION DIHDI(P*(P+1)/2), L(P*(P+1)/2), W(4*P+7)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ***  PURPOSE  ***

    !        GIVEN THE (COMPACTLY STORED) LOWER TRIANGLE OF A SCALED
    !     HESSIAN (APPROXIMATION) AND A NONZERO SCALED GRADIENT VECTOR,
    !     THIS SUBROUTINE COMPUTES A GOLDFELD-QUANDT-TROTTER STEP OF
    !     APPROXIMATE LENGTH V(RADIUS) BY THE MORE-HEBDEN TECHNIQUE.  IN
    !     OTHER WORDS, STEP IS COMPUTED TO (APPROXIMATELY) MINIMIZE
    !     PSI(STEP) = (G**T)*STEP + 0.5*(STEP**T)*H*STEP  SUCH THAT THE
    !     2-NORM OF D*STEP IS AT MOST (APPROXIMATELY) V(RADIUS), WHERE
    !     G  IS THE GRADIENT,  H  IS THE HESSIAN, AND  D  IS A DIAGONAL
    !     SCALE MATRIX WHOSE DIAGONAL IS STORED IN THE PARAMETER D.
    !     (DG7QTS ASSUMES  DIG = D**-1 * G  AND  DIHDI = D**-1 * H * D**-1.)

    !  ***  PARAMETER DESCRIPTION  ***

    !     D (IN)  = THE SCALE VECTOR, I.E. THE DIAGONAL OF THE SCALE
    !              MATRIX  D  MENTIONED ABOVE UNDER PURPOSE.
    !   DIG (IN)  = THE SCALED GRADIENT VECTOR, D**-1 * G.  IF G = 0, THEN
    !              STEP = 0  AND  V(STPPAR) = 0  ARE RETURNED.
    ! DIHDI (IN)  = LOWER TRIANGLE OF THE SCALED HESSIAN (APPROXIMATION),
    !              I.E., D**-1 * H * D**-1, STORED COMPACTLY BY ROWS., I.E.,
    !              IN THE ORDER (1,1), (2,1), (2,2), (3,1), (3,2), ETC.
    !    KA (I/O) = THE NUMBER OF HEBDEN ITERATIONS (SO FAR) TAKEN TO DETER-
    !              MINE STEP.  KA .LT. 0 ON INPUT MEANS THIS IS THE FIRST
    !              ATTEMPT TO DETERMINE STEP (FOR THE PRESENT DIG AND DIHDI)
    !              -- KA IS INITIALIZED TO 0 IN THIS CASE.  OUTPUT WITH
    !              KA = 0  (OR V(STPPAR) = 0)  MEANS  STEP = -(H**-1)*G.
    !     L (I/O) = WORKSPACE OF LENGTH P*(P+1)/2 FOR CHOLESKY FACTORS.
    !     P (IN)  = NUMBER OF PARAMETERS -- THE HESSIAN IS A  P X P  MATRIX.
    !  STEP (I/O) = THE STEP COMPUTED.
    !     V (I/O) CONTAINS VARIOUS CONSTANTS AND VARIABLES DESCRIBED BELOW.
    !     W (I/O) = WORKSPACE OF LENGTH 4*P + 6.

    !  ***  ENTRIES IN V  ***

    ! V(DGNORM) (I/O) = 2-NORM OF (D**-1)*G.
    ! V(DSTNRM) (OUTPUT) = 2-NORM OF D*STEP.
    ! V(DST0)   (I/O) = 2-NORM OF D*(H**-1)*G (FOR POS. DEF. H ONLY), OR
    !             OVERESTIMATE OF SMALLEST EIGENVALUE OF (D**-1)*H*(D**-1).
    ! V(EPSLON) (IN)  = MAX. REL. ERROR ALLOWED FOR PSI(STEP).  FOR THE
    !             STEP RETURNED, PSI(STEP) WILL EXCEED ITS OPTIMAL VALUE
    !             BY LESS THAN -V(EPSLON)*PSI(STEP).  SUGGESTED VALUE = 0.1.
    ! V(GTSTEP) (OUT) = INNER PRODUCT BETWEEN G AND STEP.
    ! V(NREDUC) (OUT) = PSI(-(H**-1)*G) = PSI(NEWTON STEP)  (FOR POS. DEF.
    !             H ONLY -- V(NREDUC) IS SET TO ZERO OTHERWISE).
    ! V(PHMNFC) (IN)  = TOL. (TOGETHER WITH V(PHMXFC)) FOR ACCEPTING STEP
    !             (MORE*S SIGMA).  THE ERROR V(DSTNRM) - V(RADIUS) MUST LIE
    !             BETWEEN V(PHMNFC)*V(RADIUS) AND V(PHMXFC)*V(RADIUS).
    ! V(PHMXFC) (IN)  (SEE V(PHMNFC).)
    !             SUGGESTED VALUES -- V(PHMNFC) = -0.25, V(PHMXFC) = 0.5.
    ! V(PREDUC) (OUT) = PSI(STEP) = PREDICTED OBJ. FUNC. REDUCTION FOR STEP.
    ! V(RADIUS) (IN)  = RADIUS OF CURRENT (SCALED) TRUST REGION.
    ! V(RAD0)   (I/O) = VALUE OF V(RADIUS) FROM PREVIOUS CALL.
    ! V(STPPAR) (I/O) IS NORMALLY THE MARQUARDT PARAMETER, I.E. THE ALPHA
    !             DESCRIBED BELOW UNDER ALGORITHM NOTES.  IF H + ALPHA*D**2
    !             (SEE ALGORITHM NOTES) IS (NEARLY) SINGULAR, HOWEVER,
    !             THEN V(STPPAR) = -ALPHA.

    !  ***  USAGE NOTES  ***

    !     IF IT IS DESIRED TO RECOMPUTE STEP USING A DIFFERENT VALUE OF
    !     V(RADIUS), THEN THIS ROUTINE MAY BE RESTARTED BY CALLING IT
    !     WITH ALL PARAMETERS UNCHANGED EXCEPT V(RADIUS).  (THIS EXPLAINS
    !     WHY STEP AND W ARE LISTED AS I/O).  ON AN INITIAL CALL (ONE WITH
    !     KA .LT. 0), STEP AND W NEED NOT BE INITIALIZED AND ONLY COMPO-
    !     NENTS V(EPSLON), V(STPPAR), V(PHMNFC), V(PHMXFC), V(RADIUS), AND
    !     V(RAD0) OF V MUST BE INITIALIZED.

    !  ***  ALGORITHM NOTES  ***

    !        THE DESIRED G-Q-T STEP (REF. 2, 3, 4, 6) SATISFIES
    !     (H + ALPHA*D**2)*STEP = -G  FOR SOME NONNEGATIVE ALPHA SUCH THAT
    !     H + ALPHA*D**2 IS POSITIVE SEMIDEFINITE.  ALPHA AND STEP ARE
    !     COMPUTED BY A SCHEME ANALOGOUS TO THE ONE DESCRIBED IN REF. 5.
    !     ESTIMATES OF THE SMALLEST AND LARGEST EIGENVALUES OF THE HESSIAN
    !     ARE OBTAINED FROM THE GERSCHGORIN CIRCLE THEOREM ENHANCED BY A
    !     SIMPLE FORM OF THE SCALING DESCRIBED IN REF. 7.  CASES IN WHICH
    !     H + ALPHA*D**2 IS NEARLY (OR EXACTLY) SINGULAR ARE HANDLED BY
    !     THE TECHNIQUE DISCUSSED IN REF. 2.  IN THESE CASES, A STEP OF
    !     (EXACT) LENGTH V(RADIUS) IS RETURNED FOR WHICH PSI(STEP) EXCEEDS
    !     ITS OPTIMAL VALUE BY LESS THAN -V(EPSLON)*PSI(STEP).  THE TEST
    !     SUGGESTED IN REF. 6 FOR DETECTING THE SPECIAL CASE IS PERFORMED
    !     ONCE TWO MATRIX FACTORIZATIONS HAVE BEEN DONE -- DOING SO SOONER
    !     SEEMS TO DEGRADE THE PERFORMANCE OF OPTIMIZATION ROUTINES THAT
    !     CALL THIS ROUTINE.

    !  ***  FUNCTIONS AND SUBROUTINES CALLED  ***

    ! DD7TPR - RETURNS INNER PRODUCT OF TWO VECTORS.
    ! DL7ITV - APPLIES INVERSE-TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX.
    ! DL7IVM - APPLIES INVERSE OF COMPACT LOWER TRIANG. MATRIX.
    ! DL7SRT  - FINDS CHOLESKY FACTOR (OF COMPACTLY STORED LOWER TRIANG.).
    ! DL7SVN - RETURNS APPROX. TO MIN. SING. VALUE OF LOWER TRIANG. MATRIX.
    ! DR7MDC - RETURNS MACHINE-DEPENDENT CONSTANTS.
    ! DV2NRM - RETURNS 2-NORM OF A VECTOR.

    !  ***  REFERENCES  ***

    ! 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), AN ADAPTIVE
    !             NONLINEAR LEAST-SQUARES ALGORITHM, ACM TRANS. MATH.
    !             SOFTWARE, VOL. 7, NO. 3.
    ! 2.  GAY, D.M. (1981), COMPUTING OPTIMAL LOCALLY CONSTRAINED STEPS,
    !             SIAM J. SCI. STATIST. COMPUTING, VOL. 2, NO. 2, PP.
    !             186-197.
    ! 3.  GOLDFELD, S.M., QUANDT, R.E., AND TROTTER, H.F. (1966),
    !             MAXIMIZATION BY QUADRATIC HILL-CLIMBING, ECONOMETRICA 34,
    !             PP. 541-551.
    ! 4.  HEBDEN, M.D. (1973), AN ALGORITHM FOR MINIMIZATION USING EXACT
    !             SECOND DERIVATIVES, REPORT T.P. 515, THEORETICAL PHYSICS
    !             DIV., A.E.R.E. HARWELL, OXON., ENGLAND.
    ! 5.  MORE, J.J. (1978), THE LEVENBERG-MARQUARDT ALGORITHM, IMPLEMEN-
    !             TATION AND THEORY, PP.105-116 OF SPRINGER LECTURE NOTES
    !             IN MATHEMATICS NO. 630, EDITED BY G.A. WATSON, SPRINGER-
    !             VERLAG, BERLIN AND NEW YORK.
    ! 6.  MORE, J.J., AND SORENSEN, D.C. (1981), COMPUTING A TRUST REGION
    !             STEP, TECHNICAL REPORT ANL-81-83, ARGONNE NATIONAL LAB.
    ! 7.  VARGA, R.S. (1965), MINIMAL GERSCHGORIN SETS, PACIFIC J. MATH. 15,
    !             PP. 719-729.

    !  ***  GENERAL  ***

    !     CODED BY DAVID M. GAY.
    !     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
    !     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
    !     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
    !     MCS-7906671.

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ***  LOCAL VARIABLES  ***

    LOGICAL :: restrt
    INTEGER :: dggdmx, diag, diag0, dstsav, emax, emin, i, im1, inc, irc,  &
        j, k, kalim, kamin, k1, lk0, phipin, q, q0, uk0, x
    DOUBLE PRECISION :: alphak, aki, akk, delta, dst, eps, gtsta, lk,  &
        oldphi, phi, phimax, phimin, psifac, rad, radsq,  &
        root, si, sk, sw, t, twopsi, t1, t2, uk, wi

    !     ***  CONSTANTS  ***
    DOUBLE PRECISION :: big, dgxfac

    !  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***

    DOUBLE PRECISION :: dd7tpr, dl7svn, dr7mdc, dv2nrm
    EXTERNAL dd7tpr, dl7itv, dl7ivm,dl7srt, dl7svn, dr7mdc, dv2nrm

    !  ***  SUBSCRIPTS FOR V  ***

    INTEGER, PARAMETER :: dgnorm=1
    INTEGER, PARAMETER :: dstnrm=2
    INTEGER, PARAMETER :: dst0=3
    INTEGER, PARAMETER :: epslon=19
    INTEGER, PARAMETER :: gtstep=4
    INTEGER, PARAMETER :: nreduc=6
    INTEGER, PARAMETER :: phmnfc=20
    INTEGER, PARAMETER :: phmxfc=21
    INTEGER, PARAMETER :: preduc=7
    INTEGER, PARAMETER :: radius=8
    INTEGER, PARAMETER :: rad0=9
    INTEGER, PARAMETER :: stppar=5

    DOUBLE PRECISION, PARAMETER :: epsfac=50.0D+0
    DOUBLE PRECISION, PARAMETER :: four=4.0D+0
    DOUBLE PRECISION, PARAMETER :: half=0.5D+0
    DOUBLE PRECISION, PARAMETER :: kappa=2.0D+0
    DOUBLE PRECISION, PARAMETER :: negone=-1.0D+0
    DOUBLE PRECISION, PARAMETER :: one=1.0D+0
    DOUBLE PRECISION, PARAMETER :: p001=1.0D-3
    DOUBLE PRECISION, PARAMETER :: six=6.0D+0
    DOUBLE PRECISION, PARAMETER :: three=3.0D+0
    DOUBLE PRECISION, PARAMETER :: two=2.0D+0
    DOUBLE PRECISION, PARAMETER :: zero=0.0D+0
    SAVE dgxfac
    DATA big/0.d+0/, dgxfac/0.d+0/

    !  ***  BODY  ***
    ! DSB NOTE the following five lines were added to address a compiler warning
    alphak = one
    kamin = 3
    phi = one
    gtsta = one
    psifac = one
    IF (big <= zero) big = dr7mdc(6)

    !     ***  STORE LARGEST ABS. ENTRY IN (D**-1)*H*(D**-1) AT W(DGGDMX).
    dggdmx = p + 1
    !     ***  STORE GERSCHGORIN OVER- AND UNDERESTIMATES OF THE LARGEST
    !     ***  AND SMALLEST EIGENVALUES OF (D**-1)*H*(D**-1) AT W(EMAX)
    !     ***  AND W(EMIN) RESPECTIVELY.
    emax = dggdmx + 1
    emin = emax + 1
    !     ***  FOR USE IN RECOMPUTING STEP, THE FINAL VALUES OF LK, UK, DST,
    !     ***  AND THE INVERSE DERIVATIVE OF MORE*S PHI AT 0 (FOR POS. DEF.
    !     ***  H) ARE STORED IN W(LK0), W(UK0), W(DSTSAV), AND W(PHIPIN)
    !     ***  RESPECTIVELY.
    lk0 = emin + 1
    phipin = lk0 + 1
    uk0 = phipin + 1
    dstsav = uk0 + 1
    !     ***  STORE DIAG OF (D**-1)*H*(D**-1) IN W(DIAG),...,W(DIAG0+P).
    diag0 = dstsav
    diag = diag0 + 1
    !     ***  STORE -D*STEP IN W(Q),...,W(Q0+P).
    q0 = diag0 + p
    q = q0 + 1
    !     ***  ALLOCATE STORAGE FOR SCRATCH VECTOR X  ***
    x = q + p
    rad = v(radius)
    radsq = rad**2
    !     ***  PHITOL = MAX. ERROR ALLOWED IN DST = V(DSTNRM) = 2-NORM OF
    !     ***  D*STEP.
    phimax = v(phmxfc) * rad
    phimin = v(phmnfc) * rad
    psifac = big
    t1 = two * v(epslon) / (three * (four * (v(phmnfc) + one) *  &
        (kappa + one)  +  kappa  +  two) * rad)
    IF (t1 < big*  MIN(rad,one)) psifac = t1 / rad
    !     ***  OLDPHI IS USED TO DETECT LIMITS OF NUMERICAL ACCURACY.  IF
    !     ***  WE RECOMPUTE STEP AND IT DOES NOT CHANGE, THEN WE ACCEPT IT.
    oldphi = zero
    eps = v(epslon)
    irc = 0
    restrt = .false.
    kalim = ka + 50

    !  ***  START OR RESTART, DEPENDING ON KA  ***

    IF (ka >= 0) GO TO 290

    !  ***  FRESH START  ***

    k = 0
    uk = negone
    ka = 0
    kalim = 50
    v(dgnorm) = dv2nrm(p, dig)
    v(nreduc) = zero
    v(dst0) = zero
    kamin = 3
    IF (v(dgnorm) == zero) kamin = 0

    !     ***  STORE DIAG(DIHDI) IN W(DIAG0+1),...,W(DIAG0+P)  ***

    j = 0
    DO  i = 1, p
      j = j + i
      k1 = diag0 + i
      w(k1) = dihdi(j)
    END DO

    !     ***  DETERMINE W(DGGDMX), THE LARGEST ELEMENT OF DIHDI  ***

    t1 = zero
    j = p * (p + 1) / 2
    DO  i = 1, j
      t =  ABS(dihdi(i))
      IF (t1 < t) t1 = t
    END DO
    w(dggdmx) = t1

    !  ***  TRY ALPHA = 0  ***

    30   CALL dl7srt(1, p, l, dihdi, irc)
    IF (irc == 0) GO TO 50
    !        ***  INDEF. H -- UNDERESTIMATE SMALLEST EIGENVALUE, USE THIS
    !        ***  ESTIMATE TO INITIALIZE LOWER BOUND LK ON ALPHA.
    j = irc*(irc+1)/2
    t = l(j)
    l(j) = one
    DO  i = 1, irc
      w(i) = zero
    END DO
    w(irc) = one
    CALL dl7itv(irc, w, l, w)
    t1 = dv2nrm(irc, w)
    lk = -t / t1 / t1
    v(dst0) = -lk
    IF (restrt) GO TO 210
    GO TO 70

    !     ***  POSITIVE DEFINITE H -- COMPUTE UNMODIFIED NEWTON STEP.  ***
    50   lk = zero
    t = dl7svn(p, l, w(q), w(q))
    IF (t >= one) GO TO 60
    IF (v(dgnorm) >= t*t*big) GO TO 70
    60   CALL dl7ivm(p, w(q), l, dig)
    gtsta = dd7tpr(p, w(q), w(q))
    v(nreduc) = half * gtsta
    CALL dl7itv(p, w(q), l, w(q))
    dst = dv2nrm(p, w(q))
    v(dst0) = dst
    phi = dst - rad
    IF (phi <= phimax) GO TO 260
    IF (restrt) GO TO 210

    !  ***  PREPARE TO COMPUTE GERSCHGORIN ESTIMATES OF LARGEST (AND
    !  ***  SMALLEST) EIGENVALUES.  ***

    70   k = 0
    DO  i = 1, p
      wi = zero
      IF (i == 1) GO TO 90
      im1 = i - 1
      DO  j = 1, im1
        k = k + 1
        t =  ABS(dihdi(k))
        wi = wi + t
        w(j) = w(j) + t
      END DO
      90      w(i) = wi
      k = k + 1
    END DO

    !  ***  (UNDER-)ESTIMATE SMALLEST EIGENVALUE OF (D**-1)*H*(D**-1)  ***

    k = 1
    t1 = w(diag) - w(1)
    IF (p <= 1) GO TO 120
    DO  i = 2, p
      j = diag0 + i
      t = w(j) - w(i)
      IF (t >= t1) CYCLE
      t1 = t
      k = i
    END DO

    120  sk = w(k)
    j = diag0 + k
    akk = w(j)
    k1 = k*(k-1)/2 + 1
    inc = 1
    t = zero
    DO  i = 1, p
      IF (i == k) GO TO 130
      aki =  ABS(dihdi(k1))
      si = w(i)
      j = diag0 + i
      t1 = half * (akk - w(j) + si - aki)
      t1 = t1 +  SQRT(t1*t1 + sk*aki)
      IF (t < t1) t = t1
      IF (i < k) GO TO 140
      130     inc = i
      140     k1 = k1 + inc
    END DO

    w(emin) = akk - t
    uk = v(dgnorm)/rad - w(emin)
    IF (v(dgnorm) == zero) uk = uk + p001 + p001*uk
    IF (uk <= zero) uk = p001

    !  ***  COMPUTE GERSCHGORIN (OVER-)ESTIMATE OF LARGEST EIGENVALUE  ***

    k = 1
    t1 = w(diag) + w(1)
    IF (p <= 1) GO TO 170
    DO  i = 2, p
      j = diag0 + i
      t = w(j) + w(i)
      IF (t <= t1) CYCLE
      t1 = t
      k = i
    END DO

    170  sk = w(k)
    j = diag0 + k
    akk = w(j)
    k1 = k*(k-1)/2 + 1
    inc = 1
    t = zero
    DO  i = 1, p
      IF (i == k) GO TO 180
      aki =  ABS(dihdi(k1))
      si = w(i)
      j = diag0 + i
      t1 = half * (w(j) + si - aki - akk)
      t1 = t1 +  SQRT(t1*t1 + sk*aki)
      IF (t < t1) t = t1
      IF (i < k) GO TO 190
      180     inc = i
      190     k1 = k1 + inc
    END DO

    w(emax) = akk + t
    lk =   MAX(lk, v(dgnorm)/rad - w(emax))

    !     ***  ALPHAK = CURRENT VALUE OF ALPHA (SEE ALG. NOTES ABOVE).  WE
    !     ***  USE MORE*S SCHEME FOR INITIALIZING IT.
    alphak = ABS(v(stppar)) * v(rad0)/rad
    alphak = MIN(uk,MAX(alphak, lk))

    IF (irc /= 0) GO TO 210

    !  ***  COMPUTE L0 FOR POSITIVE DEFINITE H  ***

    CALL dl7ivm(p, w, l, w(q))
    t = dv2nrm(p, w)
    w(phipin) = rad / t / t
    lk = MAX(lk,phi*w(phipin))

    !  ***  SAFEGUARD ALPHAK AND ADD ALPHAK*I TO (D**-1)*H*(D**-1)  ***

    210  ka = ka + 1
    IF (-v(dst0) >= alphak .OR. alphak < lk .OR. alphak >= uk)  &
        alphak = uk *   MAX(p001,  SQRT(lk/uk))
    IF (alphak <= zero) alphak = half * uk
    IF (alphak <= zero) alphak = uk
    k = 0
    DO  i = 1, p
      k = k + i
      j = diag0 + i
      dihdi(k) = w(j) + alphak
    END DO

    !  ***  TRY COMPUTING CHOLESKY DECOMPOSITION  ***

    CALL dl7srt(1, p, l, dihdi, irc)
    IF (irc == 0) GO TO 240

    !  ***  (D**-1)*H*(D**-1) + ALPHAK*I  IS INDEFINITE -- OVERESTIMATE
    !  ***  SMALLEST EIGENVALUE FOR USE IN UPDATING LK  ***

    j = (irc*(irc+1))/2
    t = l(j)
    l(j) = one
    DO  i = 1, irc
      w(i) = zero
    END DO
    w(irc) = one
    CALL dl7itv(irc, w, l, w)
    t1 = dv2nrm(irc, w)
    lk = alphak - t/t1/t1
    v(dst0) = -lk
    IF (uk < lk) uk = lk
    IF (alphak < lk) GO TO 210

    !  ***  NASTY CASE -- EXACT GERSCHGORIN BOUNDS.  FUDGE LK, UK...

    t = p001 * alphak
    IF (t <= zero) t = p001
    lk = alphak + t
    IF (uk <= lk) uk = lk + t
    GO TO 210

    !  ***  ALPHAK MAKES (D**-1)*H*(D**-1) POSITIVE DEFINITE.
    !  ***  COMPUTE Q = -D*STEP, CHECK FOR CONVERGENCE.  ***

    240  CALL dl7ivm(p, w(q), l, dig)
    gtsta = dd7tpr(p, w(q), w(q))
    CALL dl7itv(p, w(q), l, w(q))
    dst = dv2nrm(p, w(q))
    phi = dst - rad
    IF (phi <= phimax .AND. phi >= phimin) GO TO 270
    IF (phi == oldphi) GO TO 270
    oldphi = phi
    IF (phi < zero) GO TO 330

    !  ***  UNACCEPTABLE ALPHAK -- UPDATE LK, UK, ALPHAK  ***

    250  IF (ka >= kalim) GO TO 270
    !     ***  THE FOLLOWING   MIN IS NECESSARY BECAUSE OF RESTARTS  ***
    IF (phi < zero) uk =   MIN(uk, alphak)
    !     *** KAMIN = 0 ONLY IFF THE GRADIENT VANISHES  ***
    IF (kamin == 0) GO TO 210
    CALL dl7ivm(p, w, l, w(q))
    !     *** THE FOLLOWING, COMMENTED CALCULATION OF ALPHAK IS SOMETIMES
    !     *** SAFER BUT WORSE IN PERFORMANCE...
    !     T1 = DST / DV2NRM(P, W)
    !     ALPHAK = ALPHAK  +  T1 * (PHI/RAD) * T1
    t1 = dv2nrm(p, w)
    alphak = alphak  +  (phi/t1) * (dst/t1) * (dst/rad)
    lk =   MAX(lk, alphak)
    alphak = lk
    GO TO 210

    !  ***  ACCEPTABLE STEP ON FIRST TRY  ***

    260  alphak = zero

    !  ***  SUCCESSFUL STEP IN GENERAL.  COMPUTE STEP = -(D**-1)*Q  ***

    270  DO  i = 1, p
      j = q0 + i
      step(i) = -w(j)/d(i)
    END DO
    v(gtstep) = -gtsta
    v(preduc) = half * ( ABS(alphak)*dst*dst + gtsta)
    GO TO 410


    !  ***  RESTART WITH NEW RADIUS  ***

    290  IF (v(dst0) <= zero .OR. v(dst0) - rad > phimax) GO TO 310

    !     ***  PREPARE TO RETURN NEWTON STEP  ***

    restrt = .true.
    ka = ka + 1
    k = 0
    DO  i = 1, p
      k = k + i
      j = diag0 + i
      dihdi(k) = w(j)
    END DO
    uk = negone
    GO TO 30

    310  kamin = ka + 3
    IF (v(dgnorm) == zero) kamin = 0
    IF (ka == 0) GO TO 50

    dst = w(dstsav)
    alphak =  ABS(v(stppar))
    phi = dst - rad
    t = v(dgnorm)/rad
    uk = t - w(emin)
    IF (v(dgnorm) == zero) uk = uk + p001 + p001*uk
    IF (uk <= zero) uk = p001
    IF (rad > v(rad0)) GO TO 320

    !        ***  SMALLER RADIUS  ***
    lk = zero
    IF (alphak > zero) lk = w(lk0)
    lk =   MAX(lk, t - w(emax))
    IF (v(dst0) > zero) lk =   MAX(lk, (v(dst0)-rad)*w(phipin))
    GO TO 250

    !     ***  BIGGER RADIUS  ***
    320  IF (alphak > zero) uk =   MIN(uk, w(uk0))
    lk =   MAX(zero, -v(dst0), t - w(emax))
    IF (v(dst0) > zero) lk =   MAX(lk, (v(dst0)-rad)*w(phipin))
    GO TO 250

    !  ***  DECIDE WHETHER TO CHECK FOR SPECIAL CASE... IN PRACTICE (FROM
    !  ***  THE STANDPOINT OF THE CALLING OPTIMIZATION CODE) IT SEEMS BEST
    !  ***  NOT TO CHECK UNTIL A FEW ITERATIONS HAVE FAILED -- HENCE THE
    !  ***  TEST ON KAMIN BELOW.

    330  delta = alphak +   MIN(zero, v(dst0))
    twopsi = alphak*dst*dst + gtsta
    IF (ka >= kamin) GO TO 340
    !     *** IF THE TEST IN REF. 2 IS SATISFIED, FALL THROUGH TO HANDLE
    !     *** THE SPECIAL CASE (AS SOON AS THE MORE-SORENSEN TEST DETECTS
    !     *** IT).
    IF (psifac >= big) GO TO 340
    IF (delta >= psifac*twopsi) GO TO 370

    !  ***  CHECK FOR THE SPECIAL CASE OF  H + ALPHA*D**2  (NEARLY)
    !  ***  SINGULAR.  USE ONE STEP OF INVERSE POWER METHOD WITH START
    !  ***  FROM DL7SVN TO OBTAIN APPROXIMATE EIGENVECTOR CORRESPONDING
    !  ***  TO SMALLEST EIGENVALUE OF (D**-1)*H*(D**-1).  DL7SVN RETURNS
    !  ***  X AND W WITH  L*W = X.

    340  t = dl7svn(p, l, w(x), w)

    !     ***  NORMALIZE W  ***
    DO  i = 1, p
      w(i) = t*w(i)
    END DO
    !     ***  COMPLETE CURRENT INV. POWER ITER. -- REPLACE W BY (L**-T)*W.
    CALL dl7itv(p, w, l, w)
    t2 = one/dv2nrm(p, w)
    DO  i = 1, p
      w(i) = t2*w(i)
    END DO
    t = t2 * t

    !  ***  NOW W IS THE DESIRED APPROXIMATE (UNIT) EIGENVECTOR AND
    !  ***  T*X = ((D**-1)*H*(D**-1) + ALPHAK*I)*W.

    sw = dd7tpr(p, w(q), w)
    t1 = (rad + dst) * (rad - dst)
    root =  SQRT(sw*sw + t1)
    IF (sw < zero) root = -root
    si = t1 / (sw + root)

    !  ***  THE ACTUAL TEST FOR THE SPECIAL CASE...

    IF ((t2*si)**2 <= eps*(dst**2 + alphak*radsq)) GO TO 380

    !  ***  UPDATE UPPER BOUND ON SMALLEST EIGENVALUE (WHEN NOT POSITIVE)
    !  ***  (AS RECOMMENDED BY MORE AND SORENSEN) AND CONTINUE...

    IF (v(dst0) <= zero) v(dst0) =   MIN(v(dst0), t2**2 - alphak)
    lk =   MAX(lk, -v(dst0))

    !  ***  CHECK WHETHER WE CAN HOPE TO DETECT THE SPECIAL CASE IN
    !  ***  THE AVAILABLE ARITHMETIC.  ACCEPT STEP AS IT IS IF NOT.

    !     ***  IF NOT YET AVAILABLE, OBTAIN MACHINE DEPENDENT VALUE DGXFAC.
    370  IF (dgxfac == zero) dgxfac = epsfac * dr7mdc(3)

    IF (delta > dgxfac*w(dggdmx)) GO TO 250
    GO TO 270

    !  ***  SPECIAL CASE DETECTED... NEGATE ALPHAK TO INDICATE SPECIAL CASE

    380  alphak = -alphak
    v(preduc) = half * twopsi

    !  ***  ACCEPT CURRENT STEP IF ADDING SI*W WOULD LEAD TO A
    !  ***  FURTHER RELATIVE REDUCTION IN PSI OF LESS THAN V(EPSLON)/3.

    t1 = zero
    t = si*(alphak*sw - half*si*(alphak + t*dd7tpr(p,w(x),w)))
    IF (t < eps*twopsi/six) GO TO 390
    v(preduc) = v(preduc) + t
    dst = rad
    t1 = -si
    390  DO  i = 1, p
      j = q0 + i
      w(j) = t1*w(i) - w(j)
      step(i) = w(j) / d(i)
    END DO
    v(gtstep) = dd7tpr(p, dig, w(q))

    !  ***  SAVE VALUES FOR USE IN A POSSIBLE RESTART  ***

    410  v(dstnrm) = dst
    v(stppar) = alphak
    w(lk0) = lk
    w(uk0) = uk
    v(rad0) = rad
    w(dstsav) = dst

    !     ***  RESTORE DIAGONAL OF DIHDI  ***

    j = 0
    DO  i = 1, p
      j = j + i
      k = diag0 + i
      dihdi(j) = w(k)
    END DO

    !999  RETURN
    RETURN
    RETURN
    !  ***  LAST LINE OF DG7QTS FOLLOWS  ***
    END SUBROUTINE dg7qts

    ! SUBROUTINE ditsum(d, g, iv, liv, lv, p, v, x)
    SUBROUTINE ditsum(inDummy, outDummy)
    !
    ! DSB NOTE:  See note at the beginning.
    ! We are creating a dummy version of this here.
    !  ***  PRINT ITERATION SUMMARY FOR ***SOL (VERSION 2.3)  ***

    !  ***  PARAMETER DECLARATIONS  ***

    INTEGER :: inDummy, outDummy
      outDummy = inDummy
    RETURN
    !  ***  LAST LINE OF DITSUM FOLLOWS  ***
    END
    SUBROUTINE divset(alg, iv, liv, lv, v)

    !  ***  SUPPLY ***SOL (VERSION 2.3) DEFAULT VALUES TO IV AND V  ***

    !  ***  ALG = 1 MEANS REGRESSION CONSTANTS.
    !  ***  ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS.

    INTEGER :: liv, lv
    INTEGER :: alg, iv(liv)
    !      INTEGER ALG
    !      DOUBLE PRECISION IV(LIV)
    DOUBLE PRECISION :: v(lv)

    INTEGER :: i7mdcn
    EXTERNAL i7mdcn,dv7dfl
    ! I7MDCN... RETURNS MACHINE-DEPENDENT INTEGER CONSTANTS.
    ! DV7DFL.... PROVIDES DEFAULT VALUES TO V.

    INTEGER :: alg1, miv, mv
    INTEGER :: miniv(4), minv(4)

    !  ***  SUBSCRIPTS FOR IV  ***

    INTEGER :: algsav, covprt, covreq, dradpr, dtype, hc, ierr, inith,  &
        inits, ipivot, ivneed, lastiv, lastv, lmat, mxfcal,  &
        mxiter, nfcov, ngcov, nvdflt, nvsave, outlev, parprt,  &
        parsav, perm, prunit, qrtyp, rdreq, rmat, solprt, statpr, vneed, vsave, x0prt

    !  ***  IV SUBSCRIPT VALUES  ***

    PARAMETER (algsav=51, covprt=14, covreq=15, dradpr=101, dtype=16,  &
        hc=71, ierr=75, inith=25, inits=25, ipivot=76,  &
        ivneed=3, lastiv=44, lastv=45, lmat=42, mxfcal=17,  &
        mxiter=18, nfcov=52, ngcov=53, nvdflt=50, nvsave=9,  &
        outlev=19, parprt=20, parsav=49, perm=58, prunit=21,  &
        qrtyp=80, rdreq=57, rmat=78, solprt=22, statpr=23,  &
        vneed=4, vsave=60, x0prt=24)
    DATA miniv(1)/82/, miniv(2)/59/, miniv(3)/103/, miniv(4)/103/,  &
        minv(1)/98/, minv(2)/71/, minv(3)/101/, minv(4)/85/

    !-------------------------------  BODY  --------------------------------

    IF (prunit <= liv) iv(prunit) = i7mdcn(1)
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
    END SUBROUTINE divset
    SUBROUTINE dl7itv(n, x, l, y)

    !  ***  SOLVE  (L**T)*X = Y,  WHERE  L  IS AN  N X N  LOWER TRIANGULAR
    !  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
    !  ***  STORAGE.  ***

    INTEGER :: n
    DOUBLE PRECISION :: x(n), l(1), y(n)
    INTEGER :: i, ii, ij, im1, i0, j, np1
    DOUBLE PRECISION :: xi, zero
    PARAMETER (zero=0.d+0)

    DO  i = 1, n
      x(i) = y(i)
    END DO
    np1 = n + 1
    i0 = n*(n+1)/2
    DO  ii = 1, n
      i = np1 - ii
      xi = x(i)/l(i0)
      x(i) = xi
      IF (i <= 1) EXIT
      i0 = i0 - i
      IF (xi == zero) CYCLE
      im1 = i - 1
      DO  j = 1, im1
        ij = i0 + j
        x(j) = x(j) - xi*l(ij)
      END DO
    END DO
    ! 999  RETURN
    RETURN
    !  ***  LAST LINE OF DL7ITV FOLLOWS  ***
    END SUBROUTINE dl7itv
    SUBROUTINE dl7ivm(n, x, l, y)

    !  ***  SOLVE  L*X = Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR
    !  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
    !  ***  STORAGE.  ***

    INTEGER :: n
    DOUBLE PRECISION :: x(n), l(1), y(n)
    DOUBLE PRECISION :: dd7tpr
    EXTERNAL dd7tpr
    INTEGER :: i, j, k
    DOUBLE PRECISION :: t, zero
    PARAMETER (zero=0.d+0)

    DO  k = 1, n
      IF (y(k) /= zero) GO TO 20
      x(k) = zero
    END DO
    GO TO 999
    20   j = k*(k+1)/2
    x(k) = y(k) / l(j)
    IF (k >= n) GO TO 999
    k = k + 1
    DO  i = k, n
      t = dd7tpr(i-1, l(j+1), x)
      j = j + i
      x(i) = (y(i) - t)/l(j)
    END DO
    999  RETURN
    !  ***  LAST LINE OF DL7IVM FOLLOWS  ***
    END SUBROUTINE dl7ivm
    SUBROUTINE dl7mst(d, g, ierr, ipivot, ka, p, qtr, r, step, v, w)

    !  ***  COMPUTE LEVENBERG-MARQUARDT STEP USING MORE-HEBDEN TECHNIQUE  **
    !  ***  NL2SOL VERSION 2.2.  ***

    !  ***  PARAMETER DECLARATIONS  ***

    INTEGER :: ierr, ka, p
    INTEGER :: ipivot(p)
    DOUBLE PRECISION :: d(p), g(p), qtr(p), r(1), step(p), v(21), w(1)
    !     DIMENSION W(P*(P+5)/2 + 4)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ***  PURPOSE  ***

    !        GIVEN THE R MATRIX FROM THE QR DECOMPOSITION OF A JACOBIAN
    !     MATRIX, J, AS WELL AS Q-TRANSPOSE TIMES THE CORRESPONDING
    !     RESIDUAL VECTOR, RESID, THIS SUBROUTINE COMPUTES A LEVENBERG-
    !     MARQUARDT STEP OF APPROXIMATE LENGTH V(RADIUS) BY THE MORE-
    !     TECHNIQUE.

    !  ***  PARAMETER DESCRIPTION  ***

    !      D (IN)  = THE SCALE VECTOR.
    !      G (IN)  = THE GRADIENT VECTOR (J**T)*R.
    !   IERR (I/O) = RETURN CODE FROM QRFACT OR DQ7RGS -- 0 MEANS R HAS
    !             FULL RANK.
    ! IPIVOT (I/O) = PERMUTATION ARRAY FROM QRFACT OR DQ7RGS, WHICH COMPUTE
    !             QR DECOMPOSITIONS WITH COLUMN PIVOTING.
    !     KA (I/O).  KA .LT. 0 ON INPUT MEANS THIS IS THE FIRST CALL ON
    !             DL7MST FOR THE CURRENT R AND QTR.  ON OUTPUT KA CON-
    !             TAINS THE NUMBER OF HEBDEN ITERATIONS NEEDED TO DETERMINE
    !             STEP.  KA = 0 MEANS A GAUSS-NEWTON STEP.
    !      P (IN)  = NUMBER OF PARAMETERS.
    !    QTR (IN)  = (Q**T)*RESID = Q-TRANSPOSE TIMES THE RESIDUAL VECTOR.
    !      R (IN)  = THE R MATRIX, STORED COMPACTLY BY COLUMNS.
    !   STEP (OUT) = THE LEVENBERG-MARQUARDT STEP COMPUTED.
    !      V (I/O) CONTAINS VARIOUS CONSTANTS AND VARIABLES DESCRIBED BELOW.
    !      W (I/O) = WORKSPACE OF LENGTH P*(P+5)/2 + 4.

    !  ***  ENTRIES IN V  ***

    ! V(DGNORM) (I/O) = 2-NORM OF (D**-1)*G.
    ! V(DSTNRM) (I/O) = 2-NORM OF D*STEP.
    ! V(DST0)   (I/O) = 2-NORM OF GAUSS-NEWTON STEP (FOR NONSING. J).
    ! V(EPSLON) (IN) = MAX. REL. ERROR ALLOWED IN TWONORM(R)**2 MINUS
    !             TWONORM(R - J*STEP)**2.  (SEE ALGORITHM NOTES BELOW.)
    ! V(GTSTEP) (OUT) = INNER PRODUCT BETWEEN G AND STEP.
    ! V(NREDUC) (OUT) = HALF THE REDUCTION IN THE SUM OF SQUARES PREDICTED
    !             FOR A GAUSS-NEWTON STEP.
    ! V(PHMNFC) (IN)  = TOL. (TOGETHER WITH V(PHMXFC)) FOR ACCEPTING STEP
    !             (MORE*S SIGMA).  THE ERROR V(DSTNRM) - V(RADIUS) MUST LIE
    !             BETWEEN V(PHMNFC)*V(RADIUS) AND V(PHMXFC)*V(RADIUS).
    ! V(PHMXFC) (IN)  (SEE V(PHMNFC).)
    ! V(PREDUC) (OUT) = HALF THE REDUCTION IN THE SUM OF SQUARES PREDICTED
    !             BY THE STEP RETURNED.
    ! V(RADIUS) (IN)  = RADIUS OF CURRENT (SCALED) TRUST REGION.
    ! V(RAD0)   (I/O) = VALUE OF V(RADIUS) FROM PREVIOUS CALL.
    ! V(STPPAR) (I/O) = MARQUARDT PARAMETER (OR ITS NEGATIVE IF THE SPECIAL
    !             CASE MENTIONED BELOW IN THE ALGORITHM NOTES OCCURS).

    ! NOTE -- SEE DATA STATEMENT BELOW FOR VALUES OF ABOVE SUBSCRIPTS.

    !  ***  USAGE NOTES  ***

    !     IF IT IS DESIRED TO RECOMPUTE STEP USING A DIFFERENT VALUE OF
    !     V(RADIUS), THEN THIS ROUTINE MAY BE RESTARTED BY CALLING IT
    !     WITH ALL PARAMETERS UNCHANGED EXCEPT V(RADIUS).  (THIS EXPLAINS
    !     WHY MANY PARAMETERS ARE LISTED AS I/O).  ON AN INTIIAL CALL (ONE
    !     WITH KA = -1), THE CALLER NEED ONLY HAVE INITIALIZED D, G, KA, P,
    !     QTR, R, V(EPSLON), V(PHMNFC), V(PHMXFC), V(RADIUS), AND V(RAD0).

    !  ***  APPLICATION AND USAGE RESTRICTIONS  ***

    !     THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR LEAST-
    !     SQUARES) PACKAGE (REF. 1).

    !  ***  ALGORITHM NOTES  ***

    !     THIS CODE IMPLEMENTS THE STEP COMPUTATION SCHEME DESCRIBED IN
    !     REFS. 2 AND 4.  FAST GIVENS TRANSFORMATIONS (SEE REF. 3, PP. 60-
    !     62) ARE USED TO COMPUTE STEP WITH A NONZERO MARQUARDT PARAMETER.
    !        A SPECIAL CASE OCCURS IF J IS (NEARLY) SINGULAR AND V(RADIUS)
    !     IS SUFFICIENTLY LARGE.  IN THIS CASE THE STEP RETURNED IS SUCH
    !     THAT  TWONORM(R)**2 - TWONORM(R - J*STEP)**2  DIFFERS FROM ITS
    !     OPTIMAL VALUE BY LESS THAN V(EPSLON) TIMES THIS OPTIMAL VALUE,
    !     WHERE J AND R DENOTE THE ORIGINAL JACOBIAN AND RESIDUAL.  (SEE
    !     REF. 2 FOR MORE DETAILS.)

    !  ***  FUNCTIONS AND SUBROUTINES CALLED  ***

    ! DD7TPR - RETURNS INNER PRODUCT OF TWO VECTORS.
    ! DL7ITV - APPLY INVERSE-TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX.
    ! DL7IVM - APPLY INVERSE OF COMPACT LOWER TRIANG. MATRIX.
    ! DV7CPY  - COPIES ONE VECTOR TO ANOTHER.
    ! DV2NRM - RETURNS 2-NORM OF A VECTOR.

    !  ***  REFERENCES  ***

    ! 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), AN ADAPTIVE
    !             NONLINEAR LEAST-SQUARES ALGORITHM, ACM TRANS. MATH.
    !             SOFTWARE, VOL. 7, NO. 3.
    ! 2.  GAY, D.M. (1981), COMPUTING OPTIMAL LOCALLY CONSTRAINED STEPS,
    !             SIAM J. SCI. STATIST. COMPUTING, VOL. 2, NO. 2, PP.
    !             186-197.
    ! 3.  LAWSON, C.L., AND HANSON, R.J. (1974), SOLVING LEAST SQUARES
    !             PROBLEMS, PRENTICE-HALL, ENGLEWOOD CLIFFS, N.J.
    ! 4.  MORE, J.J. (1978), THE LEVENBERG-MARQUARDT ALGORITHM, IMPLEMEN-
    !             TATION AND THEORY, PP.105-116 OF SPRINGER LECTURE NOTES
    !             IN MATHEMATICS NO. 630, EDITED BY G.A. WATSON, SPRINGER-
    !             VERLAG, BERLIN AND NEW YORK.

    !  ***  GENERAL  ***

    !     CODED BY DAVID M. GAY.
    !     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
    !     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
    !     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND
    !     MCS-7906671.

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: dstsav, i, ip1, i1, j1, k, kalim, l, lk0, phipin,  &
        pp1o2, res, res0, rmat, rmat0, uk0
    DOUBLE PRECISION :: a, adi, alphak, b, dfacsq, dst, dtol, d1, d2,  &
        lk, oldphi, phi, phimax, phimin, psifac, rad,  &
        si, sj, sqrtak, t, twopsi, uk, wl

    !     ***  CONSTANTS  ***
    DOUBLE PRECISION :: dfac, eight, half, negone, one, p001, three, ttol, zero
    DOUBLE PRECISION :: big

    !  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***

    DOUBLE PRECISION :: dd7tpr, dl7svn, dr7mdc, dv2nrm
    EXTERNAL dd7tpr, dl7itv, dl7ivm, dl7svn, dr7mdc,dv7cpy, dv2nrm

    !  ***  SUBSCRIPTS FOR V  ***

    INTEGER :: dgnorm, dstnrm, dst0, epslon, gtstep, nreduc, phmnfc,  &
        phmxfc, preduc, radius, rad0, stppar
    PARAMETER (dgnorm=1, dstnrm=2, dst0=3, epslon=19, gtstep=4,  &
        nreduc=6, phmnfc=20, phmxfc=21, preduc=7, radius=8, rad0=9, stppar=5)

    PARAMETER (dfac=256.d+0, eight=8.d+0, half=0.5D+0, negone=-1.d+0,  &
        one=1.d+0, p001=1.d-3, three=3.d+0, ttol=2.5D+0, zero=0.d+0)
    SAVE big
    DATA big/0.d+0/

    ! DSB NOTE:  The following lines were added to get rid of compiler warnings.
    alphak = one
    psifac = one
    !  ***  BODY  ***
    !     ***  FOR USE IN RECOMPUTING STEP, THE FINAL VALUES OF LK AND UK,
    !     ***  THE INVERSE DERIVATIVE OF MORE*S PHI AT 0 (FOR NONSING. J)
    !     ***  AND THE VALUE RETURNED AS V(DSTNRM) ARE STORED AT W(LK0),
    !     ***  W(UK0), W(PHIPIN), AND W(DSTSAV) RESPECTIVELY.
    lk0 = p + 1
    phipin = lk0 + 1
    uk0 = phipin + 1
    dstsav = uk0 + 1
    rmat0 = dstsav
    !     ***  A COPY OF THE R-MATRIX FROM THE QR DECOMPOSITION OF J IS
    !     ***  STORED IN W STARTING AT W(RMAT), AND A COPY OF THE RESIDUAL
    !     ***  VECTOR IS STORED IN W STARTING AT W(RES).  THE LOOPS BELOW
    !     ***  THAT UPDATE THE QR DECOMP. FOR A NONZERO MARQUARDT PARAMETER
    !     ***  WORK ON THESE COPIES.
    rmat = rmat0 + 1
    pp1o2 = p * (p + 1) / 2
    res0 = pp1o2 + rmat0
    res = res0 + 1
    rad = v(radius)
    IF (rad > zero)  &
        psifac = v(epslon)/((eight*(v(phmnfc) + one) + three) * rad**2)
    IF (big <= zero) big = dr7mdc(6)
    phimax = v(phmxfc) * rad
    phimin = v(phmnfc) * rad
    !     ***  DTOL, DFAC, AND DFACSQ ARE USED IN RESCALING THE FAST GIVENS
    !     ***  REPRESENTATION OF THE UPDATED QR DECOMPOSITION.
    dtol = one/dfac
    dfacsq = dfac*dfac
    !     ***  OLDPHI IS USED TO DETECT LIMITS OF NUMERICAL ACCURACY.  IF
    !     ***  WE RECOMPUTE STEP AND IT DOES NOT CHANGE, THEN WE ACCEPT IT.
    oldphi = zero
    lk = zero
    uk = zero
    kalim = ka + 12

    !  ***  START OR RESTART, DEPENDING ON KA  ***

    IF (ka < 0) THEN
      GO TO    10
    ELSE IF (ka == 0) THEN
      GO TO    20
    ELSE
      GO TO   370
    END IF

    !  ***  FRESH START -- COMPUTE V(NREDUC)  ***

    10   ka = 0
    kalim = 12
    k = p
    IF (ierr /= 0) k = IABS(ierr) - 1
    v(nreduc) = half*dd7tpr(k, qtr, qtr)

    !  ***  SET UP TO TRY INITIAL GAUSS-NEWTON STEP  ***

    20   v(dst0) = negone
    IF (ierr /= 0) GO TO 90
    t = dl7svn(p, r, step, w(res))
    IF (t >= one) GO TO 30
    IF (dv2nrm(p, qtr) >= big*t) GO TO 90

    !  ***  COMPUTE GAUSS-NEWTON STEP  ***

    !     ***  NOTE -- THE R-MATRIX IS STORED COMPACTLY BY COLUMNS IN
    !     ***  R(1), R(2), R(3), ...  IT IS THE TRANSPOSE OF A
    !     ***  LOWER TRIANGULAR MATRIX STORED COMPACTLY BY ROWS, AND WE
    !     ***  TREAT IT AS SUCH WHEN USING DL7ITV AND DL7IVM.
    30   CALL dl7itv(p, w, r, qtr)
    !     ***  TEMPORARILY STORE PERMUTED -D*STEP IN STEP.
    DO  i = 1, p
      j1 = ipivot(i)
      step(i) = d(j1)*w(i)
    END DO
    dst = dv2nrm(p, step)
    v(dst0) = dst
    phi = dst - rad
    IF (phi <= phimax) GO TO 410
    !     ***  IF THIS IS A RESTART, GO TO 110  ***
    IF (ka > 0) GO TO 110

    !  ***  GAUSS-NEWTON STEP WAS UNACCEPTABLE.  COMPUTE L0  ***

    DO  i = 1, p
      j1 = ipivot(i)
      step(i) = d(j1)*(step(i)/dst)
    END DO
    CALL dl7ivm(p, step, r, step)
    t = one / dv2nrm(p, step)
    w(phipin) = (t/rad)*t
    lk = phi*w(phipin)

    !  ***  COMPUTE U0  ***

    90   DO  i = 1, p
      w(i) = g(i)/d(i)
    END DO
    v(dgnorm) = dv2nrm(p, w)
    uk = v(dgnorm)/rad
    IF (uk <= zero) GO TO 390

    !     ***  ALPHAK WILL BE USED AS THE CURRENT MARQUARDT PARAMETER.  WE
    !     ***  USE MORE*S SCHEME FOR INITIALIZING IT.

    alphak =  ABS(v(stppar)) * v(rad0)/rad
    alphak =   MIN(uk,   MAX(alphak, lk))


    !  ***  TOP OF LOOP -- INCREMENT KA, COPY R TO RMAT, QTR TO RES  ***

    110  ka = ka + 1
    CALL dv7cpy(pp1o2, w(rmat), r)
    CALL dv7cpy(p, w(res), qtr)

    !  ***  SAFEGUARD ALPHAK AND INITIALIZE FAST GIVENS SCALE VECTOR.  ***

    IF (alphak <= zero .OR. alphak < lk .OR. alphak >= uk)  &
        alphak = uk *   MAX(p001,  SQRT(lk/uk))
    IF (alphak <= zero) alphak = half * uk
    sqrtak =  SQRT(alphak)
    DO  i = 1, p
      w(i) = one
    END DO

    !  ***  ADD ALPHAK*D AND UPDATE QR DECOMP. USING FAST GIVENS TRANS.  ***

    DO  i = 1, p
    !        ***  GENERATE, APPLY 1ST GIVENS TRANS. FOR ROW I OF ALPHAK*D.
    !        ***  (USE STEP TO STORE TEMPORARY ROW)  ***
      l = i*(i+1)/2 + rmat0
      wl = w(l)
      d2 = one
      d1 = w(i)
      j1 = ipivot(i)
      adi = sqrtak*d(j1)
      IF (adi >=  ABS(wl)) GO TO 150
      130     a = adi/wl
      b = d2*a/d1
      t = a*b + one
      IF (t > ttol) GO TO 150
      w(i) = d1/t
      d2 = d2/t
      w(l) = t*wl
      a = -a
      DO  j1 = i, p
        l = l + j1
        step(j1) = a*w(l)
      END DO
      GO TO 170

      150     b = wl/adi
      a = d1*b/d2
      t = a*b + one
      IF (t > ttol) GO TO 130
      w(i) = d2/t
      d2 = d1/t
      w(l) = t*adi
      DO  j1 = i, p
        l = l + j1
        wl = w(l)
        step(j1) = -wl
        w(l) = a*wl
      END DO

      170     IF (i == p) GO TO 280

    !        ***  NOW USE GIVENS TRANS. TO ZERO ELEMENTS OF TEMP. ROW  ***

      ip1 = i + 1
      DO  i1 = ip1, p
        l = i1*(i1+1)/2 + rmat0
        wl = w(l)
        si = step(i1-1)
        d1 = w(i1)

    !             ***  RESCALE ROW I1 IF NECESSARY  ***

        IF (d1 >= dtol) GO TO 190
        d1 = d1*dfacsq
        wl = wl/dfac
        k = l
        DO  j1 = i1, p
          k = k + j1
          w(k) = w(k)/dfac
        END DO

    !             ***  USE GIVENS TRANS. TO ZERO NEXT ELEMENT OF TEMP. ROW

        190          IF ( ABS(si) >  ABS(wl)) GO TO 220
        IF (si == zero) CYCLE
        200          a = si/wl
        b = d2*a/d1
        t = a*b + one
        IF (t > ttol) GO TO 220
        w(l) = t*wl
        w(i1) = d1/t
        d2 = d2/t
        DO  j1 = i1, p
          l = l + j1
          wl = w(l)
          sj = step(j1)
          w(l) = wl + b*sj
          step(j1) = sj - a*wl
        END DO
        GO TO 240

        220          b = wl/si
        a = d1*b/d2
        t = a*b + one
        IF (t > ttol) GO TO 200
        w(i1) = d2/t
        d2 = d1/t
        w(l) = t*si
        DO  j1 = i1, p
          l = l + j1
          wl = w(l)
          sj = step(j1)
          w(l) = a*wl + sj
          step(j1) = b*sj - wl
        END DO

    !             ***  RESCALE TEMP. ROW IF NECESSARY  ***

        240          IF (d2 >= dtol) CYCLE
        d2 = d2*dfacsq
        DO  k = i1, p
          step(k) = step(k)/dfac
        END DO
      END DO
    END DO

    !  ***  COMPUTE STEP  ***

    280  CALL dl7itv(p, w(res), w(rmat), w(res))
    !     ***  RECOVER STEP AND STORE PERMUTED -D*STEP AT W(RES)  ***
    DO  i = 1, p
      j1 = ipivot(i)
      k = res0 + i
      t = w(k)
      step(j1) = -t
      w(k) = t*d(j1)
    END DO
    dst = dv2nrm(p, w(res))
    phi = dst - rad
    IF (phi <= phimax .AND. phi >= phimin) GO TO 430
    IF (oldphi == phi) GO TO 430
    oldphi = phi

    !  ***  CHECK FOR (AND HANDLE) SPECIAL CASE  ***

    IF (phi > zero) GO TO 310
    IF (ka >= kalim) GO TO 430
    twopsi = alphak*dst*dst - dd7tpr(p, step, g)
    IF (alphak >= twopsi*psifac) GO TO 310
    v(stppar) = -alphak
    GO TO 440

    !  ***  UNACCEPTABLE STEP -- UPDATE LK, UK, ALPHAK, AND TRY AGAIN  ***

    300  IF (phi < zero) uk =   MIN(uk, alphak)
    GO TO 320
    310  IF (phi < zero) uk = alphak
    320  DO  i = 1, p
      j1 = ipivot(i)
      k = res0 + i
      step(i) = d(j1) * (w(k)/dst)
    END DO
    CALL dl7ivm(p, step, w(rmat), step)
    DO  i = 1, p
      step(i) = step(i) /  SQRT(w(i))
    END DO
    t = one / dv2nrm(p, step)
    alphak = alphak + t*phi*t/rad
    lk =   MAX(lk, alphak)
    alphak = lk
    GO TO 110

    !  ***  RESTART  ***

    370  lk = w(lk0)
    uk = w(uk0)
    IF (v(dst0) > zero .AND. v(dst0) - rad <= phimax) GO TO 20
    alphak =  ABS(v(stppar))
    dst = w(dstsav)
    phi = dst - rad
    t = v(dgnorm)/rad
    IF (rad > v(rad0)) GO TO 380

    !        ***  SMALLER RADIUS  ***
    uk = t
    IF (alphak <= zero) lk = zero
    IF (v(dst0) > zero) lk =   MAX(lk, (v(dst0)-rad)*w(phipin))
    GO TO 300

    !     ***  BIGGER RADIUS  ***
    380  IF (alphak <= zero .OR. uk > t) uk = t
    lk = zero
    IF (v(dst0) > zero) lk =   MAX(lk, (v(dst0)-rad)*w(phipin))
    GO TO 300

    !  ***  SPECIAL CASE -- RAD .LE. 0 OR (G = 0 AND J IS SINGULAR)  ***

    390  v(stppar) = zero
    dst = zero
    lk = zero
    uk = zero
    v(gtstep) = zero
    v(preduc) = zero
    DO  i = 1, p
      step(i) = zero
    END DO
    GO TO 450

    !  ***  ACCEPTABLE GAUSS-NEWTON STEP -- RECOVER STEP FROM W  ***

    410  alphak = zero
    DO  i = 1, p
      j1 = ipivot(i)
      step(j1) = -w(i)
    END DO

    !  ***  SAVE VALUES FOR USE IN A POSSIBLE RESTART  ***

    430  v(stppar) = alphak
    440  v(gtstep) =   MIN(dd7tpr(p,step,g), zero)
    v(preduc) = half * (alphak*dst*dst - v(gtstep))
    450  v(dstnrm) = dst
    w(dstsav) = dst
    w(lk0) = lk
    w(uk0) = uk
    v(rad0) = rad

    ! 999  RETURN
    RETURN

    !  ***  LAST LINE OF DL7MST FOLLOWS  ***
    END SUBROUTINE dl7mst
    SUBROUTINE dl7sqr(n, a, l)

    !  ***  COMPUTE  A = LOWER TRIANGLE OF  L*(L**T),  WITH BOTH
    !  ***  L  AND  A  STORED COMPACTLY BY ROWS.  (BOTH MAY OCCUPY THE
    !  ***  SAME STORAGE.

    !  ***  PARAMETERS  ***

    INTEGER :: n
    DOUBLE PRECISION :: a(1), l(1)
    !     DIMENSION A(N*(N+1)/2), L(N*(N+1)/2)

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: i, ii, ij, ik, ip1, i0, j, jj, jk, j0, k, np1
    DOUBLE PRECISION :: t

    np1 = n + 1
    i0 = n*(n+1)/2
    DO  ii = 1, n
      i = np1 - ii
      ip1 = i + 1
      i0 = i0 - i
      j0 = i*(i+1)/2
      DO  jj = 1, i
        j = ip1 - jj
        j0 = j0 - j
        t = 0.0D0
        DO  k = 1, j
          ik = i0 + k
          jk = j0 + k
          t = t + l(ik)*l(jk)
        END DO
        ij = i0 + j
        a(ij) = t
      END DO
    END DO
    ! 999  RETURN
    RETURN
    END SUBROUTINE dl7sqr
    SUBROUTINE dl7srt(n1, n, l, a, irc)

    !  ***  COMPUTE ROWS N1 THROUGH N OF THE CHOLESKY FACTOR  L  OF
    !  ***  A = L*(L**T),  WHERE  L  AND THE LOWER TRIANGLE OF  A  ARE BOTH
    !  ***  STORED COMPACTLY BY ROWS (AND MAY OCCUPY THE SAME STORAGE).
    !  ***  IRC = 0 MEANS ALL WENT WELL.  IRC = J MEANS THE LEADING
    !  ***  PRINCIPAL  J X J  SUBMATRIX OF  A  IS NOT POSITIVE DEFINITE --
    !  ***  AND  L(J*(J+1)/2)  CONTAINS THE (NONPOS.) REDUCED J-TH DIAGONAL.

    !  ***  PARAMETERS  ***

    INTEGER :: n1, n, irc
    DOUBLE PRECISION :: l(1), a(1)
    !     DIMENSION L(N*(N+1)/2), A(N*(N+1)/2)

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: i, ij, ik, im1, i0, j, jk, jm1, j0, k
    DOUBLE PRECISION :: t, td, zero

    PARAMETER (zero=0.d+0)

    !  ***  BODY  ***

    i0 = n1 * (n1 - 1) / 2
    DO  i = n1, n
      td = zero
      IF (i == 1) GO TO 40
      j0 = 0
      im1 = i - 1
      DO  j = 1, im1
        t = zero
        IF (j == 1) GO TO 20
        jm1 = j - 1
        DO  k = 1, jm1
          ik = i0 + k
          jk = j0 + k
          t = t + l(ik)*l(jk)
        END DO
        20           ij = i0 + j
        j0 = j0 + j
        t = (a(ij) - t) / l(j0)
        l(ij) = t
        td = td + t*t
      END DO
      40      i0 = i0 + i
      t = a(i0) - td
      IF (t <= zero) GO TO 60
      l(i0) =  SQRT(t)
    END DO

    irc = 0
    GO TO 999

    60   l(i0) = t
    irc = i

    999  RETURN

    !  ***  LAST LINE OF DL7SRT  ***
    END SUBROUTINE dl7srt
    DOUBLE PRECISION FUNCTION dl7svn(p, l, x, y)

    !  ***  ESTIMATE SMALLEST SING. VALUE OF PACKED LOWER TRIANG. MATRIX L

    !  ***  PARAMETER DECLARATIONS  ***

    INTEGER :: p
    DOUBLE PRECISION :: l(1), x(p), y(p)
    !     DIMENSION L(P*(P+1)/2)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ***  PURPOSE  ***

    !     THIS FUNCTION RETURNS A GOOD OVER-ESTIMATE OF THE SMALLEST
    !     SINGULAR VALUE OF THE PACKED LOWER TRIANGULAR MATRIX L.

    !  ***  PARAMETER DESCRIPTION  ***

    !  P (IN)  = THE ORDER OF L.  L IS A  P X P  LOWER TRIANGULAR MATRIX.
    !  L (IN)  = ARRAY HOLDING THE ELEMENTS OF  L  IN ROW ORDER, I.E.
    !             L(1,1), L(2,1), L(2,2), L(3,1), L(3,2), L(3,3), ETC.
    !  X (OUT) IF DL7SVN RETURNS A POSITIVE VALUE, THEN X IS A NORMALIZED
    !             APPROXIMATE LEFT SINGULAR VECTOR CORRESPONDING TO THE
    !             SMALLEST SINGULAR VALUE.  THIS APPROXIMATION MAY BE VERY
    !             CRUDE.  IF DL7SVN RETURNS ZERO, THEN SOME COMPONENTS OF X
    !             ARE ZERO AND THE REST RETAIN THEIR INPUT VALUES.
    !  Y (OUT) IF DL7SVN RETURNS A POSITIVE VALUE, THEN Y = (L**-1)*X IS AN
    !             UNNORMALIZED APPROXIMATE RIGHT SINGULAR VECTOR CORRESPOND-
    !             ING TO THE SMALLEST SINGULAR VALUE.  THIS APPROXIMATION
    !             MAY BE CRUDE.  IF DL7SVN RETURNS ZERO, THEN Y RETAINS ITS
    !             INPUT VALUE.  THE CALLER MAY PASS THE SAME VECTOR FOR X
    !             AND Y (NONSTANDARD FORTRAN USAGE), IN WHICH CASE Y OVER-
    !             WRITES X (FOR NONZERO DL7SVN RETURNS).

    !  ***  ALGORITHM NOTES  ***

    !     THE ALGORITHM IS BASED ON (1), WITH THE ADDITIONAL PROVISION THAT
    !     DL7SVN = 0 IS RETURNED IF THE SMALLEST DIAGONAL ELEMENT OF L
    !     (IN MAGNITUDE) IS NOT MORE THAN THE UNIT ROUNDOFF TIMES THE
    !     LARGEST.  THE ALGORITHM USES A RANDOM NUMBER GENERATOR PROPOSED
    !     IN (4), WHICH PASSES THE SPECTRAL TEST WITH FLYING COLORS -- SEE
    !     (2) AND (3).

    !  ***  SUBROUTINES AND FUNCTIONS CALLED  ***

    !        DV2NRM - FUNCTION, RETURNS THE 2-NORM OF A VECTOR.

    !  ***  REFERENCES  ***

    !     (1) CLINE, A., MOLER, C., STEWART, G., AND WILKINSON, J.H.(1977),
    !         AN ESTIMATE FOR THE CONDITION NUMBER OF A MATRIX, REPORT
    !         TM-310, APPLIED MATH. DIV., ARGONNE NATIONAL LABORATORY.

    !     (2) HOAGLIN, D.C. (1976), THEORETICAL PROPERTIES OF CONGRUENTIAL
    !         RANDOM-NUMBER GENERATORS --  AN EMPIRICAL VIEW,
    !         MEMORANDUM NS-340, DEPT. OF STATISTICS, HARVARD UNIV.

    !     (3) KNUTH, D.E. (1969), THE ART OF COMPUTER PROGRAMMING, VOL. 2
    !         (SEMINUMERICAL ALGORITHMS), ADDISON-WESLEY, READING, MASS.

    !     (4) SMITH, C.S. (1971), MULTIPLICATIVE PSEUDO-RANDOM NUMBER
    !         GENERATORS WITH PRIME MODULUS, J. ASSOC. COMPUT. MACH. 18,
    !         PP. 586-593.

    !  ***  HISTORY  ***

    !     DESIGNED AND CODED BY DAVID M. GAY (WINTER 1977/SUMMER 1978).

    !  ***  GENERAL  ***

    !     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
    !     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
    !     MCS-7600324, DCR75-10143, 76-14311DSS, AND MCS76-11989.

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: i, ii, ix, j, ji, jj, jjj, jm1, j0, pm1
    DOUBLE PRECISION :: b, sminus, splus, t, xminus, xplus

    !  ***  CONSTANTS  ***

    DOUBLE PRECISION :: half, one, r9973, zero

    !  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***

    DOUBLE PRECISION :: dd7tpr, dv2nrm
    EXTERNAL dd7tpr, dv2nrm,dv2axy

    PARAMETER (half=0.5D+0, one=1.d+0, r9973=9973.d+0, zero=0.d+0)

    !  ***  BODY  ***

    ix = 2
    pm1 = p - 1

    !  ***  FIRST CHECK WHETHER TO RETURN DL7SVN = 0 AND INITIALIZE X  ***

    ii = 0
    j0 = p*pm1/2
    jj = j0 + p
    IF (l(jj) == zero) GO TO 110
    ix = MOD(3432*ix, 9973)
    b = half*(one + FLOAT(ix)/r9973)
    xplus = b / l(jj)
    x(p) = xplus
    IF (p <= 1) GO TO 60
    DO  i = 1, pm1
      ii = ii + i
      IF (l(ii) == zero) GO TO 110
      ji = j0 + i
      x(i) = xplus * l(ji)
    END DO

    !  ***  SOLVE (L**T)*X = B, WHERE THE COMPONENTS OF B HAVE RANDOMLY
    !  ***  CHOSEN MAGNITUDES IN (.5,1) WITH SIGNS CHOSEN TO MAKE X LARGE.

    !     DO J = P-1 TO 1 BY -1...
    DO  jjj = 1, pm1
      j = p - jjj
    !       ***  DETERMINE X(J) IN THIS ITERATION. NOTE FOR I = 1,2,...,J
    !       ***  THAT X(I) HOLDS THE CURRENT PARTIAL SUM FOR ROW I.
      ix = MOD(3432*ix, 9973)
      b = half*(one + FLOAT(ix)/r9973)
      xplus = (b - x(j))
      xminus = (-b - x(j))
      splus =  ABS(xplus)
      sminus =  ABS(xminus)
      jm1 = j - 1
      j0 = j*jm1/2
      jj = j0 + j
      xplus = xplus/l(jj)
      xminus = xminus/l(jj)
      IF (jm1 == 0) GO TO 30
      DO  i = 1, jm1
        ji = j0 + i
        splus = splus +  ABS(x(i) + l(ji)*xplus)
        sminus = sminus +  ABS(x(i) + l(ji)*xminus)
      END DO
      30      IF (sminus > splus) xplus = xminus
      x(j) = xplus
    !       ***  UPDATE PARTIAL SUMS  ***
      IF (jm1 > 0) CALL dv2axy(jm1, x, xplus, l(j0+1), x)
    END DO

    !  ***  NORMALIZE X  ***

    60   t = one/dv2nrm(p, x)
    DO  i = 1, p
      x(i) = t*x(i)
    END DO

    !  ***  SOLVE L*Y = X AND RETURN DL7SVN = 1/TWONORM(Y)  ***

    DO  j = 1, p
      jm1 = j - 1
      j0 = j*jm1/2
      jj = j0 + j
      t = zero
      IF (jm1 > 0) t = dd7tpr(jm1, l(j0+1), y)
      y(j) = (x(j) - t) / l(jj)
    END DO

    dl7svn = one/dv2nrm(p, y)
    GO TO 999

    110  dl7svn = zero
    999  RETURN
    !  ***  LAST LINE OF DL7SVN FOLLOWS  ***
    END FUNCTION dl7svn
    DOUBLE PRECISION FUNCTION dl7svx(p, l, x, y)

    !  ***  ESTIMATE LARGEST SING. VALUE OF PACKED LOWER TRIANG. MATRIX L

    !  ***  PARAMETER DECLARATIONS  ***

    INTEGER :: p
    DOUBLE PRECISION :: l(1), x(p), y(p)
    !     DIMENSION L(P*(P+1)/2)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ***  PURPOSE  ***

    !     THIS FUNCTION RETURNS A GOOD UNDER-ESTIMATE OF THE LARGEST
    !     SINGULAR VALUE OF THE PACKED LOWER TRIANGULAR MATRIX L.

    !  ***  PARAMETER DESCRIPTION  ***

    !  P (IN)  = THE ORDER OF L.  L IS A  P X P  LOWER TRIANGULAR MATRIX.
    !  L (IN)  = ARRAY HOLDING THE ELEMENTS OF  L  IN ROW ORDER, I.E.
    !             L(1,1), L(2,1), L(2,2), L(3,1), L(3,2), L(3,3), ETC.
    !  X (OUT) IF DL7SVX RETURNS A POSITIVE VALUE, THEN X = (L**T)*Y IS AN
    !             (UNNORMALIZED) APPROXIMATE RIGHT SINGULAR VECTOR
    !             CORRESPONDING TO THE LARGEST SINGULAR VALUE.  THIS
    !             APPROXIMATION MAY BE CRUDE.
    !  Y (OUT) IF DL7SVX RETURNS A POSITIVE VALUE, THEN Y = L*X IS A
    !             NORMALIZED APPROXIMATE LEFT SINGULAR VECTOR CORRESPOND-
    !             ING TO THE LARGEST SINGULAR VALUE.  THIS APPROXIMATION
    !             MAY BE VERY CRUDE.  THE CALLER MAY PASS THE SAME VECTOR
    !             FOR X AND Y (NONSTANDARD FORTRAN USAGE), IN WHICH CASE X
    !             OVER-WRITES Y.

    !  ***  ALGORITHM NOTES  ***

    !     THE ALGORITHM IS BASED ON ANALOGY WITH (1).  IT USES A
    !     RANDOM NUMBER GENERATOR PROPOSED IN (4), WHICH PASSES THE
    !     SPECTRAL TEST WITH FLYING COLORS -- SEE (2) AND (3).

    !  ***  SUBROUTINES AND FUNCTIONS CALLED  ***

    !        DV2NRM - FUNCTION, RETURNS THE 2-NORM OF A VECTOR.

    !  ***  REFERENCES  ***

    !     (1) CLINE, A., MOLER, C., STEWART, G., AND WILKINSON, J.H.(1977),
    !         AN ESTIMATE FOR THE CONDITION NUMBER OF A MATRIX, REPORT
    !         TM-310, APPLIED MATH. DIV., ARGONNE NATIONAL LABORATORY.

    !     (2) HOAGLIN, D.C. (1976), THEORETICAL PROPERTIES OF CONGRUENTIAL
    !         RANDOM-NUMBER GENERATORS --  AN EMPIRICAL VIEW,
    !         MEMORANDUM NS-340, DEPT. OF STATISTICS, HARVARD UNIV.

    !     (3) KNUTH, D.E. (1969), THE ART OF COMPUTER PROGRAMMING, VOL. 2
    !         (SEMINUMERICAL ALGORITHMS), ADDISON-WESLEY, READING, MASS.

    !     (4) SMITH, C.S. (1971), MULTIPLICATIVE PSEUDO-RANDOM NUMBER
    !         GENERATORS WITH PRIME MODULUS, J. ASSOC. COMPUT. MACH. 18,
    !         PP. 586-593.

    !  ***  HISTORY  ***

    !     DESIGNED AND CODED BY DAVID M. GAY (WINTER 1977/SUMMER 1978).

    !  ***  GENERAL  ***

    !     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
    !     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
    !     MCS-7600324, DCR75-10143, 76-14311DSS, AND MCS76-11989.

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: i, ix, j, ji, jj, jjj, jm1, j0, pm1, pplus1
    DOUBLE PRECISION :: b, blji, sminus, splus, t, yi

    !  ***  CONSTANTS  ***

    DOUBLE PRECISION :: half, one, r9973, zero

    !  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***

    DOUBLE PRECISION :: dd7tpr, dv2nrm
    EXTERNAL dd7tpr, dv2nrm,dv2axy

    PARAMETER (half=0.5D+0, one=1.d+0, r9973=9973.d+0, zero=0.d+0)

    !  ***  BODY  ***

    ix = 2
    pplus1 = p + 1
    pm1 = p - 1

    !  ***  FIRST INITIALIZE X TO PARTIAL SUMS  ***

    j0 = p*pm1/2
    jj = j0 + p
    ix = MOD(3432*ix, 9973)
    b = half*(one + FLOAT(ix)/r9973)
    x(p) = b * l(jj)
    IF (p <= 1) GO TO 40
    DO  i = 1, pm1
      ji = j0 + i
      x(i) = b * l(ji)
    END DO

    !  ***  COMPUTE X = (L**T)*B, WHERE THE COMPONENTS OF B HAVE RANDOMLY
    !  ***  CHOSEN MAGNITUDES IN (.5,1) WITH SIGNS CHOSEN TO MAKE X LARGE.

    !     DO J = P-1 TO 1 BY -1...
    DO  jjj = 1, pm1
      j = p - jjj
    !       ***  DETERMINE X(J) IN THIS ITERATION. NOTE FOR I = 1,2,...,J
    !       ***  THAT X(I) HOLDS THE CURRENT PARTIAL SUM FOR ROW I.
      ix = MOD(3432*ix, 9973)
      b = half*(one + FLOAT(ix)/r9973)
      jm1 = j - 1
      j0 = j*jm1/2
      splus = zero
      sminus = zero
      DO  i = 1, j
        ji = j0 + i
        blji = b * l(ji)
        splus = splus +  ABS(blji + x(i))
        sminus = sminus +  ABS(blji - x(i))
      END DO
      IF (sminus > splus) b = -b
      x(j) = zero
    !        ***  UPDATE PARTIAL SUMS  ***
      CALL dv2axy(j, x, b, l(j0+1), x)
    END DO

    !  ***  NORMALIZE X  ***

    40   t = dv2nrm(p, x)
    IF (t <= zero) GO TO 80
    t = one / t
    DO  i = 1, p
      x(i) = t*x(i)
    END DO

    !  ***  COMPUTE L*X = Y AND RETURN SVMAX = TWONORM(Y)  ***

    DO  jjj = 1, p
      j = pplus1 - jjj
      ji = j*(j-1)/2 + 1
      y(j) = dd7tpr(j, l(ji), x)
    END DO

    !  ***  NORMALIZE Y AND SET X = (L**T)*Y  ***

    t = one / dv2nrm(p, y)
    ji = 1
    DO  i = 1, p
      yi = t * y(i)
      x(i) = zero
      CALL dv2axy(i, x, yi, l(ji), x)
      ji = ji + i
    END DO
    dl7svx = dv2nrm(p, x)
    GO TO 999

    80   dl7svx = zero

    999  RETURN
    !  ***  LAST LINE OF DL7SVX FOLLOWS  ***
    END FUNCTION dl7svx
    SUBROUTINE dl7tvm(n, x, l, y)

    !  ***  COMPUTE  X = (L**T)*Y, WHERE  L  IS AN  N X N  LOWER
    !  ***  TRIANGULAR MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY
    !  ***  OCCUPY THE SAME STORAGE.  ***

    INTEGER :: n
    DOUBLE PRECISION :: x(n), l(1), y(n)
    !     DIMENSION L(N*(N+1)/2)
    INTEGER :: i, ij, i0, j
    DOUBLE PRECISION :: yi, zero
    PARAMETER (zero=0.d+0)

    i0 = 0
    DO  i = 1, n
      yi = y(i)
      x(i) = zero
      DO  j = 1, i
        ij = i0 + j
        x(j) = x(j) + yi*l(ij)
      END DO
      i0 = i0 + i
    END DO
    ! 999  RETURN
    RETURN
    !  ***  LAST LINE OF DL7TVM FOLLOWS  ***
    END SUBROUTINE dl7tvm
    SUBROUTINE dl7vml(n, x, l, y)

    !  ***  COMPUTE  X = L*Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR
    !  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
    !  ***  STORAGE.  ***

    INTEGER :: n
    DOUBLE PRECISION :: x(n), l(1), y(n)
    !     DIMENSION L(N*(N+1)/2)
    INTEGER :: i, ii, ij, i0, j, np1
    DOUBLE PRECISION :: t, zero
    PARAMETER (zero=0.d+0)

    np1 = n + 1
    i0 = n*(n+1)/2
    DO  ii = 1, n
      i = np1 - ii
      i0 = i0 - i
      t = zero
      DO  j = 1, i
        ij = i0 + j
        t = t + l(ij)*y(j)
      END DO
      x(i) = t
    END DO
    ! 999  RETURN
    RETURN
    !  ***  LAST LINE OF DL7VML FOLLOWS  ***
    END SUBROUTINE dl7vml
    SUBROUTINE do7prd(l, ls, p, s, w, y, z)

    !  ***  FOR I = 1..L, SET S = S + W(I)*Y(.,I)*(Z(.,I)**T), I.E.,
    !  ***        ADD W(I) TIMES THE OUTER PRODUCT OF Y(.,I) AND Z(.,I).

    INTEGER :: l, ls, p
    DOUBLE PRECISION :: s(ls), w(l), y(p,l), z(p,l)
    !     DIMENSION S(P*(P+1)/2)

    INTEGER :: i, j, k, m
    DOUBLE PRECISION :: wk, yi, zero
    DATA zero/0.d+0/

    DO  k = 1, l
      wk = w(k)
      IF (wk == zero) CYCLE
      m = 1
      DO  i = 1, p
        yi = wk * y(i,k)
        DO  j = 1, i
          s(m) = s(m) + yi*z(j,k)
          m = m + 1
        END DO
      END DO
    END DO

    ! 999  RETURN
    RETURN
    !  ***  LAST LINE OF DO7PRD FOLLOWS  ***
    END SUBROUTINE do7prd
    ! DSB NOTE
    ! SUBROUTINE dparck removed
    SUBROUTINE dq7adr(p, qtr, rmat, w, y)

    !  ***  ADD ROW W TO QR FACTORIZATION WITH R MATRIX RMAT AND
    !  ***  Q**T * RESIDUAL = QTR.  Y = NEW COMPONENT OF RESIDUAL
    !  ***  CORRESPONDING TO W.

    INTEGER :: p
    DOUBLE PRECISION :: qtr(p), rmat(1), w(p), y
    !     DIMENSION RMAT(P*(P+1)/2)

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: i, ii, ij, ip1, j
    DOUBLE PRECISION :: ri, rw, t, u1, u2, v, wi, wr

    DOUBLE PRECISION :: one, zero
    PARAMETER (one=1.d+0, zero=0.d+0)

    !------------------------------ BODY -----------------------------------

    ii = 0
    DO  i = 1, p
      ii = ii+i
      wi = w(i)
      IF (wi == zero) CYCLE
      ri = rmat(ii)
      IF (ri /= zero) GO TO 20
      ij = ii
    !           *** SWAP W AND ROW I OF RMAT ***
      DO  j = i, p
        t = rmat(ij)
        rmat(ij) = w(j)
        w(j) = t
        ij = ij+j
      END DO
      t = qtr(i)
      qtr(i) = y
      y = t
      CYCLE
      20      ip1 = i+1
      ij = ii+i
      IF ( ABS(wi) <=  ABS(ri)) GO TO 40
      rw = ri/wi
      t =  SQRT(one+rw**2)
      IF (rw > zero) t = -t
      v = rw-t
      u1 = one/t
      u2 = one/(t*v)
      rmat(ii) = wi*t
      t = y+v*qtr(i)
      qtr(i) = qtr(i)+t*u1
      y = y+t*u2
      IF (ip1 > p) CYCLE
      DO  j = ip1, p
        t = w(j)+v*rmat(ij)
        rmat(ij) = rmat(ij)+t*u1
        w(j) = w(j)+t*u2
        ij = ij+j
      END DO
      CYCLE

    !        *** AT THIS POINT WE MUST HAVE ABS(WI) .LE. ABS(RI)...

      40      wr = wi/ri
      t = - SQRT(one+wr**2)
      v = wr/(one-t)
      u1 = one/t-one
      u2 = v*u1
      rmat(ii) = ri*t
      t = qtr(i)+v*y
      qtr(i) = qtr(i)+t*u1
      y = y+t*u2
      IF (ip1 > p) CYCLE
      DO  j = ip1, p
        t = rmat(ij)+v*w(j)
        rmat(ij) = rmat(ij)+t*u1
        w(j) = w(j)+t*u2
        ij = ij+j
      END DO
    END DO
    ! 999  RETURN
    RETURN
    END SUBROUTINE dq7adr
    DOUBLE PRECISION FUNCTION drldst(p, d, x, x0)

    !  ***  COMPUTE AND RETURN RELATIVE DIFFERENCE BETWEEN X AND X0  ***
    !  ***  NL2SOL VERSION 2.2  ***

    INTEGER :: p
    DOUBLE PRECISION :: d(p), x(p), x0(p)

    INTEGER :: i
    DOUBLE PRECISION :: emax, t, xmax, zero
    PARAMETER (zero=0.d+0)

    !  ***  BODY  ***

    emax = zero
    xmax = zero
    DO  i = 1, p
      t =  ABS(d(i) * (x(i) - x0(i)))
      IF (emax < t) emax = t
      t = d(i) * ( ABS(x(i)) +  ABS(x0(i)))
      IF (xmax < t) xmax = t
    END DO
    drldst = zero
    IF (xmax > zero) drldst = emax / xmax
    ! 999  RETURN
    RETURN
    !  ***  LAST LINE OF DRLDST FOLLOWS  ***
    END FUNCTION drldst
    SUBROUTINE ds7lup(a, cosmin, p, size, step, u, w, wchmtd, wscale, y)

    !  ***  UPDATE SYMMETRIC  A  SO THAT  A * STEP = Y  ***
    !  ***  (LOWER TRIANGLE OF  A  STORED ROWWISE       ***

    !  ***  PARAMETER DECLARATIONS  ***

    INTEGER :: p
    DOUBLE PRECISION :: a(1), cosmin, size, step(p), u(p), w(p),  &
        wchmtd(p), wscale, y(p)
    !     DIMENSION A(P*(P+1)/2)

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: i, j, k
    DOUBLE PRECISION :: denmin, sdotwm, t, ui, wi

    !     ***  CONSTANTS  ***
    DOUBLE PRECISION :: half, one, zero

    !  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***

    DOUBLE PRECISION :: dd7tpr, dv2nrm
    EXTERNAL dd7tpr, ds7lvm, dv2nrm

    PARAMETER (half=0.5D+0, one=1.d+0, zero=0.d+0)

    !-----------------------------------------------------------------------

    sdotwm = dd7tpr(p, step, wchmtd)
    denmin = cosmin * dv2nrm(p,step) * dv2nrm(p,wchmtd)
    wscale = one
    IF (denmin /= zero) wscale =   MIN(one,  ABS(sdotwm/denmin))
    t = zero
    IF (sdotwm /= zero) t = wscale / sdotwm
    DO  i = 1, p
      w(i) = t * wchmtd(i)
    END DO
    CALL ds7lvm(p, u, a, step)
    t = half * (size * dd7tpr(p, step, u)  -  dd7tpr(p, step, y))
    DO  i = 1, p
      u(i) = t*w(i) + y(i) - size*u(i)
    END DO

    !  ***  SET  A = A + U*(W**T) + W*(U**T)  ***

    k = 1
    DO  i = 1, p
      ui = u(i)
      wi = w(i)
      DO  j = 1, i
        a(k) = size*a(k) + ui*w(j) + wi*u(j)
        k = k + 1
      END DO
    END DO

    ! 999  RETURN
    RETURN
    !  ***  LAST LINE OF DS7LUP FOLLOWS  ***
    END SUBROUTINE ds7lup
    SUBROUTINE ds7lvm(p, y, s, x)

    !  ***  SET  Y = S * X,  S = P X P SYMMETRIC MATRIX.  ***
    !  ***  LOWER TRIANGLE OF  S  STORED ROWWISE.         ***

    !  ***  PARAMETER DECLARATIONS  ***

    INTEGER :: p
    DOUBLE PRECISION :: s(1), x(p), y(p)
    !     DIMENSION S(P*(P+1)/2)

    !  ***  LOCAL VARIABLES  ***

    INTEGER :: i, im1, j, k
    DOUBLE PRECISION :: xi


    !  ***  EXTERNAL FUNCTION  ***

    DOUBLE PRECISION :: dd7tpr
    EXTERNAL dd7tpr

    !-----------------------------------------------------------------------

    j = 1
    DO  i = 1, p
      y(i) = dd7tpr(i, s(j), x)
      j = j + i
    END DO

    IF (p <= 1) GO TO 999
    j = 1
    DO  i = 2, p
      xi = x(i)
      im1 = i - 1
      j = j + 1
      DO  k = 1, im1
        y(k) = y(k) + s(j)*xi
        j = j + 1
      END DO
    END DO

    999  RETURN
    !  ***  LAST LINE OF DS7LVM FOLLOWS  ***
    END SUBROUTINE ds7lvm
    SUBROUTINE dv2axy(p, w, a, x, y)

    !  ***  SET W = A*X + Y  --  W, X, Y = P-VECTORS, A = SCALAR  ***

    INTEGER :: p
    DOUBLE PRECISION :: a, w(p), x(p), y(p)

    INTEGER :: i

    DO  i = 1, p
      w(i) = a*x(i) + y(i)
    END DO
    RETURN
    END SUBROUTINE dv2axy
    DOUBLE PRECISION FUNCTION dv2nrm(p, x)

    !  ***  RETURN THE 2-NORM OF THE P-VECTOR X, TAKING  ***
    !  ***  CARE TO AVOID THE MOST LIKELY UNDERFLOWS.    ***

    INTEGER :: p
    DOUBLE PRECISION :: x(p)

    INTEGER :: i, j
    DOUBLE PRECISION :: one, r, scale, sqteta, t, xi, zero
    DOUBLE PRECISION :: dr7mdc
    EXTERNAL dr7mdc

    PARAMETER (one=1.d+0, zero=0.d+0)
    SAVE sqteta
    DATA sqteta/0.d+0/

    IF (p > 0) GO TO 10
    dv2nrm = zero
    GO TO 999
    10   DO  i = 1, p
      IF (x(i) /= zero) GO TO 30
    END DO
    dv2nrm = zero
    GO TO 999

    30   scale =  ABS(x(i))
    IF (i < p) GO TO 40
    dv2nrm = scale
    GO TO 999
    40   t = one
    IF (sqteta == zero) sqteta = dr7mdc(2)

    !     ***  SQTETA IS (SLIGHTLY LARGER THAN) THE SQUARE ROOT OF THE
    !     ***  SMALLEST POSITIVE FLOATING POINT NUMBER ON THE MACHINE.
    !     ***  THE TESTS INVOLVING SQTETA ARE DONE TO PREVENT UNDERFLOWS.

    j = i + 1
    DO  i = j, p
      xi =  ABS(x(i))
      IF (xi > scale) GO TO 50
      r = xi / scale
      IF (r > sqteta) t = t + r*r
      CYCLE
      50           r = scale / xi
      IF (r <= sqteta) r = zero
      t = one  +  t * r*r
      scale = xi
    END DO

    dv2nrm = scale *  SQRT(t)
    999  RETURN
    !  ***  LAST LINE OF DV2NRM FOLLOWS  ***
    END FUNCTION dv2nrm
    SUBROUTINE dv7cpy(p, y, x)

    !  ***  SET Y = X, WHERE X AND Y ARE P-VECTORS  ***

    INTEGER :: p
    DOUBLE PRECISION :: x(p), y(p)

    INTEGER :: i

    DO  i = 1, p
      y(i) = x(i)
    END DO
    RETURN
    END SUBROUTINE dv7cpy
    SUBROUTINE dv7dfl(alg, lv, v)

    !  ***  SUPPLY ***SOL (VERSION 2.3) DEFAULT VALUES TO V  ***

    !  ***  ALG = 1 MEANS REGRESSION CONSTANTS.
    !  ***  ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS.

    INTEGER :: alg, lv
    DOUBLE PRECISION :: v(lv)

    DOUBLE PRECISION :: dr7mdc
    EXTERNAL dr7mdc
    ! DR7MDC... RETURNS MACHINE-DEPENDENT CONSTANTS

    DOUBLE PRECISION :: machep, mepcrt, one, sqteps, three

    !  ***  SUBSCRIPTS FOR V  ***

    INTEGER :: afctol, bias, cosmin, decfac, delta0, dfac, dinit, dltfdc,  &
        dltfdj, dtinit, d0init, epslon, eta0, fuzz, huberc,  &
        incfac, lmax0, lmaxs, phmnfc, phmxfc, rdfcmn, rdfcmx,  &
        rfctol, rlimit, rsptol, sctol, sigmin, tuner1, tuner2,  &
        tuner3, tuner4, tuner5, xctol, xftol

    PARAMETER (one=1.d+0, three=3.d+0)

    !  ***  V SUBSCRIPT VALUES  ***

    PARAMETER (afctol=31, bias=43, cosmin=47, decfac=22, delta0=44,  &
        dfac=41, dinit=38, dltfdc=42, dltfdj=43, dtinit=39,  &
        d0init=40, epslon=19, eta0=42, fuzz=45, huberc=48,  &
        incfac=23, lmax0=35, lmaxs=36, phmnfc=20, phmxfc=21,  &
        rdfcmn=24, rdfcmx=25, rfctol=32, rlimit=46, rsptol=49,  &
        sctol=37, sigmin=50, tuner1=26, tuner2=27, tuner3=28,  &
        tuner4=29, tuner5=30, xctol=33, xftol=34)

    !-------------------------------  BODY  --------------------------------

    machep = dr7mdc(3)
    v(afctol) = 1.d-20
    IF (machep > 1.d-10) v(afctol) = machep**2
    v(decfac) = 0.5D+0
    sqteps = dr7mdc(4)
    v(dfac) = 0.6D+0
    v(dtinit) = 1.d-6
    mepcrt = machep ** (one/three)
    v(d0init) = 1.d+0
    v(epslon) = 0.1D+0
    v(incfac) = 2.d+0
    v(lmax0) = 1.d+0
    v(lmaxs) = 1.d+0
    v(phmnfc) = -0.1D+0
    v(phmxfc) = 0.1D+0
    v(rdfcmn) = 0.1D+0
    v(rdfcmx) = 4.d+0
    v(rfctol) =   MAX(1.d-10, mepcrt**2)
    v(sctol) = v(rfctol)
    v(tuner1) = 0.1D+0
    v(tuner2) = 1.d-4
    v(tuner3) = 0.75D+0
    v(tuner4) = 0.5D+0
    v(tuner5) = 0.75D+0
    v(xctol) = sqteps
    v(xftol) = 1.d+2 * machep

    IF (alg >= 2) GO TO 10

    !  ***  REGRESSION  VALUES

    v(cosmin) =   MAX(1.d-6, 1.d+2 * machep)
    v(dinit) = 0.d+0
    v(delta0) = sqteps
    v(dltfdc) = mepcrt
    v(dltfdj) = sqteps
    v(fuzz) = 1.5D+0
    v(huberc) = 0.7D+0
    v(rlimit) = dr7mdc(5)
    v(rsptol) = 1.d-3
    v(sigmin) = 1.d-4
    GO TO 999

    !  ***  GENERAL OPTIMIZATION VALUES

    10   v(bias) = 0.8D+0
    v(dinit) = -1.0D+0
    v(eta0) = 1.0D+3 * machep

    999  RETURN
    !  ***  LAST LINE OF DV7DFL FOLLOWS  ***
    END SUBROUTINE dv7dfl
    SUBROUTINE dv7scl(n, x, a, y)

    !  ***  SET X(I) = A*Y(I), I = 1(1)N  ***

    INTEGER :: n
    DOUBLE PRECISION :: a, x(n), y(n)

    INTEGER :: i

    DO  i = 1, n
      x(i) = a * y(i)
    END DO
    !999    RETURN
    RETURN
    !  ***  LAST LINE OF DV7SCL FOLLOWS  ***
    END SUBROUTINE dv7scl
    SUBROUTINE dv7scp(p, y, s)

    !  ***  SET P-VECTOR Y TO SCALAR S  ***

    INTEGER :: p
    DOUBLE PRECISION :: s, y(p)

    INTEGER :: i

    DO  i = 1, p
      y(i) = s
    END DO
    RETURN
    END SUBROUTINE dv7scp
    DOUBLE PRECISION FUNCTION dvsum(n, x)
    INTEGER :: n
    DOUBLE PRECISION :: x(n)
    INTEGER :: i

    dvsum = 0.d+0
    DO  i = 1, n
      dvsum = dvsum + x(i)
    END DO
    END FUNCTION dvsum
    LOGICAL FUNCTION stopx(idummy)
    !     *****PARAMETERS...
    INTEGER :: idummy

    !     ..................................................................

    !     *****PURPOSE...
    !     THIS FUNCTION MAY SERVE AS THE STOPX (ASYNCHRONOUS INTERRUPTION)
    !     FUNCTION FOR THE NL2SOL (NONLINEAR LEAST-SQUARES) PACKAGE AT
    !     THOSE INSTALLATIONS WHICH DO NOT WISH TO IMPLEMENT A
    !     DYNAMIC STOPX.

    !     *****ALGORITHM NOTES...
    !     AT INSTALLATIONS WHERE THE NL2SOL SYSTEM IS USED
    !     INTERACTIVELY, THIS DUMMY STOPX SHOULD BE REPLACED BY A
    !     FUNCTION THAT RETURNS .TRUE. IF AND ONLY IF THE INTERRUPT
    !     (BREAK) KEY HAS BEEN PRESSED SINCE THE LAST CALL ON STOPX.

    !     ..................................................................
    ! DSB NOTE: This stopx thing is a real blast from the past
    ! We can investigate what to do about this at a later time
    ! The following was added just to eliminate a compiler warning
    idummy = 0
    stopx = .false.
    RETURN
    END FUNCTION stopx
