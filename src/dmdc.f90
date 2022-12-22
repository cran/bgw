


DOUBLE PRECISION FUNCTION dr7mdc(k)

! Code converted using TO_F90 by Alan Miller
! Date: 2021-07-13  Time: 16:54:41

!  ***  RETURN MACHINE DEPENDENT CONSTANTS USED BY NL2SOL  ***

! +++  COMMENTS BELOW CONTAIN DATA STATEMENTS FOR VARIOUS MACHINES.  +++
! +++  TO CONVERT TO ANOTHER MACHINE, PLACE A C IN COLUMN 1 OF THE   +++
! +++  DATA STATEMENT LINE(S) THAT CORRESPOND TO THE CURRENT MACHINE +++
! +++  AND REMOVE THE C FROM COLUMN 1 OF THE DATA STATEMENT LINE(S)  +++
! +++  THAT CORRESPOND TO THE NEW MACHINE.                           +++


INTEGER, INTENT(IN)                  :: k


!  ***  THE CONSTANT RETURNED DEPENDS ON K...

!  ***        K = 1... SMALLEST POS. ETA SUCH THAT -ETA EXISTS.
!  ***        K = 2... SQUARE ROOT OF ETA.
!  ***        K = 3... UNIT ROUNDOFF = SMALLEST POS. NO. MACHEP SUCH
!  ***                 THAT 1 + MACHEP .GT. 1 .AND. 1 - MACHEP .LT. 1.
!  ***        K = 4... SQUARE ROOT OF MACHEP.
!  ***        K = 5... SQUARE ROOT OF BIG (SEE K = 6).
!  ***        K = 6... LARGEST MACHINE NO. BIG SUCH THAT -BIG EXISTS.

! DSB NOTE: In this version of dmdc, we have replaced hard-coded constants
! with calls to internal functions.
DOUBLE PRECISION :: machep
!  DOUBLE PRECISION :: big, eta, machep
! INTEGER :: bigi(2), etai(2), machei(2)
! EQUIVALENCE (big,bigi(1)), (eta,etai(1)), (machep,machei(1))
DOUBLE PRECISION, PARAMETER :: zero=0.d+0
DOUBLE PRECISION, PARAMETER :: one=1.d+0
INTRINSIC epsilon, huge, tiny
! CHARACTER (LEN=132) :: output_string

!  +++ IEEE ARITHMETIC MACHINES IN WHICH THE MOST SIGNIFICANT BYTE
!  +++ IS STORED FIRST, SUCH AS THE AT&T 3B SERIES AND MACHINES
!  +++ BASED ON SPARC, MIPS, AND MOTOROLA 68XXX PROCESSORS.

!    DATA BIGI(1),BIGI(2)     / 2146435071,         -1 /
!      DATA ETAI(1),ETAI(2)     /    1048576,          0 /
!      DATA MACHEI(1),MACHEI(2) / 1017118720,          0 /

!  +++ IEEE ARITHMETIC MACHINES IN WHICH THE LEAST SIGNIFICANT BYTE
!  +++ IS STORED FIRST, SUCH AS MACHINES BASED ON INTEL PROCESSORS,
!  +++ E.G. PERSONAL COMPUTERS WITH AN INTEL 80X87.

!!!DATA bigi(1),bigi(2)     / -1, 2146435071 /
!!! etai(1),etai(2)     /  0,    1048576 /
!!!DATA machei(1),machei(2) /  0, 1017118720 /

!  +++  IBM, AMDAHL, OR XEROX MAINFRAME  +++

!      DATA BIGI(1),BIGI(2)/2147483647, -1/
!      DATA ETAI(1),ETAI(2)/1048576, 0/
!      DATA MACHEI(1),MACHEI(2)/873463808,0/

!  +++  VAX  +++

!      DATA BIGI(1),BIGI(2)     / -32769, -1 /
!      DATA ETAI(1),ETAI(2)     /    128,  0 /
!      DATA MACHEI(1),MACHEI(2) /   9344,  0 /

!  +++  CRAY  +++

!      DATA BIGI(1)/6917247552664371199/
!      DATA BIGI(2)/128891879815246481/
!      DATA ETAI(1)/2332160919536140288/
!      DATA ETAI(2)/0/
!      DATA MACHEI(1)/4585931058058362880/
!      DATA MACHEI(2)/0/

!  +++  PORT LIBRARY -- REQUIRES MORE THAN JUST A DATA STATEMENT, +++
!  +++                  BUT HAS CONSTANTS FOR MANY MORE MACHINES. +++

!  To get the current D1MACH, which has constants for many more
!  machines, ask netlib@research.att.com to
!                    send d1mach from cor
!  For machines with rounded arithmetic (e.g., IEEE or VAX arithmetic),
!  use MACHEP = 0.5D0 * D1MACH(4) below.

!      DOUBLE PRECISION D1MACH
!      EXTERNAL D1MACH
!      DATA BIG/0.D+0/, ETA/0.D+0/, MACHEP/0.D+0/, ZERO/0.D+0/
!      IF (BIG .GT. ZERO) GO TO 1
!         BIG = D1MACH(2)
!         ETA = D1MACH(1)
!         MACHEP = D1MACH(4)
!1     CONTINUE

!  +++ END OF PORT +++

!-------------------------------  BODY  --------------------------------
machep = epsilon(1.d0)
! output_string = ' epsilon = '
! CALL dblepr1(output_string,11,machep)
IF (machep <= zero) THEN
  WRITE(*,*) 'Edit DR7MDC to activate the appropriate statements'
  STOP 987
END IF
SELECT CASE ( k )
  CASE (    1)
    GO TO 10
  CASE (    2)
    GO TO  20
  CASE (    3)
    GO TO  30
  CASE (    4)
    GO TO  40
  CASE (    5)
    GO TO  50
  CASE (    6)
    GO TO  60
END SELECT

! 10   dr7mdc = eta
10   dr7mdc = tiny(1.d0)
GO TO 999

! 20   dr7mdc = SQRT(256.d+0*eta)/16.d+0
20   dr7mdc = SQRT(256.d+0*tiny(1.d0))/16.d+0
GO TO 999

! 30   dr7mdc = machep
30   dr7mdc = epsilon(1.d0)
! output_string = ' dr7mdc = '
! CALL dblepr1(output_string,10,dr7mdc)
GO TO 999

! 40   dr7mdc = SQRT(machep)
40   dr7mdc = SQRT(epsilon(1.d0))
GO TO 999

! 50   dr7mdc = SQRT(big/256.d+0)*16.d+0
50   dr7mdc = SQRT(huge(1.d0)/256.d+0)*16.d+0
GO TO 999

! 60   dr7mdc = big
60   dr7mdc = huge(1.d0)

999  RETURN
!  ***  LAST LINE OF DR7MDC FOLLOWS  ***
END FUNCTION dr7mdc

INTEGER FUNCTION i7mdcn(k)


INTEGER, INTENT(IN OUT)                  :: k


!  ***  RETURN INTEGER MACHINE-DEPENDENT CONSTANTS  ***

!     ***  K = 1 MEANS RETURN STANDARD OUTPUT UNIT NUMBER.   ***
!     ***  K = 2 MEANS RETURN ALTERNATE OUTPUT UNIT NUMBER.  ***
!     ***  K = 3 MEANS RETURN  INPUT UNIT NUMBER.            ***
!          (NOTE -- K = 2, 3 ARE USED ONLY BY TEST PROGRAMS.)

!  +++  PORT VERSION FOLLOWS...
!      INTEGER I1MACH
!      EXTERNAL I1MACH
!      INTEGER MDPERM(3)
!      DATA MDPERM(1)/2/, MDPERM(2)/4/, MDPERM(3)/1/
!      I7MDCN = I1MACH(MDPERM(K))
!  +++  END OF PORT VERSION  +++

!  +++  NON-PORT VERSION FOLLOWS...
INTEGER :: mdcon(3)
DATA mdcon(1)/6/, mdcon(2)/8/, mdcon(3)/5/

i7mdcn = mdcon(k)
!  +++  END OF NON-PORT VERSION  +++

! 999  RETURN
RETURN
!  ***  LAST LINE OF I7MDCN FOLLOWS  ***
END FUNCTION i7mdcn
