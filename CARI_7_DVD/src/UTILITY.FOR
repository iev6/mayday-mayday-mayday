      FUNCTION CHAR2REAL(CHAR_IN,D)
      IMPLICIT NONE
! 
! CAHR2REAL CONVERTS SHORT STRINGS TO REAL NUMBERS
!
      INTEGER(4) :: D
      REAL(8) :: CHAR2REAL,ZED
      CHARACTER(D)::CHAR_IN
          
      IF (CHAR_IN.EQ.'') THEN 
         CHAR2REAL=0.0
         RETURN
      ENDIF       

      OPEN (UNIT=99,STATUS='SCRATCH') 
      WRITE(99,*) CHAR_IN
      REWIND 99
      READ(99,199) ZED
      CLOSE(99) 
      CHAR2REAL=ZED
199   FORMAT (G16.7E2)
      END FUNCTION CHAR2REAL
!                                                                      7
!----6-----------------------------------------------------------------2
! SUB TO CONVERT  SMALL INTEGERS TO SHORT CHARACTER STRINGS 
! FOR PRINTING DATE FORMATS
!
      SUBROUTINE INT2CHAR(I,CC,L)
        INTEGER, INTENT(IN)::I
        CHARACTER(L), INTENT(OUT)::CC
        
        IF (I.LT.10 .AND. L.EQ.1) THEN
           IF (I.EQ.0) CC='0'
         IF (I.EQ.1) CC='1'
          IF (I.EQ.2) CC='2'
          IF (I.EQ.3) CC='3'
          IF (I.EQ.4) CC='4'
           IF (I.EQ.5) CC='5'
          IF (I.EQ.6) CC='6'
          IF (I.EQ.7) CC='7'
          IF (I.EQ.8) CC='8'
         IF (I.EQ.9) CC='9'
      ELSEIF (I.LT.10 .AND. L.EQ.2) THEN
           IF (I.EQ.0) CC='00'
         IF (I.EQ.1) CC='01'
          IF (I.EQ.2) CC='02'
          IF (I.EQ.3) CC='03'
          IF (I.EQ.4) CC='04'
           IF (I.EQ.5) CC='05'
          IF (I.EQ.6) CC='06'
          IF (I.EQ.7) CC='07'
          IF (I.EQ.8) CC='08'
         IF (I.EQ.9) CC='09'
      ELSEIF (I.GT.9 .AND. I.LT.100 .AND. L.EQ.2) THEN
         OPEN(UNIT=99,STATUS='SCRATCH')
          WRITE(99,20010) I
          REWIND(99)
          READ(99,20011) CC
          CLOSE(99)       
      ELSE
         CALL EPITATH(' NUMBER IN INT2CHAR TOO LARGE ', 29)      
      ENDIF
20010 FORMAT(I2)
20011 FORMAT(A2)    
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2

      FUNCTION CHAR2INT(CHAR_IN,D)
      IMPLICIT NONE
! 
! CHAR2INT CONVERTS SHORT STRINGS TO INTEGERS
!     
      INTEGER(4)::CHAR2INT,I,D
      CHARACTER(D)::CHAR_IN

         IF (TRIM(CHAR_IN).EQ.'') THEN 
            CHAR2INT=0
            RETURN
         ENDIF
            DO I=1,D
               IF (CHAR_IN(I:I).EQ.'-') THEN
                  CHAR_IN=CHAR_IN(I:D)
                  EXIT
               ENDIF
            ENDDO     
         OPEN (UNIT=99,STATUS='SCRATCH') 
         WRITE(99,*) CHAR_IN
         REWIND 99
         READ(99,299) I
         CLOSE(99) 
         CHAR2INT=I
299      FORMAT (I8)   
      END FUNCTION CHAR2INT
!----------------------------------------------------------------------
! Linear interpolation function	  
      FUNCTION linterp(x1,x2,y1,y2,x)
      IMPLICIT NONE
      REAL(8)::m,x1,x2,y1,y2,x,linterp,y

!	WRITE(130,*)'x1,x2,y1,y2,x',x1,x2,y1,y2,x  
!	Write(1,*)'x1,x2,y1,y2,x',x1,x2,y1,y2,x  
      m=(y2-y1)/(x2-x1)
      y=m*(x-x2)+y2
!	WRITE(130,*)'y,m',y,m  
!	Write(1,*)'y,m',y,m  
      linterp = y
	  
      END FUNCTION linterp 	  
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE TDSPLINE(N,X,Y,YP1,YPN,XO,YO)
      IMPLICIT NONE
!
!       THIS SUBROUTINE WILL CALCULATE THE CUBIC-SPLINE INTEPOLATED VALUE OF
!	    GIVEN ANY FUNCTION FOR SINGLE VARIABLE OF ANY DIMENSION
!      
!       INPUTS/OUTPUT:
!       N : (INPUT) [SCALAR] DIMENSION OF INPUT ARRAY
!       X : (INPUT) [VECTOR OF LENGTH 'N'], INDEPENDENT VARIABLE OF DIMENSION 'N'
!       Y : (INPUT) [VECTOR OF LENGTH 'N']' THE FUNCTION OF X (AN ARRAY OF DIMENSION 'N')
!       YP1,YPN : (INPUT) [SCALAR] =1.E+30 :: THE ROUTINE IS SIGNALED TO SET
!            THE CORRESPONDING BOUNDARY CONDITION
!            FOR A NATURAL SPLINE, WITH ZERO
!            SECOND DERIVATIVE ON THAT BOUNDARY
!       XO : (INPUT) [SCALAR] 'X' VALUE AT WHICH YOU WANT TO SPLINE
!       YO : (OUTPUT) [SCALAR] SPLINED VALUE
!      
!       POSTED BY TANMOY DAS 
!       TO
!       HTTP://ENCODEDCODE.BLOGSPOT.COM/2008/10/SPLINE-INTERPOLATION-FUNCTION-FOR.HTML
!       FEB 13, 2008.
!
!       ADAPTED FOR USE IN CARI BY KC  5 NOV 2013
      INTEGER(4)::N,NO,I,IO,K,KLO,KHI
      REAL(8),DIMENSION(1:N)::X,Y,Y2,U
      REAL(8)::XO,YO
      REAL(8)::YP1,YPN,A,B,C,H,SIG,P,QN,UN
      CHARACTER*3::DIAGNOSE

      DIAGNOSE='YES'
      DIAGNOSE='NO!'
!     KC MODIFIED FOR X0 OUTISDE BOUNDARIES OF X(N)
      IF (DIAGNOSE.eq.'YES') THEN
         OPEN(43,file='DIAGNOSE/TDspline.chk',status='unknown')
      ENDIF    
      IF (DIAGNOSE.eq.'YES') THEN
         DO I=1,N
            WRITE(43,*) X(I),Y(I)
         ENDDO
         WRITE(43,*) 'Test value is ',XO
      ENDIF
      IF ((XO.LT.X(1)).AND.(XO.LT.X(N))) THEN
      !   YO=Y AT XMIN
         IF (X(1).LT.X(N)) THEN
            YO=Y(1)
         ELSE
            YO=Y(N)
         ENDIF
         IF (DIAGNOSE.eq.'YES') CLOSE(43)
         RETURN
      ENDIF
      IF ((XO.GT.X(1)).AND.(XO.GT.X(N))) THEN
      !   YO=Y AT XMAX
         IF (X(1).LT.X(N)) THEN
            YO=Y(N)
         ELSE
            YO=Y(1)
         ENDIF
         IF (DIAGNOSE.eq.'YES') CLOSE(43)
         RETURN
      ENDIF
! ELSE, IF X IS IN RANGE, USE THE SPLINE 
      IF (YP1.GT.0.99E30) THEN
         Y2(1)=0.
         U(1)=0.
      ELSE
         Y2(1)=-0.5
         U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.
         Y2(I)=(SIG-1.)/P
         U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     &        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      END DO
      IF (YPN.GT.0.99E30) THEN
         QN=0.
         UN=0.
      ELSE
         QN=0.5
         UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO K=N-1,1,-1
         Y2(K)=Y2(K)*Y2(K+1)+U(K)
      END DO
      KLO=1
      KHI=N
10010 IF (KHI-KLO.GT.1) THEN
         K=(KHI+KLO)/2.
         IF (X(K).GT.XO) THEN
            KHI=K
         ELSE
            KLO=K
         ENDIF
         GOTO 10010
      ENDIF
      H=X(KHI)-X(KLO)
      A=(X(KHI)-XO)/H
      B=(XO-X(KLO))/H
      C=(A*(A*A-1.)*Y2(KLO)+B*(B*B-1.)*Y2(KHI))*H*H/6.
      YO=A*Y(KLO)+B*Y(KHI)+C
      IF (DIAGNOSE.eq.'YES') WRITE(43,*) 'TDspline answer is ',YO 
      IF (DIAGNOSE.eq.'YES') CLOSE(43) 
      END SUBROUTINE TDSPLINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE KCSPLINE(N,X,Y,YP1,YPN,XO,YO)
      IMPLICIT NONE
!
!       THIS SUBROUTINE WILL CALCULATE THE LINEAR INTERPOLATED VALUE OF
!	    GIVEN ANY FUNCTION FOR SINGLE VARIABLE OF ANY DIMENSION
!
!       IT IS USED WHEREVER TDSPLINE IS TOO UNSTABLE FOR THE DATA  
!      
!       INPUTS/OUTPUT:
!       N : (INPUT) [SCALAR] DIMENSION OF INPUT ARRAY
!       X : (INPUT) [VECTOR OF LENGTH 'N'], INDEPENDENT VARIABLE OF DIMENSION 'N'
!       Y : (INPUT) [VECTOR OF LENGTH 'N']' THE FUNCTION OF X (AN ARRAY OF DIMENSION 'N')
!       YP1,YPN : (INPUT) [SCALAR] =1.E+30 :: NOT USED 
!       XO : (INPUT) [SCALAR] 'X' VALUE AT WHICH YOU WANT 'Y'
!       YO : (OUTPUT) [SCALAR] INTERPOLATED VALUE
!
!       BUILT FOR USE IN CARI BY KC  31 DEC 2013
      INTEGER(4)::N,NO,I,IO
      REAL(8),DIMENSION(1:N)::X,Y
      REAL(8)::XO,YO,YP1,YPN,XLO,XHI,YLO,YHI,LINTERP
      logical::A,B,C,D
      CHARACTER*3::DIAGNOSE

      DIAGNOSE='YES'
      DIAGNOSE='NO!'
!     KC MODIFIED FOR X0 OUTISDE BOUNDARIES OF X(N)
      IF (DIAGNOSE.eq.'YES') THEN
         OPEN(42,file='DIAGNOSE/kcspline.chk',status='unknown')
      ENDIF    
      IF (DIAGNOSE.eq.'YES') THEN
         DO I=1,N
            WRITE(42,*) X(I),Y(I)
         ENDDO
         WRITE(42,*) 'Test value is ',XO
      ENDIF
      IF ((XO.LT.X(1)).AND.(XO.LT.X(N))) THEN
      !   YO=Y AT XMIN
         IF (X(1).LT.X(N)) THEN
            YO=Y(1)
         ELSE
            YO=Y(N)
         ENDIF
         IF (DIAGNOSE.eq.'YES') CLOSE(42)
         RETURN
      ENDIF
      IF ((XO.GT.X(1)).AND.(XO.GT.X(N))) THEN
      !   YO=Y AT XMAX
         IF (X(1).LT.X(N)) THEN
            YO=Y(N)
         ELSE
            YO=Y(1)
         ENDIF
         IF (DIAGNOSE.eq.'YES') CLOSE(42)
         RETURN
      ENDIF
! ELSE, IF X IS IN RANGE, USE LINEAR INTERPOLATION
      DO I=1,N-1
         A=(XO.GE.X(I))
         B=(XO.LT.X(I+1))
         C=(XO.GE.X(I+1))
         D=(XO.LT.X(I)) 
         IF ((A.AND.B).OR.(C.AND.D)) THEN
            XLO=X(I)
            XHI=X(I+1)
            YLO=Y(I)
            YHI=Y(I+1) 
            YO=LINTERP(XLO,XHI,YLO,YHI,XO)
            IF (Y(I).EQ.0.0) THEN 
               Y(I)=1E-30
            ENDIF    
            IF (Y(I+1).EQ.0.0) THEN 
               Y(I+1)=1E-30
            ENDIF    
            EXIT
         ENDIF
      ENDDO 
      ! Corrects for computer rounding errors when interpolating close to 0 
      ! between small positive numbers and zero
      IF ((YO.LT.0.).AND.(YLO.GE.0.).AND.(YHI.GE.0.)) YO=0. 
      IF (DIAGNOSE.eq.'YES') WRITE(42,*) 'KCspline answer is ',YO 
      IF (DIAGNOSE.eq.'YES') CLOSE(42)
296   FORMAT(6ES12.4)
      END SUBROUTINE KCspline
!                                                                      7
!----6-----------------------------------------------------------------2
!      FUNCTION SPLINE_4PT(DELTA, VI, VJ, VK, VL, X)
!
! A CATMUL-ROM CUBIC SPLINE
! A CUBIC SPLINE REQUIRING 4 DATA POINTS, 1 ABOVE X, 2 BRACKETING X, 
!   AND 1 BELOW X. ASSUMES  0<=x<1, EVENLY SPACED POINTS ON X AXIS
!   B AND C ARE BRACKETING X COORDS, VA-VD ARE VALUES OF F(X) AT 
!   POINTS A-D, XN IS NORMALIZED X
!
!      INTEGER :: J
!      REAL :: DELTA, X, XN, VI, VJ, VK, VL, F, A0, A1, A2, A3
!   	  REAL :: SPLINE_4PT
!      CHARACTER(3) :: DIAGNOSE='NO!' 
!
! 1. NORMALIZE X TO A VALUE FROM 0 TO 1 [INCLUSIVE]
!     DELTA = DISTANCE BETWEEN POINTS  
!      XN = MOD(X,DELTA)/DELTA 
!      WRITE(40,*) 'XN=',XN,' X=',X,' DELTA=',DELTA 
	  !
! 2. CALCULATE F(X)
!      A0 = -0.5*VI +1.5*VJ -1.5*VK+0.5* VL
!      A1 = VI - 2.5*VJ + 2.0*VK - 0.5*VL
!      A2 = 0.5*VK - 0.5*VI
!      A3 = VJ
!      F= A0*XN**3+A1*XN**2+A2*XN+A3
!      SPLINE_4PT = F
!      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'SPLINE INPUT', VI,VJ,VK,VL,X
!      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) 'SPLINE INPUT', VI,VJ,VK,VL,X
!      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'COEFFICIENTS', A0,A1,A2,A3
!      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) 'COEFFICIENTS', A0,A1,A2,A3
!      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'SPLINE STEP AND OUTPUT'       &
!     &   ,DELTA, F
!      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) 'SPLINE STEP AND OUTPUT'        &
!     &   ,DELTA, F
!      END FUNCTION SPLINE_4PT
!                                                                      7
!----6-----------------------------------------------------------------2
! Convert year, month, day, and hour to fractional year 
      SUBROUTINE YMDH2T(Y,M,D,H,T)
      IMPLICIT NONE
      REAL(8),INTENT(OUT)::T
      INTEGER(4),INTENT(IN)::Y,M,D,H
      INTEGER(4)::DIY
      REAL(8)::DUMMY,DSF,LEAPYR

      LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE

      IF (M.EQ.0) YEARLYAVE=.TRUE. 
      IF (D.EQ.0) MONTHLYAVE=.TRUE.
      IF (H.EQ.0) DAILYAVE=.TRUE.

      DUMMY=REAL(Y,kind=8)
      DUMMY=DUMMY+DSF(Y,M,D,H)/(365.+LEAPYR(Y))
      T=DUMMY

      END SUBROUTINE YMDH2T
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION DSF(Y,M,D,H)
      !Days so far
      IMPLICIT NONE
      REAL(8)::DSF,DUMMY,DOM,DIN,LEAPYR
      INTEGER(4)::Y,M,D,H

      IF (D.EQ.0) THEN
          DIN = 15 !THIS IS A MONTHLY AVERAGE
      ELSE
          DIN = D
      ENDIF
         DUMMY=DIN+(H/24.)

      IF (M .EQ. 1) DOM=  0.0
      IF (M .EQ. 2) DOM= 31.0
      IF (M .EQ. 3) DOM= 59.0+LEAPYR(Y)
      IF (M .EQ. 4) DOM= 90.0+LEAPYR(Y)
      IF (M .EQ. 5) DOM=120.0+LEAPYR(Y)
      IF (M .EQ. 6) DOM=151.0+LEAPYR(Y)
      IF (M .EQ. 7) DOM=181.0+LEAPYR(Y)
      IF (M .EQ. 8) DOM=212.0+LEAPYR(Y)
      IF (M .EQ. 9) DOM=243.0+LEAPYR(Y)
      IF (M .EQ.10) DOM=273.0+LEAPYR(Y)
      IF (M .EQ.11) DOM=304.0+LEAPYR(Y)
      IF (M .EQ.12) DOM=334.0+LEAPYR(Y)

      DUMMY=DUMMY+DOM
      DSF=DUMMY
      END FUNCTION DSF
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION LEAPYR(Y)
      IMPLICIT NONE
      INTEGER(4)::Y
      REAL(8)::LEAPYR,L

!     GREGORIAN CALENDAR SCHEME

      IF (MOD(Y,400).EQ.0) THEN 
         L=1.0  
      ELSEIF (MOD(Y,100).EQ.0) THEN 
         L=0.0! not_leap_year
      ELSEIF (MOD(Y,4).EQ.0) THEN 
         L=1.0
      ELSE 
         L=0.0
      ENDIF
      LEAPYR=L
      END FUNCTION LEAPYR
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE CLS
! MIMICS THE CLEAR SCREEN COMMAND CLS IN QBASIC
      DO I = 1 , 50 ! 50 FOR LARGE DOS BOXES
         PRINT*,' ' 
      ENDDO
!      CALL SYSTEM('CLS')
      END SUBROUTINE CLS
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE OOPS(MESSAGE,L)
!
! Equivalent to the obsolete pause statement, with an added message to 
! the user
! 
      INTEGER::L
      CHARACTER(L)::MESSAGE  
      CHARACTER(1)::GOON
      CHARACTER(12)::VIEWER
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
                   
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      
      WRITE(*,10100) MESSAGE
      IF (DIAGNOSE.EQ.'YES') THEN
        WRITE(40,*) MESSAGE
      ENDIF
      IF (MENUS.EQ.'YES') THEN
        PRINT*,' '
        WRITE(*,10101)'Press a character then press <ENTER> to continue'
        READ*, GOON 
      ENDIF
10100 FORMAT(10X,A60)  
10101 FORMAT(10X,A50)
      END SUBROUTINE 
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE EPITATH(MESSAGE,L)

!
! Error message to user, then kill program 
! 
      INTEGER::L
      CHARACTER(L)::MESSAGE  
      CHARACTER(1)::GOON
      CHARACTER(3)::DIAGNOSE='YES'

      WRITE(*,10102) MESSAGE
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) MESSAGE
      PRINT*,' '
      WRITE(*,*) 'Fatal error, aborting program!'
!      WRITE(*,10103)'Press a character then press <ENTER> to continue'
!      READ*, GOON 
      STOP
!      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'GOODBYE'
10102 FORMAT(10X,A60)  
10103 FORMAT(10X,A50)
      END SUBROUTINE 
!                                                                      7
!----6-----------------------------------------------------------------2
! LOADS Doses and fluences from MCNPX 2.7.0. shower data
! data read has been patched, see details below, KC 30OCT2014 
      SUBROUTINE LOADFTDCCS(OS)
      IMPLICIT NONE

      CHARACTER*72::HEADER2
      CHARACTER*48::HEADER1
      CHARACTER*8::SUBDIR
      CHARACTER*5::OS 
      CHARACTER*12, DIMENSION(29)::FTDCCF
      INTEGER(4)::IZ,IA,IE,IS
      INTEGER(4)::J1,J2,J3,J4
!     FTDCC matrices for fluence, E103, E60, H*(10), and whole body dose
      REAL(8), DIMENSION(27,19,19,38)::PHI,E103,E60,H10,GRAY
!     Variance information matrices for ftdccs 
      REAL(8), DIMENSION(27,19,19,38)::PHIS,E103S,E60S,H10S,GRAYS
      REAL(8)::patch
       
      COMMON /ftdcc/PHI,E103,E60,H10,GRAY
      COMMON /ftdccs/PHIS,E103S,E60S,H10S,GRAYS
     

      DATA FTDCCF/'nushower.dat','H_shower.dat','Heshower.dat',         &
     & 'Lishower.dat','Beshower.dat','B_shower.dat','C_shower.dat',     &
     & 'N_shower.dat','O_shower.dat','F_shower.dat','Neshower.dat',     &
     & 'Nashower.dat','Mgshower.dat','Alshower.dat','Sishower.dat',     &
     & 'P_shower.dat','S_shower.dat','Clshower.dat','Arshower.dat',     &
     & 'K_shower.dat','Cashower.dat','Scshower.dat','Tishower.dat',     &
     & 'V_shower.dat','Crshower.dat','Mnshower.dat','Feshower.dat',     &
     & 'Coshower.dat','Nishower.dat'/

      WRITE(*,*) 'OS is ',OS(1:3)
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
      ! USE DOS/WINDOWS PATHS
      SUBDIR='SHOWERS\'
      ELSE 
      ! USE LINUX/UNIX PATHS
      SUBDIR='SHOWERS/'
      ENDIF


      ! LOAD ALL THE FTDCCS INTO MEMORY, SHARE IN ftdccs COMMON BLOCK
      ! open GCR shower data files
      OPEN(unit=100,file=SUBDIR//FTDCCF(1),status='OLD')
      OPEN(unit=101,file=SUBDIR//FTDCCF(2),status='OLD')
      OPEN(unit=102,file=SUBDIR//FTDCCF(3),status='OLD')
      OPEN(unit=103,file=SUBDIR//FTDCCF(4),status='OLD')
      OPEN(unit=104,file=SUBDIR//FTDCCF(5),status='OLD')
      OPEN(unit=105,file=SUBDIR//FTDCCF(6),status='OLD')
      OPEN(unit=106,file=SUBDIR//FTDCCF(7),status='OLD')
      OPEN(unit=107,file=SUBDIR//FTDCCF(8),status='OLD')
      OPEN(unit=108,file=SUBDIR//FTDCCF(9),status='OLD')
      OPEN(unit=109,file=SUBDIR//FTDCCF(10),status='OLD')
      OPEN(unit=110,file=SUBDIR//FTDCCF(11),status='OLD')
      OPEN(unit=111,file=SUBDIR//FTDCCF(12),status='OLD')
      OPEN(unit=112,file=SUBDIR//FTDCCF(13),status='OLD')
      OPEN(unit=113,file=SUBDIR//FTDCCF(14),status='OLD')
      OPEN(unit=114,file=SUBDIR//FTDCCF(15),status='OLD')
      OPEN(unit=115,file=SUBDIR//FTDCCF(16),status='OLD')
      OPEN(unit=116,file=SUBDIR//FTDCCF(17),status='OLD')
      OPEN(unit=117,file=SUBDIR//FTDCCF(18),status='OLD')
      OPEN(unit=118,file=SUBDIR//FTDCCF(19),status='OLD')
      OPEN(unit=119,file=SUBDIR//FTDCCF(20),status='OLD')
      OPEN(unit=120,file=SUBDIR//FTDCCF(21),status='OLD')
      OPEN(unit=121,file=SUBDIR//FTDCCF(22),status='OLD')
      OPEN(unit=122,file=SUBDIR//FTDCCF(23),status='OLD')
      OPEN(unit=123,file=SUBDIR//FTDCCF(24),status='OLD')
      OPEN(unit=124,file=SUBDIR//FTDCCF(25),status='OLD')
      OPEN(unit=125,file=SUBDIR//FTDCCF(26),status='OLD')
      OPEN(unit=126,file=SUBDIR//FTDCCF(27),status='OLD')
!      OPEN(unit=127,file=SUBDIR//FTDCCF(28),status='OLD')
!      OPEN(unit=128,file=SUBDIR//FTDCCF(29),status='OLD')

      DO IZ=1,27 !IZ=1+PRIMARY NUCLEAR CHARGE
        DO IE=1,19
          READ(UNIT=99+IZ,FMT=10001) HEADER1
          READ(UNIT=99+IZ,FMT=10002) HEADER2
          DO IA=1,19
            DO IS=1,38
              READ(UNIT=99+IZ,FMT=10000) J1,J2,J3,J4,PHI(IZ,IE,IA,IS),  &
     &        PHIS(IZ,IE,IA,IS),E103(IZ,IE,IA,IS),E103S(IZ,IE,IA,IS),   &
     &        E60(IZ,IE,IA,IS),E60S(IZ,IE,IA,IS),H10(IZ,IE,IA,IS),      &
     &        H10S(IZ,IE,IA,IS),GRAY(IZ,IE,IA,IS),GRAYS(IZ,IE,IA,IS)     
            END DO
          END DO
        END DO
        CLOSE(IZ+99)
      END DO

! KC 30OCT2014 Below is my patch to a database writing/processing error I havent found  
! yet and to the discovery that the highest altitude fluences (Surface 101 in the MCNPX 
! files) are missing the primary particle fluence, i.e. they are wrong in that they 
! do not contain the primary incoming fluence, ony the secondary fluence. This is 
! clear because particle fluence for primary HZEs increase with depth, despite having
! no sources within the HZE cascade. Also, the fluence at the top of the atmoshpere
! is orders of magnitude too low (e.g. 8E-22, it should be about 2E-19 based on 
! geometry).  
!
      GO TO 9999
      DO IZ=1,27 !IZ=1+PRIMARY NUCLEAR CHARGE
        DO IE=1,19
          IA=19
          DO IS=1,38
          PHI(IZ,IE,IA,IS) = PHI(IZ,IE,IA,IS)+ 1.9E-19
          PHIS(IZ,IE,IA,IS) = PHI(IZ,IE,IA,IS)*.001
          END DO
        END DO
      END DO
      DO IZ=1,27 !IZ=1+PRIMARY NUCLEAR CHARGE
        DO IE=1,19
         DO IA=19,2,-1
          DO IS=1,38
             E103(IZ,IE,IA,IS)=E103(IZ,IE,IA-1,IS)
             E103S(IZ,IE,IA,IS)=E103S(IZ,IE,IA-1,IS) 
             E60(IZ,IE,IA,IS)=E60(IZ,IE,IA-1,IS)             
             E60S(IZ,IE,IA,IS)=E60S(IZ,IE,IA-1,IS) 
             H10(IZ,IE,IA,IS)=H10(IZ,IE,IA-1,IS)
             H10S(IZ,IE,IA,IS)=H10S(IZ,IE,IA-1,IS)
             GRAY(IZ,IE,IA,IS)=GRAY(IZ,IE,IA-1,IS)
             GRAYS(IZ,IE,IA,IS)=GRAYS(IZ,IE,IA-1,IS)     
          END DO
         END DO
        END DO
      END DO
      DO IZ=1,27 !IZ=1+PRIMARY NUCLEAR CHARGE
        DO IE=1,19
         IA=1
          DO IS=1,38
             IF (PHI(IZ,IE,IA+1,IS).GT.0.0) THEN 
                patch=PHI(IZ,IE,IA,IS)/PHI(IZ,IE,IA+1,IS)
             ELSE
                patch=1. 
             ENDIF
             E103(IZ,IE,IA,IS)=E103(IZ,IE,IA+1,IS)*patch
             E103S(IZ,IE,IA,IS)=E103S(IZ,IE,IA+1,IS)*patch 
             E60(IZ,IE,IA,IS)=E60(IZ,IE,IA+1,IS)*patch             
             E60S(IZ,IE,IA,IS)=E60S(IZ,IE,IA+1,IS)*patch 
             H10(IZ,IE,IA,IS)=H10(IZ,IE,IA+1,IS)*patch
             H10S(IZ,IE,IA,IS)=H10S(IZ,IE,IA+1,IS)*patch
             GRAY(IZ,IE,IA,IS)=GRAY(IZ,IE,IA+1,IS)*patch
             GRAYS(IZ,IE,IA,IS)=GRAYS(IZ,IE,IA+1,IS)*patch     
          END DO
!        END DO
        END DO
      END DO
9999  CONTINUE
! end 30OCT2014 patch
10000 FORMAT (A2,I10,2I4,10ES10.3)
10001 FORMAT (A48)
10002 FORMAT (A72)
10003 FORMAT (2I3,1X,5ES10.3)

      END SUBROUTINE 
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION KM2DEPTH(km)
!     convert km to gpcms 
      IMPLICIT NONE
      REAL(8)::KM2DEPTH,gpcms,km,feet
       
! recall km=kilofeet*0.3048       

      feet = km / 0.0003048
      CALL FT2GPCMS(feet,gpcms)
      KM2DEPTH = gpcms

!      WRITE(40,*) " converted km = ", km, " to depth = ",gpcms
      
      END FUNCTION 
!                                                                      7
!----6-----------------------------------------------------------------2
!______________________________________________________________________
! P(x) calculates momentum in GV from KE in GeV
      FUNCTION P(X,Z)
      IMPLICIT NONE
      REAL(8)::P,X,Z
      REAL(8), DIMENSION(28) :: A
      INTEGER(4), DIMENSION(28)::N
      COMMON /ATOMIC/N,A
!
      P = SQRT(X*X+2.*X*.938272) !GV from GeV/Nucleon
!
      END FUNCTION P
!______________________________________________________________________
! R(x,Z) calculates rigidity in GV from KEPN in GeV
      FUNCTION R(X,Z)
      IMPLICIT NONE
      REAL(8)::R,X
      INTEGER(4)::Z
      REAL(8), DIMENSION(28) :: A
      INTEGER(4), DIMENSION(28)::N
      COMMON /ATOMIC/N,A
!
      R = SQRT(X*X+2.*X*.938272)/Z !GV from GeV/nucleon!
      END FUNCTION R
!______________________________________________________________________
! LOG10(x) calculates LOG base 10 of X
      FUNCTION LOG10(X)
      IMPLICIT NONE
      REAL(8)::LOG10,X
!
      LOG10 = LOG(X)/LOG(10.0)
!
      END FUNCTION LOG10
!______________________________________________________________________
      FUNCTION FINDU(t)
      IMPLICIT NONE
      REAL(8)::t,FINDU
      INTEGER(4)::yyyy,mm,dd,h,U

      CALL Dates(t,yyyy,mm,dd,h)
!      PRINT*, yyyy,mm,dd,h
      CALL DATE2HP(yyyy,mm,dd,U)
!      print*,U
      FINDU=REAL(U,kind=8)
      END FUNCTION FINDU 
!_______________________________________________________________________
! HP RELATED SUBS 
!
! THESE ARE ADAPTED FROM CARI
! PERMANENT DATABASE IS UNIT 28, USER DBASE IS UNIT 31
      SUBROUTINE DATE2HP(YEAR,MONTH,DAY,HP)
        ! TAKES A DATE AND ASSIGNS AN HP
           
         INTEGER(4), INTENT(IN)::MONTH,DAY,YEAR
         INTEGER(4), INTENT(OUT)::HP
         INTEGER(4)::PHP,CHAR2INT
         INTEGER(4)::MM,YYYY,DD
         CHARACTER(10)::DATESTR
         CHARACTER(3)::DIAGNOSE='NO '
         CHARACTER(10)::INIVAR
         CHARACTER(12)::INIVAL
         CHARACTER(12)::VIEWER 
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT

      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
         OPEN(UNIT=28,FILE='SOLARMOD\MV-DATES.L99',STATUS='OLD')
         OPEN(UNIT=31,FILE='SOLARMOD\MORDATES.2K',STATUS='OLD')
      ELSE   
         OPEN(UNIT=28,FILE='SOLARMOD/MV-DATES.L99',STATUS='OLD')
         OPEN(UNIT=31,FILE='SOLARMOD/MORDATES.2K',STATUS='OLD')
      ENDIF
!       REWIND (28)
!       REWIND (31)

         DO
            READ(31,FMT=2802,ERR=2800,END=2800) DATESTR, PHP
            CALL DATE2YMD(DATESTR,YYYY,MM,DD)
! CONVERT DATESTRING
            IF (MONTH.EQ.MM .AND. YEAR.EQ.YYYY .AND. DAY.EQ.DD) THEN
               HP=PHP
               Close(28)
               Close(31)
            RETURN
         ENDIF       
         ENDDO
2800     CONTINUE !DATE NOT FOUND YET
         DO
            HHP = PHP
            READ(28,FMT=2802,ERR=2801,END=2801) DATESTR, PHP
            CALL DATE2YMD(DATESTR,YYYY,MM,DD) 
! SINCE THERE ARE NO DAYS IN MV-DATES, ACCEPT MATCHING MONTH         
            IF (MONTH.EQ.MM .AND. YEAR.EQ.YYYY) THEN
               HP=PHP
               CLose(28)
               Close(31)
            RETURN
         ENDIF       
         ENDDO
2801     CONTINUE !DATE STILL NOT FOUND, USE MOST RECENT POTENTIAL IN 
!                  THE PERMANENT DBASE
!
         HP=HHP
         Close(28)
         Close(31)
2802     FORMAT(A10,I10)
      END SUBROUTINE
!_______________________________________________________________________&
      SUBROUTINE DATE2YMD(DSTR,Y,M,D)

         CHARACTER(10), INTENT(IN)::DSTR
         INTEGER(4), INTENT(OUT)::Y,M,D 
         INTEGER(4)::CHAR2INT

! CONVERT DATESTRING TO INTEGERS

            IF (DSTR(8:8).EQ.',') THEN !DATE IS MM/YYYY
               M=CHAR2INT(DSTR(1:2),2)
               Y=CHAR2INT(DSTR(4:7),4)
               D=0
            ELSE                       !DATE IS YYYY/MM/DD
               M=CHAR2INT(DSTR(6:7),2)
               Y=CHAR2INT(DSTR(1:4),4)
               D=CHAR2INT(DSTR(9:10),2)
            ENDIF
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
