!
!      Program CARI-7/7A 
!      Fortran coding by Kyle Copeland (KC), Civil Aerospace Medical Insitute
!
!      Flight routines are based on the C6W Fortran code set developed by 
!      Kyle Copeland for the US FAA and for Royal Military College of Canada 
!      course CC-503 (Winter 2012 term)
!
!                        Other Contributions to source code
!      DAP
!      The main lookup engine: uses the pre-analyzed CARI-7 shower data files
!      Coded by KC as part of RMC course TH600 and in support of FAA 2013 NARP 
!      goals. 
!      
!      SOLAR_CR  
!      Solar Event CR spectrum generator library  
!
!      GCR_BO11
!      The Badhwar and O'Neil 2011 GCR spectrum code, modified by KC for 
!      inclusion here, original is from O'Neil.
!
!      GCR_BO14
!      The Badhwar and O'Neil 2014 GCR spectrum code, modified by KC from the BO11 
!      code using NASA/TP-2015-218569 by PM O'Neill, S Golge, and TC Slaba.
!      As of OCT 2016 this is still a work in progress.
!
!      GCR_ISO
!      The 2004/2013 ISO GCR spectrum, coded entirely by KC using data and the 
!      algorithm from the ISO report.
!
!      GEODESIC (combination of NOAA's Forward.for and Reverse.for)
!      Source codes are freely available from NOAA and are slightly modified for 
!      use as libraries by KC here.
!
!      SKY_LIBS
!      Contains: 
!      Hani Anid's vertical magnetic rigidity cutoff correction based on Kp. 
!      Shea and Smarts' calculations of nonvertical magnetic rigidity cutoffs, 
!      Altitude based corrections to cutoffs, limited to 1500 g/sq.cc, by KC. 
!      Aircraft horizon, i.e., which parts of the sky are physically within 
!      line-of-site of the aircraft, all coded by KC. 
!
!      UTILITY
!      Some utility functions and subroutines, coded by KC, and a spline 
!      by Tanmoy Das.   
!
!      CARI-6 
!      Some legacy code has been ported from BASIC to FORTRAN from CARI-6.
!      See Help file for complete list of developers of CARI-6
!
!      To have access to all airports and room to add more, for Windows
!      compile with extra stack space. 
!      The /F3000000 option is enough with ifort: compile with, e.g., 
!      "ifort /F3000000 cari7.f"
!      This does not seem to be an issue in Linux 
!
!      BUG fixes since thesis publication (starting with the most recent)
!      KC 20170404 Fixed a bug in SKYLIBS:FT2GPCMS preventing conversion of alts below 
!          sea level 
!      KC 20161103 Fixed unix/linux/posix path name errors
!      KC 20161030 Added a geometric normalization factor derived from direct comparison 
!          of incident flux to MCNPX generated flux per incident particle at the highest
!          possible altitude
!      KC 20161028 
!         -Fixed shower databases for neutron and photon flux summing error
!         -updated shower dose data to ICRP 116 for particles up through He++
!         -optimized NZ and NA selections for zenith and azimuth dependent calcs
!      KC 20160124 eliminated errors and simplified code for normalizations   
!      KC 20160801 Fixed only East longitudes allowed for new airports bug       
!      KC 20160405 Changed STEPFEET and STEPMIN FROM dim 12 to 20, now matches Cari-6
!      KC 20150217 Elim character spillage from older runs 
!
!      OTHER MODS OF NOTE
!      KC 20170413 remove uncertainty/sigma from normal printed output for locations, 
!         it is still in the diagnostic print of UNIT 40. Since I am not using the 
!         full uncertainty calculation, it is rather pointless and misleading to 
!         include the easily misinterpreted numbers in the normal output. 
!         NOTE: this is difference from CARI-7A in subs ONESPOT and READIT.
!      KC 20161117-1209 added capability to read DEG files like CARI-6M 
!      There are several variables/functions/subroutines that allow the user to play 
!      with the handling of the various inputs:
!         DAP fluxmod: this can be modified to change the incident GCR for any
!             ion
!         DAP zfilter: turns on or off input of primaries of ions Z  
!         DAP FHF: option to convert the highest energy bin to integral flux  
!
!      CARI-7 is a reduced option version of CARI-7A designed for easy maintenance and 
!      use:
!        GCR Spectrum is always the HP modulated ISO spectrum corrected to match ICRU data
!        at high latitudes (CARI-6 HP + 250MV (ADDED INTERNALLY).  
!        Zenith and azimuth resolved cuttoffs and shower depths corrections are not used.

      PROGRAM CARI7

      IMPLICIT NONE
         INTEGER(4)::PARTICLE,DOSEKIND,HOUR
         CHARACTER(10)::INIVAR,YMD1
         CHARACTER(30)::PNAME
         CHARACTER(12)::INIVAL
         CHARACTER(12)::VIEWER 
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE,PTYPE
         INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR
 
         CALL GETINI 
         CALL LOADFTDCCS(OS)
         CALL LOADFORBUSH 
         CALL OPENDATABASES 
         IF (MENUS.EQ.'NO!') THEN !batch mode 
            !all run options preset in DEFAULT.INP and CARI.INI
            OPEN(UNIT=45,FILE='DEFAULT.INP',STATUS='OLD')
            READ(45,'(A10)') YMD1
!            PRINT*, YMD1
            READ(45,'(I4)') HOUR
!            PRINT*, HOUR
            READ(45,'(I4)') PARTICLE
!            PRINT*, PARTICLE
            READ(45,'(I4)') DOSEKIND
!            PRINT*, DOSEKIND
            READ(45,'(A30)') PNAME
!            PRINT*, PNAME
            PTYPE=PNAME(LEN_TRIM(PNAME)-2:LEN_TRIM(PNAME))
!            PRINT*, PTYPE
            CLOSE(45)
            CALL LC2UC(PTYPE,PTYPE,3)
            IF ((PTYPE.EQ.'BIG').OR.(PTYPE.EQ.'big')) THEN 
               CALL RUNBIG !Run a set of flight profiles using default settings
            ELSEIF (PTYPE.EQ.'DEG'.OR.(PTYPE.EQ.'deg')) THEN 
               CALL RUNDEG !Run a waypoint based flight profile using default settings
            ELSE
               CALL READIT !Try to read it as a locations data file
            ENDIF
            STOP ! After analysis of the named file, stop this job
         ENDIF   
         CALL MENU_MAIN !User wants to select options at runtime 
                 
      END PROGRAM CARI7
! END OF PROGRAM CARI7 MAIN MODULE      
!
! Include source code for geodesics
! These code are freely available from NOAA
! These subs use double precision i/o
!  
      INCLUDE 'GEODESIC.FOR' 
!              inverse(in:: lat1, lon1, lat2, lon2
!                     out:: faz, baz, meters)  

!              forward(in:: lat1, lon1, faz, meters
!                     out:: lat2, lon2)
!
      INCLUDE 'DAP.for'
!          Dose at point        
!          FUNCTION DOSE(lat,lon,alt,time,dosekind,gcrmodel,testrun,superposition) 
      INCLUDE 'SKY_LIBS.FOR'  
!          SUB FINDCUT(in::dyear,lat,lon,alt
!                     out::VC,SkyCut(18,18),SkyPass(18,18),SkyWeight(18,18))
!          FUN LINTERP(x1,x2,y1,y2,x)
!          SUB DATES(in:dyear
!                   out:y,m,d,h)
!          SUB dt2dyear(in:y,m,d,uth,utm
!                      out:dyear)           
!                     
!     Include GCR model (ISO LIS modulated by heliocentric potential method)
      INCLUDE 'GCR.FOR'
!
!     A Static library of subs and functions       
      INCLUDE 'Utility.for' ! misc useful subs and functions
!-----------------------------------------------------------------------XXXXX                   
!----6-----------------------------------------------------------------2
! START OF FUNCTIONS AND SUBS TO GET DOSE RATE 
! Written by Kyle Copeland
      SUBROUTINE OPENDATABASES

      IMPLICIT NONE

      CHARACTER(10)::INIVAR
      CHARACTER(12)::INIVAL
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
      INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR

      OPEN(UNIT=1, FILE='DEBUGGER.DAT',STATUS='UNKNOWN') 

!      UNIT 11, Cutoffs\WGRC1965.1X1  -  SKY_LIBS
!      UNIT 12, Cutoffs\IGRF1980.1X1  -  SKY_LIBS
!      UNIT 13, Cutoffs\DGRF1990.1X1  -  SKY_LIBS
!      UNIT 14, Cutoffs\IGRF1995.1X1  -  SKY_LIBS
!      UNIT 15, Cutoffs\IGRF2000.1X1  -  SKY_LIBS
!      UNIT 16, Cutoffs\IGRF2010.1X1  -  SKY_LIBS

!      UNIT 19, DATAIN
!      UNIT 20, DATAOUT
!      UNIT 21, kp_index\kp_index.dat  -  SKY_LIBS

! OPEN PERMANENT AND TRANSIENT DATABASES

      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
!     SET WINDOWS PATHS
        OPEN(UNIT=27,FILE='AIRPORTS\AIRPORTS.DAT',STATUS='OLD')
        OPEN(UNIT=28,FILE='SOLARMOD\MV-DATES.L99',STATUS='OLD')
        OPEN(UNIT=30,FILE='AIRPORTS\CODES',STATUS='OLD')
!     OPEN TRANSIENT DATABASES
        OPEN(UNIT=31,FILE='SOLARMOD\MORDATES.2K',STATUS='OLD')
        OPEN(UNIT=32, FILE='AIRPORTS\NEWPORTS.DAT',STATUS='OLD')
        OPEN(UNIT=40, FILE='DIAGNOSE\CARI7CHK.DAT',STATUS='UNKNOWN') 
      ELSE
! SET LINUX PATHS
        OPEN(UNIT=27,FILE='AIRPORTS/AIRPORTS.DAT',STATUS='OLD')
        OPEN(UNIT=28,FILE='SOLARMOD/MV-DATES.L99',STATUS='OLD')
        OPEN(UNIT=30,FILE='AIRPORTS/CODES',STATUS='OLD')
!     TRANSIENT DATABASES
        OPEN(UNIT=31,FILE='SOLARMOD/MORDATES.2K',STATUS='OLD')
        OPEN(UNIT=32, FILE='AIRPORTS/NEWPORTS.DAT',STATUS='OLD')
        OPEN(UNIT=40, FILE='DIAGNOSE/CARI7CHK.DAT',STATUS='UNKNOWN')
      ENDIF

! OTHER UNITS USED IN THE PROGRAM
!      UNIT 33,'CITY.NDX'
!      UNIT 34,'PORT.NDX'
!      UNIT 96, INVERSE AND FORWARD FILES  
!      UNIT 99, SCRATCH

      CALL MAKE_NDX(91,7000,OS) 
                         ! Without special adjustment, this number is 
                         ! machine dependent and may lead to stack 
                         ! space errors if the airport number exceeds 
                         ! 7000 (seems typical) at runtime. See note at  
                         ! start for IFORT compiling fix using Windows.
                         ! Still need to find the right option in linux        
!
      END SUBROUTINE OPENDATABASES
!                                                                      7
!----6-----------------------------------------------------------------2
!0-9999 
! FUNCTION TO FIND DOSE RATE OF TYPE DT, FROM PARTICLES RAD, AT ALTITUDE DEPTH,
! AT LONGITUDE LON, AT LATITUDE LAT, AT HELIOCENTRIC POTENTIAL HP, ON DATE YMD
!
!
!                RAD=1, NEUTRONS  
!                RAD=2, PHOTONS
!                RAD=3, POSITRONS  
!                RAD=4, ELECTRONS  
!                RAD=5, +MUONS  
!                RAD=6, -MUONS  
!                RAD=7, PROTONS (H-1)+ 
!                RAD=8, +PIONS  
!                RAD=9, -PIONS
!                RAD=10, DEUTERONS (H-2)+
!                RAD=11, TRITONS (H-3)+
!                RAD=12, HELIONS (He-3)++ 
!                RAD=13, ALPHAS (He-4)++
!                RAD=14, LI  
!                RAD=15, B  
!                RAD=16, BE  
!                RAD=17, C  
!                RAD=18, N  
!                RAD=19, O  
!                RAD=20, F  
!                RAD=21, NE  
!                RAD=22, NA  
!                RAD=23, MG  
!                RAD=24, AL  
!                RAD=25, SI  
!                RAD=26, P  
!                RAD=27, S
!                RAD=28, CL  
!                RAD=29, AR 
!                RAD=30, K
!                RAD=31, CA  
!                RAD=32, SC 
!                RAD=33, TI  
!                RAD=34, V  
!                RAD=35, CR 
!                RAD=36, MN 
!                RAD=37, FE 
!                RAD=38, TOTAL For all particles  
! 
!                DT=1, Secondary Particle flux
!                DT=2, ICRP103 EFFECTIVE DOSE
!                DT=3, ICRP_60 EFFECTIVE DOSE
!                DT=4, ICRU AMBIENT DOSE EQUIVALENT, H*(10)
!                DT=5, ESTIMATED WHOLE-BODY ABSORBED DOSE 
!
!                                                                      7
!----6-----------------------------------------------------------------2
! END OF SUBS AND FUNCTIONS TO GET DOSE RATE
      SUBROUTINE GETINI
         CHARACTER(10)::INIVAR
         CHARACTER(12)::VIEWER,INIVAL 
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE,SPA,PTYPE
         CHARACTER(1)::JUNK,GCRA,NVRC
         INTEGER(4)::GCR,SP 
         REAL(8)::DEFAULTMV
         LOGICAL::TF                   
          
         COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
         COMMON /INI/SP,GCR
         COMMON /NVC_TF/TF,NVRC
 
!         TF=.FALSE. !DO NOT USE NONVERTICAL CUTOFFS
!         TF=.TRUE. !USE NONVERTICAL CUTOFFS
!         DIAGNOSE='NO!' 
         DIAGNOSE='YES'
          
           OPEN (UNIT=99,FILE='CARI.INI',STATUS='OLD',ACTION='READ')
       DO 
            READ(99,9101) INIVAR,JUNK,INIVAL
!            READ(99,9100,ADVANCE='NO') INIVAR
!           IF (INIVAR(1:2).EQ.'HP') THEN 
!               READ(99,9103,ADVANCE='YES') INIVAR,JUNK,DEFAULTMV
!
         IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR, JUNK, INIVAL
!            ENDIF            
         IF (INIVAR(1:2).EQ.'VI') THEN 
               VIEWER=INIVAL
               IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR,'=',VIEWER
         ENDIF          
         IF (INIVAR(1:2).EQ.'OU') THEN 
               OUTPUT=INIVAL(1:4)
               IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR,'=',OUTPUT
         ENDIF          
         IF (INIVAR(1:2).EQ.'OS') THEN 
               OS=INIVAL(1:5)
               IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR,'=',OS
         ENDIF          
         IF (INIVAR(1:2).EQ.'DI') THEN 
               DISPLAY=INIVAL(1:3)
               IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR,'=',DISPLAY
         ENDIF          
         IF (INIVAR(1:2).EQ.'ME') THEN 
               MENUS = INIVAL(1:3)
               IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR, '=',MENUS
         ENDIF                   
         IF (INIVAR(1:2).EQ.'NV') THEN 
            NVRC = INIVAL(1:1)
            IF ((NVRC.EQ.'1').OR.(NVRC.EQ.'3')) THEN 
               TF=.TRUE.
            ELSE
               TF=.FALSE.
            ENDIF      
            IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR, '=',NVRC
         ENDIF                   
         IF (INIVAR(1:2).EQ.'SU') THEN 
               SPA = INIVAL(1:3)
               IF (SPA.EQ.'YES') THEN
                  SP=1
               ELSE
                  SP=0
               ENDIF
               IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR, '=', SP
         ENDIF                   
         IF (INIVAR(1:2).EQ.'GC') THEN 
               GCRA = INIVAL(1:1)
               IF (GCRA.EQ.'7') THEN       !Custom spectrum 
                  GCR=7
               ELSEIF (GCRA.EQ.'6') THEN   ! Sep 89 SPE
                  GCR=6 
               ELSEIF (GCRA.EQ.'5') THEN   ! Feb 56 SPE
                  GCR=5 
               ELSEIF (GCRA.EQ.'4') THEN   ! GCR: ISO LIS  w/HP Modulation 
                  GCR=4 
               ELSEIF (GCRA.EQ.'3') THEN   ! GCR: BO14
                  GCR=3 
               ELSEIF (GCRA.EQ.'2') THEN   ! GCR: BO11
                  GCR=2 
               ELSE
                  GCR=1                    ! ISO 2004/2013 GCR   
               ENDIF    
               IF (DIAGNOSE.EQ.'YES')  WRITE(40,*) INIVAR, '=', GCR
         ENDIF                   
         IF (INIVAR(1:2).EQ.'EN') EXIT
       ENDDO     
         DIAGNOSE = 'NO!'
!CARI-7 ignores some choices and always uses these values
         GCR=4
         NVRC='0'
         SP=0
         TF=.FALSE.
9101     FORMAT(A10,A2,A12)
         CLOSE(99) 
      END SUBROUTINE      

!                                                                      7
!----6-----------------------------------------------------------------2
! END OF SUBS AND FUNCTIONS TO GET DOSE RATE
!
!10000-19999
! BEGIN SUBS AND FUNCTIONS TO MANIPULATE/EVALUATE FLIGHT PROFILE DATA
! A. BIG FILES (SEVERAL FLIGHTS, ASSUMES GREAT CIRCLE ROUTES, AIRPORTS)
!     1. CREATE NEW BIG FILE
!     2. CREATE NEW SHORT LIST FOR AN EXISTING BIG FILE
!     3. DELETE EXISTING BIG FILE
!     4. DELETE EXISTING SHORT LIST
!     5. ADD FLIGHTS TO EXISTING BIG FILE
!     6. REVIEW FLIGHTS IN AN EXISTING BIG FILE
!     7. REVIEW FLIGHTS IN AN EXISTING SHORT LIST
!     8. REMOVE FLIGHTS FROM AN EXISTING BIG FILE
! B. DEG FILES (WAYPOINT FILES FOR SINGLE FLIGHTS)
!     1. CREATE A NEW DEG FILE
!     2. DELETE AN EXISTING DEG FILE
!     3. CREATE/EDIT/DELETE A FLIGHT LIST OF DEG FILES
!        A. ADD ALL DEG FILES IN CURRENT DIRECTORY TO LIST
!        B. ADD ONLY SELECTED DEG FILES
!        C. REMOVE FILES FROM AN EXISTING LIST 
! C. EVALUATE FLIGHT PROFILES OUTPUT
!     1. RESULTS TO ARCHIVE 
!     2. RESULTS TO SCREEN
!     3. RESULTS TO ARCHIVE AND SCREEN
! D. EVALUATE FLIGHTS
! RDBIGFLIGHT READS THE WHOLE FLIGHT PROFILE FOR ANALYSYS BY SUBROUTINE
! BIG_FLT_DOSE 
      SUBROUTINE RDBIGFLT(FLTNAME,YMD,OPORT,DPORT,CLIMBMIN,NOSTEPS,     & 
     &   STEPFEET,STEPMIN,DESCMIN,CRUISEMIN,TRIPMIN,BAD)    
      CHARACTER(30)::FLTNAME
      CHARACTER(10)::FLTDATE,YMD
      CHARACTER(6)::OPORT,DPORT,ENDSTR   
      INTEGER(4)::  NOSTEPS, STEPFEET(12), STEPMIN(12), DESCMIN 
      INTEGER(4):: CLIMBMIN, CRUISEMIN, TRIPMIN, I, J, BAD
      REAL(8) :: DIST,TRIPMILES,MAXFEET
      
      CHARACTER(3)::DIAGNOSE, LOCALDIAG, DIAG

      DIAGNOSE='no!'
      LOCALDIAG='YES'
      IF ((DIAGNOSE.EQ.'YES').OR.(LOCALDIAG.EQ.'YES')) DIAG='YES' 

      BAD=0
      MAXFEET = 0
        CLIMBMIN = 0
        CRUISEMIN = 0 
      DESCMIN = 0
      TRIPMIN = 0
      NOSTEPS = 0
      DO I = 1,12
         STEPFEET(I)=0
         STEPMIN(I)=0       
      ENDDO
      READ(18,78801,ERR=78810,END=78820) FLTNAME
      READ(18,78802,ERR=78810,END=78820) FLTDATE
      IF (DIAG.EQ.'YES') WRITE(40,*) 'READING', FLTNAME
      IF (DIAG.EQ.'YES') WRITE(40,*) 'READING DATE', FLTDATE
      FLTDATE=ADJUSTL(FLTDATE)
      IF ((LEN(TRIM(FLTDATE))).EQ.7) THEN
           IF(FLTDATE(3:3).EQ.'/') THEN !EXPECTED FORMAT
             YMD(1:4)=FLTDATE(4:7)
            YMD(5:5)='/'
            YMD(8:8)='/'
            YMD(6:7)=FLTDATE(1:2)
            YMD(9:10)='00'
          ELSEIF (FLTDATE(5:5).EQ.'/') THEN !EUROPEAN FORMAT
             YMD(1:4)=FLTDATE(4:7)
            YMD(5:5)='/'
            YMD(8:8)='/'
            YMD(6:7)=FLTDATE(1:2)
            YMD(9:10)='00'
          ELSE
            BAD = 1
          ENDIF
      ELSEIF ((LEN(TRIM(FLTDATE))).EQ.10) THEN
           IF (FLTDATE(5:5).EQ.'/') THEN !EXPECTED LONG FORMAT
             YMD(1:10)=FLTDATE(1:10)
          ELSE
            BAD = 1
          ENDIF
      ELSE
         BAD = 1 
      ENDIF
      IF (BAD.EQ.1) THEN
         DO 
            READ(18,*,ERR=78810,END=78820) ENDSTR
            IF (SCAN(ENDSTR,'-').NE.0 .OR. SCAN(ENDSTR,'_').NE.0) EXIT 
         ENDDO     
         IF (DIAG.EQ.'YES') WRITE(40,*) 'EXITING RDBIGFLT, BAD=',BAD 
         RETURN
      ENDIF   
      READ(18,*) OPORT
      OPORT = ADJUSTL(OPORT)
      IF (DIAG.EQ.'YES') WRITE(40,*) 'OPORT ', OPORT
      READ(18,*) DPORT
      DPORT = ADJUSTL(DPORT)
      IF (DIAG.EQ.'YES') WRITE(40,*) 'DPORT ', DPORT
      READ(18,*) NOSTEPS 
      READ(18,*) CLIMBMIN
      IF (DIAG.EQ.'YES') WRITE(40,*) 'STEPS ', NOSTEPS
      IF (DIAG.EQ.'YES') WRITE(40,*) 'CLIMB MIN ', CLIMBMIN
      DO I = 1, NOSTEPS
          READ(18,*,ERR=78810,END=78820) STEPFEET(I), STEPMIN(I)
          IF(STEPFEET(I) > MAXFEET ) MAXFEET = STEPFEET(I) 
            ! store highest altitude
          IF (DIAG.EQ.'YES') WRITE(40,*) STEPFEET(I), STEPMIN(I) 
      END DO
      READ(18,*,ERR=78810,END=78820) DESCMIN
      READ(18,*,ERR=78810,END=78820) ENDSTR
      IF (DIAG.EQ.'YES') WRITE(40,*) DESCMIN, ENDSTR 
!     Altitude limit is removed in CARI-7, KC        
      DO I = 1, NOSTEPS
         CRUISEMIN = CRUISEMIN + STEPMIN(I)
      END DO
!          *** TripTime is ground to ground time in minutes ***
      TRIPMIN = CRUISEMIN + DESCMIN + CLIMBMIN
78801 FORMAT(A30)
78802 FORMAT(A10)
      IF (DIAG.EQ.'YES') WRITE(40,*) 'EXITING RDBIGFLT, BAD=',BAD 
      RETURN      
78810 BAD=2
      IF (DIAG.EQ.'YES') WRITE(40,*) 'EXITING RDBIGFLT, BAD=',BAD 
      RETURN     
78820 BAD=3
      IF (DIAG.EQ.'YES') WRITE(40,*) 'EXITING RDBIGFLT, BAD=',BAD 
      END SUBROUTINE 
!--------1---------2---------3---------4---------5---------6---------7--
      FUNCTION BADSTR(I)
      CHARACTER(45)::BADSTR
      INTEGER(4)::I

      SELECT CASE (I)     
           CASE(1)
             BADSTR='  CANNOT READ DATE...SKIPPING THIS FLIGHT    '
           CASE(2)
             BADSTR='  CANNOT READ ALTS...SKIPPING THIS FLIGHT    '
           CASE(3)
             BADSTR='  REACHED END OF BIG FILE                    '
           CASE(4)
             BADSTR='  ALT > 60,000 FT...SKIPPING THIS FLIGHT     '
           CASE(5)
             BADSTR='  TRIP TIME >1996 MIN...SKIPPING THIS FLIGHT '
           CASE(6)
             BADSTR='  BAD ORIGIN PORT...SKIPPING THIS FLIGHT     '
           CASE(7)
             BADSTR='  BAD DESTINATION PORT...SKIPPING THIS FLIGHT'
         CASE DEFAULT
             BADSTR='  BAD PROFILE DATA...SKIPPING THIS FLIGHT    '
      END SELECT
      END FUNCTION
!--------1---------2---------3---------4---------5---------6---------7--
!+--------------- use Geodesic survey method for trip distance ------+
      SUBROUTINE USE_INVERSE (la1, lo1, la2, lo2, miles, FAZ)
              
          REAL(8) :: METERS2MILES,meters,miles
          REAL(8), INTENT(IN)::la1,lo1,la2,lo2
          REAL(8) :: FAZ,dlat1,dlat2,dlon1,dlon2,BAZ,metres
 
         CHARACTER(3)::DIAGNOSE
         DIAGNOSE='no!'
 
         dlat1=DBLE(la1)       
         dlat2=DBLE(la2)       
         dlon1=DBLE(lo1)       
         dlon2=DBLE(lo2)       

         call inverse(dlat1,dlon1,dlat2,dlon2,FAZ,BAZ,metres)

         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'FROM LAT= ',la1,' LON=',lo1 
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'TO LAT= ',la2,' LON=',lo2 
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'ALONG FAZ= ',FAZ, ' IS M=', &
     &    metres, ' METERS' 
        
         meters=real(metres,kind=4)     
         METERS2MILES = 1. / 1852.
           miles = meters * METERS2MILES
          
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) ' ',miles,' NAUTICAL MILES' 

      END SUBROUTINE USE_INVERSE
!--------1---------2---------3---------4---------5---------6---------7--
!+--use Geodesic survey method for intermediate coordinates--+
!   for best accurracy, miles should be less than 100, BUT 
!   for CARI a few millimeters is not an issue, results are 
!   reversible at all reasonable distances            

      SUBROUTINE USE_FORWARD (la1, lo1, la2, lo2, miles, FAZ)
        
!        LATITUDE IS NORTH=POSITIVE, SOUTH=NEGATIVE
!        LONGITUDE IS EAST POSITIVE, WEST=NEGATIVE
         
         REAL(8), INTENT(IN) :: miles, la1, lo1
         REAL(8), INTENT(OUT) ::  la2, lo2
         REAL(8) :: meters
         REAL(8), INTENT(IN) :: FAZ
         DOUBLE PRECISION :: dlat1,dlon1,dlat2,dlon2,metres       
 
         CHARACTER(3)::DIAGNOSE

         DIAGNOSE='NO!'

         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'FORWARD MILES=',miles,     & 
     &       'ALONG FAZ= ',FAZ,  ' FROM '
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'LAT= ',la1,' LON= ',lo1 
 
         meters = 1852.* miles
         metres=DBLE(meters)
         dlat1=DBLE(la1)
         dlon1=DBLE(lo1)
         call forward(dlat1,dlon1,FAZ,metres,dlat2,dlon2)
         lo2=real(dlon2, kind=8)
         IF (lo2<-180.0) lo2=lo2+360.0
         la2=real(dlat2, kind=8)

         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'IS LAT= ',la2,' LON= ',lo2 

      END SUBROUTINE USE_FORWARD
!-----------------------------------------------------------------------&        
      SUBROUTINE BIG_FLT_DOSE(FLTNAME,OPORT,DPORT,YMD2,HR,FLIGHTDOSE,   &
     &                             RADIAT,DOSEKIND,BAD)
      ! SUB GENERATES A SET OF WAYPOINTS BASED ON THE FLIGHT PROFILE
      ! THEN CALCULATES A FLIGHT DOSE
         IMPLICIT NONE
         REAL(8)::LAT(2000),LON(2000),DEPTH(2000),DR,DT(2000)

         CHARACTER(30)::DAPORT,OAPORT,DCITY,OCITY
         CHARACTER(10)::YMD
         CHARACTER(6)::OPORT,DPORT,ENDSTR
         CHARACTER(1)::ONS,OEW,DEW,DNS     

         CHARACTER(10), INTENT(IN)::YMD2 
         INTEGER(4), INTENT(IN)::HR, RADIAT, DOSEKIND

         INTEGER(4)::Y,M,D,H,HP,DK
         INTEGER(4)::ALLSTEPS, NUMLOCS
         INTEGER(4)::CLIMBSTEPS, DESCSTEPS 
         INTEGER(4)::NOSTEPS, STEPFEET(20), STEPMIN(20), DESCMIN 
! Changed STEPFEET and STEPMIN FROM dim 12 to 20, NOW matchES Cari-6, KC 20160405
         INTEGER(4)::CLIMBMIN, CRUISEMIN, TRIPMIN, I, J 
          
         CHARACTER(30), INTENT(OUT)::FLTNAME
         INTEGER(4), INTENT(OUT)::BAD
         REAL(8), INTENT(OUT)::FLIGHTDOSE
          
         REAL(8) :: OALT, DALT, NSM, EWM, SPEED, STEPDIST(12)
         REAL(8) :: DIST,TRIPMILES, RHP,TOTALDOSE,ALT,CAS
         REAL(8) :: OLAT, OLON, DLAT, DLON, miles, HOUR
         REAL(8) :: FAZ, DOSE, T, CLIMBDIST, DOWNDIST, CRUISEDIST
          
         CHARACTER(3)::DIAGNOG
         CHARACTER(12)::VIEWER 
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
         INTEGER(4)::SP,GCR
         LOGICAL::DIAG,FRACK
         LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE

      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE                   
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR

      DIAGNOG='YES'
!      DIAGNOG='no!'
         BAD=0 
         IF ((DIAGNOSE=='YES').or.(DIAGNOG=='YES')) THEN
            DIAG=.TRUE.
         ELSE
            DIAG=.FALSE.
         ENDIF
!         PRINT*, DIAGNOG, DIAGNOSE, DIAG !diagnostic
         DK=RADIAT+38*(DOSEKIND-1) 
! fyi     
!                RADIAT=1, NEUTRONS  
!                RADIAT=2, PROTONS  
!                RADIAT=3, PHOTONS
!                RADIAT=4, POSITRONS  
!                RADIAT=5, ELECTRONS  
!                RADIAT=6, +MUONS  
!                RADIAT=7, -MUONS  
!                RADIAT=8, +PIONS  
!                RADIAT=9, -PIONS
!                ...  
!                RADIAT=37, Fe  
!                RADIAT=38, TOTAL  
! 
!
!                DOSEKIND=1, Secondary Particle Flux
!                DOSEKIND=2, ICRP103 EFFECTIVE DOSE
!                DOSEKIND=3, ICRP_60 EFFECTIVE DOSE
!                DOSEKIND=4, ICRU AMBIENT DOSE EQUIVALENT, H*(10)
!                DOSEKIND=5, ESTIMATED WHOLE-BODY ABSORBED DOSE 

         IF (DIAGNOSE.EQ.'YES') THEN
             WRITE(40,*) 'ANALYZING A FLIGHT' 
         ENDIF    
! INITIALIZE ARRAYS
         DO I=1,2000
            DT(I)=0 
            DEPTH(I)=0        
            LAT(I)=0
            LON(I)=0
         ENDDO
         DR=0
         FLIGHTDOSE=0
         NUMLOCS=0
         ALLSTEPS=0
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'CALLING RDBIGFLT ' 
         CALL RDBIGFLT(FLTNAME,YMD,OPORT,DPORT,CLIMBMIN,NOSTEPS,STEPFEET& 
     &                 ,STEPMIN,DESCMIN,CRUISEMIN,TRIPMIN,BAD)

         IF (BAD.NE.0) RETURN 
         IF (DIAG) THEN 
            WRITE(40,*) 'RDBIGFLT RETURNS FLIGHTNAME: ',FLTNAME
            WRITE(40,*) 'RDBIGFLT RETURNS ORIGIN PORT: ', OPORT
            WRITE(40,*) 'RDBIGFLT RETURNS DESTINATION: ', DPORT
            WRITE(40,*) 'RDBIGFLT RETURNS CLIMB TIME (MINS): ',CLIMBMIN
            WRITE(40,*) 'RDBIGFLT RETURNS NOSTEPS: ',NOSTEPS
            DO I = 1, NOSTEPS
                WRITE(40,*) 'STEP ',I, STEPMIN(I), STEPFEET(I)
            ENDDO
            WRITE(40,*) 'RDBIGFLT RETURNS CRUISE TIME (MIN): ',CRUISEMIN
            WRITE(40,*) 'RDBIGFLT RETURNS DESCENT TIME (MIN): ',DESCMIN
            WRITE(40,*) 'RDBIGFLT RETURNS FLT TIME (MIN): ' ,TRIPMIN
            WRITE(40,*) ' '
         ENDIF
                   
         FRACK=.FALSE.

         CALL PORT_INFO(OPORT,OAPORT,OCITY,ONS,OLAT,OEW,OLON,OALT,FRACK)
         IF (FRACK) THEN !Abort this profile, missing airport data
            IF (DIAG) WRITE(40,*)'Missing airport data for', OPORT
            IF (DIAG) WRITE(40,*)'Aborting flight'
            BAD=6
            RETURN
         ENDIF
         CALL PORT_INFO(DPORT,DAPORT,DCITY,DNS,DLAT,DEW,DLON,DALT,FRACK)
         IF (FRACK) THEN !Abort this profile, missing airport data
            IF (DIAG) WRITE(40,*)'Missing airport data for', DPORT
            IF (DIAG) WRITE(40,*)'Aborting flight'
            BAD=7
            RETURN
         ENDIF

                
         IF (DIAG) THEN 
            WRITE(40,*) 'ORIGIN INFORMATION'
            WRITE(40,*) 'CITY AND PORTNAMES ', OCITY, OAPORT
            WRITE(40,*) 'READING LAT AND LON ', ONS,OLAT,OEW,OLON
            WRITE(40,*) 'READING ALTITUDE ', OALT, 'FEET'
            WRITE(40,*) ' '
            WRITE(40,*) 'DESTINATION INFORMATION'
            WRITE(40,*) 'CITY AND PORTNAMES ', DCITY, DAPORT
            WRITE(40,*) 'READING LAT AND LON ', DNS,DLAT,DEW,DLON
            WRITE(40,*) 'READING ALTITUDE ', DALT, 'FEET'
         ENDIF          
         IF (ONS.EQ.'N') THEN
             NSM=1.
         ELSE  
             NSM=-1.
         ENDIF
         IF (OEW.EQ.'E') THEN
             EWM=1.
         ELSE  
             EWM=-1.
         ENDIF
         OLAT=OLAT*NSM
         OLON=OLON*EWM
         LAT(1)=OLAT
         LON(1)=OLON
         IF (DNS.EQ.'N') THEN
             NSM=1.
         ELSE  
             NSM=-1.
         ENDIF
         IF (DEW.EQ.'E') THEN 
             EWM=1.
         ELSE  
             EWM=-1.
         ENDIF
         DLAT=DLAT*NSM
         DLON=DLON*EWM
          
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'CALLING INVERSE'
         CALL USE_INVERSE(OLAT,OLON,DLAT,DLON,miles,FAZ)
!        SPEED FOR ESTIMATING SEGMENT DISTANCES           
         SPEED = miles/(CRUISEMIN+0.5*(DESCMIN+CLIMBMIN))
         IF (DIAG) WRITE(40,*) 'CALCULATED SPEED miles/min ', SPEED
         
         CLIMBSTEPS = CLIMBMIN+1
         DESCSTEPS = DESCMIN+1

!        CLIMBING STEPS, USE ONE STEP PER MINUTE, CENTERED ON 1/2 STEP
         IF (DIAG) WRITE(40,*) 'CALCULATING CLIMB'
         CALL FT2GPCMS(OALT,DEPTH(1))
         DT(1)=0.5
         CAS=(STEPFEET(1)-OALT)/CLIMBMIN
         DO I = 1,CLIMBMIN  !ALT,TIME,LOCS DURING CLIMB
            DT(I+1)=1.
            CLIMBDIST = 0.5*SPEED*I
            CALL FT2GPCMS(OALT+CAS*I,DEPTH(I+1))
            CALL USE_FORWARD(LAT(1),LON(1),LAT(I+1),LON(I+1),CLIMBDIST, &
     &                       FAZ)
            IF(I.EQ.1) THEN  
               IF(DIAG) WRITE(40,*)I,DEPTH(I),0.0,FAZ,LAT(1),LON(1)
            ENDIF
            IF(DIAG) WRITE(40,*)I+1,DEPTH(I+1),CLIMBDIST,FAZ,LAT(I+1),  &
     &       LON(I+1)
         ENDDO
         DT(CLIMBMIN+1)=0.5
         
!        CRUISING STEPS [CLIMBMIN+2,CLIMBMIN+CRUISEMIN+2,]

         IF (DIAG) WRITE(40,*) 'CALCULATING CRUISE COORDS'
         IF (TRIPMIN>1996) THEN 
            BAD=5  !CAN'T USE MINUTES AS STEPS, TRIP TOO LONG
            TOTALDOSE= -5.0
            RETURN
         ENDIF 
         DO I = CLIMBMIN+2,CLIMBMIN+CRUISEMIN+2  !ALT,TIME,LOCS DURING CRUISES
            DT(I)=1.
            IF(I.EQ.CLIMBMIN+2) DT(I)=0.0 
            ! must not weight the transition step, time already used in climb 
            CRUISEDIST = CLIMBDIST+SPEED*(I-(CLIMBMIN+1)-0.5) !AVE POSITION FOR EACH MIN
            CALL ALTNOW(CLIMBMIN,NOSTEPS,STEPFEET,STEPMIN,DESCMIN,      &
     &                 CRUISEMIN,I,ALT)         
            CALL FT2GPCMS(ALT,DEPTH(I))
            CALL USE_FORWARD(LAT(1),LON(1),LAT(I),LON(I),CRUISEDIST,FAZ)
            IF (DIAG) WRITE(40,*)I,DEPTH(I),CRUISEDIST,FAZ,LAT(I),LON(I)
         ENDDO
         J=CLIMBMIN+CRUISEMIN+3
         ALLSTEPS=TRIPMIN+4 
!        DESCENDING STEPS, USE ONE STEP PER MINUTE, CENTERED ON 1/2 STEP
         IF (DIAG) WRITE(40,*) 'CALCULATING DESCENT'
         CALL FT2GPCMS(DALT,DEPTH(ALLSTEPS))
         DT(J)=0.5
         LAT(J)=LAT(J-1) !descent starts at location and alt of end of cruise 
         LON(J)=LON(J-1)
         DEPTH(J)=DEPTH(J-1)
         CAS=(STEPFEET(NOSTEPS)-DALT)/(DESCMIN)
         DO I = 1,DESCMIN  !ALT,TIME,LOCS DURING DESCENT
            J=J+1 
            DT(J)=(REAL(DESCMIN,kind=8)-1.)/DESCMIN !avoid overweighting descent segments
            DOWNDIST = CRUISEDIST + 0.5*SPEED*I
            CALL FT2GPCMS(STEPFEET(NOSTEPS)-CAS*I,DEPTH(J))
            CALL USE_FORWARD(LAT(1),LON(1),LAT(J),LON(J),DOWNDIST,FAZ)
            IF (DIAG) WRITE(40,*)J, DEPTH(J),DOWNDIST,FAZ,LAT(J),LON(J)
         ENDDO           
         DT(J+1)=0.5
         CALL FT2GPCMS(DALT,DEPTH(J+1))
         LAT(J+1)=DLAT
         LON(J+1)=DLON
         IF (DIAG) WRITE(40,*)J+1, DEPTH(J+1),0.5*DOWNDIST,FAZ,LAT(J+1),&
     &                  LON(J+1), ALLSTEPS
             
! WHICH DATE TO USE? YMD2='0000/00/00' FOR USE PROFILE INFO
         IF (YMD2.NE.'0000/00/00') THEN
            YMD=YMD2   
            H=HR/100
            HOUR=REAL(HR-H*100)*1.9026E-06 
            !minutes to add to T at start, expressed as fraction of a year
         ELSE
            H=0
            HOUR=0.
         ENDIF
         
         CALL DATE2YMD(YMD,Y,M,D)
         
         CALL YMDH2T(Y,M,D,H,T)

!         CALL DATE2HP(Y,M,D,HP)
                
         IF(M.EQ.0) then 
            YEARLYAVE=.TRUE.
         ELSE
            YEARLYAVE=.FALSE.
         ENDIF
         IF(D.EQ.0) then 
            MONTHLYAVE=.TRUE.
         ELSE
            MONTHLYAVE=.FALSE.
         ENDIF
         IF((H.EQ.0).AND.(HOUR.EQ.0.0))then 
            DAILYAVE=.TRUE.
         ELSE
            DAILYAVE=.FALSE.
         ENDIF
         T=T+HOUR
         IF (DIAG) WRITE(40,*) 'CALCULATING FLIGHT DOSE ', T,DK,GCR,SP
         IF (DIAG) WRITE(40,*)'STEP, LAT, LON, ALT, DOSERATE, '         &
     &                                       //'TIME, CUMDOSE' 
!         IF (DIAG) WRITE(*,*) 'STEP, LAT, LON, ALT, DOSERATE, '         &
!     &                                       //'TIME, CUMDOSE' 
!        STOP !DIAGNOSTIC
         DO I = 1, ALLSTEPS 
            IF (DIAGNOSE.EQ.'YES') THEN !GLOBAL DIAGNOSTIC OUTPUT IS ON
               DR=DOSE(LAT(I),LON(I),DEPTH(I),T,DK,GCR,1,SP)
               
            ELSE ! LOCAL ONLY DIAGNOSTICS ONLY
               DR=DOSE(LAT(I),LON(I),DEPTH(I),T,DK,GCR,0,SP)
            ENDIF
            T=T+1.9026E-06*DT(I) 
            ! add DT minutes each time for accurate Forbush and kp
            FLIGHTDOSE=DR*DT(I)/60.+FLIGHTDOSE
            IF (DIAG) WRITE(40,*)I, LAT(I), LON(I),                     &
     &                              DEPTH(I), DR,  DT(I), FLIGHTDOSE    
!            IF (DIAG) WRITE(*,*) I, LAT(I), LON(I),                     &
!     &                              DEPTH(I), DR,  DT(I), FLIGHTDOSE    
         ENDDO
         IF (DIAG) WRITE(40,*) FLIGHTDOSE,RADIAT,DOSEKIND,SP
                
        END SUBROUTINE BIG_FLT_DOSE
!_______________________________________________________________________&
!        
! ALTNOW FINDS ALTITUDE DURING THE FLIGHT AT ANY TIME DURING THE CRUISE      
      SUBROUTINE ALTNOW(CLIMBMIN,N,STEPFEET,STEPMIN,DESCMIN,CRUISEMIN,I,&
     &                  ALT)
      INTEGER(4)::N
      INTEGER(4)::CLIMBMIN,STEPFEET(N),STEPMIN(N),DESCMIN,CRUISEMIN    
      INTEGER(4)::I,J,K,M
      REAL(8)::ALT
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
      INTEGER(4)::SP,GCR
                   
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR
           
      K=1
      J=STEPMIN(K)
      DO M=1,CRUISEMIN
         IF(M.LT.J) THEN 
            ALT=REAL(STEPFEET(K),KIND=8)
            ELSE
               ALT=REAL(STEPFEET(K),KIND=8)
               K=K+1
               J=J+STEPMIN(K)
            ENDIF          
            IF (M+CLIMBMIN+2.EQ.I) EXIT        
         ENDDO
         IF (DIAGNOSE=='YES') WRITE(40,*) I, ALT
      END SUBROUTINE ALTNOW
!_----------------------------------------------------------------------&
!     1. EVALUATE AT BIG FILE
!     USER CAN OVERRIDE INTERNAL DATE & SELECT DOSE OUTPUT
      SUBROUTINE RUNBIG!(FILENAME,YMD1,HP1,PARTICLE,DOSEKIND)
!
      IMPLICIT NONE

      INTEGER(4)::PARTICLE,DOSEKIND,HOUR,BAD,WP,I
        
      REAL(8) :: FLIGHTDOSE
        
      CHARACTER(1)::C
      CHARACTER(10)::YMD1,PARTSTR
      CHARACTER(30)::OUTFILENAME
      CHARACTER(39)::DKSTR
      CHARACTER(30)::FILENAME,FLIGHTNAME
      CHARACTER(45)::BADSTR
      CHARACTER(6)::OPORT,DPORT
      LOGICAL::FRACK
! COMMON BLOCK VARIABLES

         CHARACTER(10)::INIVAR
         CHARACTER(12)::VIEWER,INIVAL 
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE,DIAGALL
         INTEGER(4)::SP,GCR 

      LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR

! LOCALLY OVERRIDE DIAGNOSE
          DIAGALL=DIAGNOSE
!         DIAGNOSE='YES'
          DIAGNOSE='NO!'
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'MENUS,OS,DISPLAY,VIEWER,OUTPUT'
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) MENUS,OS,DISPLAY,VIEWER,OUTPUT
      IF (MENUS.EQ.'NO!') THEN
          OPEN(UNIT=45,FILE='DEFAULT.INP',STATUS='OLD')
            READ(45,'(A10)') YMD1
            READ(45,*) HOUR
            READ(45,*) PARTICLE
            READ(45,*) DOSEKIND
            READ(45,'(A30)') FILENAME
          CLOSE(45)
      ELSE
          CALL PICKBIG(FILENAME)
          CALL PICK_FLT_OPTIONS(PARTICLE,GCR,DOSEKIND,SP,YMD1,HOUR)
      ENDIF
        
      OPEN (UNIT=18,FILE=FILENAME,STATUS='OLD')
!
!     READ PROFILES FROM BIG FILE 1 AT A TIME, PRINT DOSES TO OUTPUT
!
      WP=SCAN(FILENAME,'.',BACK=.TRUE.)
      OUTFILENAME(WP:WP+4)='.OUT'
      OUTFILENAME(1:WP-1)=FILENAME(1:WP-1)
      DO I=12,(WP+5),-1
         !pad unused characters in OUTFILENAME with blanks to avoid  
         !character spillage from older runs, KC 20150217
         OUTFILENAME(I:I)=" " 
      ENDDO
      OPEN (UNIT=17,FILE=OUTFILENAME,STATUS='UNKNOWN')
      WRITE(17,*)'FLIGHTS FROM '//FILENAME
      WRITE(17,*)'PARTICLE ',PARTICLE,' DOSEKIND',DOSEKIND
      WRITE(17,*)'GCR MODEL ',GCR,' SUPERPOSITION ',SP,' YMD ',YMD1

      DO
       IF(DIAGNOSE.EQ.'YES') WRITE(40,*)'YMD1, HOUR, PARTICLE, DOSEKIND'     
       IF(DIAGNOSE.EQ.'YES') WRITE(40,*) YMD1,HOUR,PARTICLE,DOSEKIND
       CALL BIG_FLT_DOSE(FLIGHTNAME,OPORT,DPORT,YMD1,HOUR,FLIGHTDOSE,   &
     &                   PARTICLE,DOSEKIND,BAD) 
       IF (BAD.EQ.0) THEN 
          WRITE(17,18922)FLIGHTNAME,OPORT,DPORT,FLIGHTDOSE,             &
     &                   DKSTR(DOSEKIND),PARTSTR(PARTICLE)
          IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'FINISHED ', FLIGHTNAME    
          WRITE(*,*) 'FINISHED ', FLIGHTNAME    
       ELSE
          WRITE(17,18923) BADSTR(BAD)
          IF (BAD.EQ.3) THEN
             DIAGNOSE=DIAGALL ! AVOID SWITCH OF GLOBAL DIAGNOSTIC PRINTS
             CLOSE(17)
             CLOSE(18)
             WRITE(40,*) 'AT END OF FILE'      
             EXIT !QUIT RUN AT END OF FILE
          ENDIF
       ENDIF         
      ENDDO 
      CLOSE(17)
      CLOSE(18)
      WRITE(40,*) 'runbig finished'      
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) 'runbig finished'
      DIAGNOSE=DIAGALL
18922 FORMAT(A30,1X,A6,1X,A6,1X,ES11.4,1X,A42,1X,A10)
18923 FORMAT(A45)
      END SUBROUTINE RUNBIG
!_______________________________________________________________________&  
      SUBROUTINE PICK_FLT_OPTIONS(PARTICLE,GCR,DOSEKIND,SP,YMD1,HOUR)

      IMPLICIT NONE

      INTEGER(4)::PARTICLE,DOSEKIND,GCR,HOUR,SP
      CHARACTER(1)::C
      CHARACTER(10)::YMD1

! COMMON BLOCK VARIABLES

         CHARACTER(12)::VIEWER 
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE

      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT

      CALL CLS
      CALL MENUHEADER
        !     TYPE, 10 (n,p,gamma,e+,e-,mu+,mu-,pi+,pi-,total)  
18902 WRITE(*,*)'      SELECT RADIATION.'
      PRINT*,   '      <0> TOTAL       <10> DEUTERONS   <20>F    <30>K '
      PRINT*,   '      <1> NEUTRONS    <11> TRITONS     <21>Ne   <31>Ca'
      PRINT*,   '      <2> PHOTONS     <12> HELIONS     <22>Na   <32>Sc'
      PRINT*,   '      <3> ELECTRONS   <13> ALPHAS      <23>Mg   <33>Ti'
      PRINT*,   '      <4> POSITRONS   <14> Li          <24>Al   <34>V '
      PRINT*,   '      <5> NEG. MUONS  <15> Be          <25>Si   <35>Cr'
      PRINT*,   '      <6> POS. MUONS  <16> B           <26>P    <36>Mn'
      PRINT*,   '      <7> PROTONS     <17> C           <27>S    <37>Fe'
      PRINT*,   '      <8> POS. PIONS  <18> N           <28>Cl         '
      PRINT*,   '      <9> NEG. PIONS  <19> O           <29>Ar         '
      PRINT*,   '      Enter 0-37 and press <enter>.'
      READ(*,*) PARTICLE
      IF (PARTICLE.LT.0 .OR. PARTICLE.GT.37) THEN
         CALL CLS   
           PRINT*,   '      Entry must be 0 to 37 '
           GOTO 18902
      ENDIF
      IF (PARTICLE.EQ.0) PARTICLE=38
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'PARTICLE',PARTICLE
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) 'PARTICLE',PARTICLE
18903 PRINT*, ' '
!      WRITE(*,*)'      SELECT COSMIC RAY MODEL              '
!      PRINT*,   '      <1> GCR: ISO TS15390:2004-MSU-NYMMIK '
!      PRINT*,   '      <2> GCR: BADHWAR-ONEILL 2011 USING INTL SUNSPOTS'
!      PRINT*,   '      <3> USE DEFAULT  '
!      PRINT*,   '      <4> GCR: HP MODULATED ISO  '
!      PRINT*,   '      <5> SPE: LaRC SEP 1989 EVENT TOTAL '
!      PRINT*,   '      <6> SPE: LaRC FEB 1956 EVENT TOTAL '
!      PRINT*,   '      <7> USER SELECTED: GCR_MODELS\MY_MODEL.OUT'
!      READ(*,*) GCR
       GCR=4 !CARI-7 always uses 4
!      IF (GCR.LT.1 .OR. GCR.GT.7) THEN
!         CALL CLS   
!           PRINT*,   '      Entry must be 1 to 7 '
!           GOTO 18903
!      ENDIF
18904 PRINT*, ' '
      WRITE(*,*)'      SELECT DOSE TYPE'
      PRINT*,   '      <1> SECONDARY PARTICLE FLUENCE'
      PRINT*,   '      <2> ICRP PUB 103 EFFECTIVE DOSE'
      PRINT*,   '      <3> ICRP PUB 60 EFFECTIVE DOSE'
      PRINT*,   '      <4> ICRU H*(10) AMBIENT DOSE EQUIVALENT'
      PRINT*,   '      <5> WHOLE BODY ABSORBED DOSE'
      READ*, DOSEKIND
      IF (DOSEKIND.LT.1 .OR. DOSEKIND.GT.5) THEN
         CALL CLS   
           PRINT*,   '      Entry must be 1 to 5 '
           GOTO 18904
      ENDIF
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'DOSEKIND',DOSEKIND 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) 'DOSEKIND',DOSEKIND 

!18905 WRITE(*,*)'      USE THE SUPERPOSITION APPROXIMATION <Y,N> ?'
!      WRITE(*,*)'                   (DEFAULT IS N)'
!      CALL READKB(C,1)
!      IF (C.EQ.'Y') then
!         SP=1
!      ELSE
         SP=0
!      ENDIF
!    Assume using dates in profiles until otherwise indicated
      HOUR=0
      YMD1='0000/00/00'

18906 WRITE(*,*)'      USE DATES FROM FLIGHT DATA FILE <Y,N> ?'
      CALL READKB(C,1)
      IF (C.EQ.'N') THEN
18910     WRITE(*,*)
          WRITE(*,*) '      ENTER A DATE FOR FLIGHT(S) <YYYY/MM/DD>'
          WRITE(*,*) '      (USE MM=00 AND DD=00 TO USE YEARLY AVERAGE)'
          WRITE(*,*) '      (USE DD=00 TO USE MONTHLY AVERAGE)         '
          READ(*,18921) YMD1
          IF (YMD1(5:5).NE.'/' .OR. YMD1(8:8).NE.'/' ) THEN
             WRITE(*,*) 'WRONG FORMAT: RE-ENTER'
             GOTO 18910
          ENDIF
          IF (YMD1(9:10).NE.'00') THEN !User entered a specific date
18912     WRITE(*,*) ''
          WRITE(*,*) '      ENTER A UT TAKE-OFF TIME FOR FLIGHT(S)?'
          WRITE(*,*) '      (0 for no FORBUSH or GEOMAGNETIC STORMS, or'
          WRITE(*,*) '      range 1-2400. NOTE: if not 0, calculation  '
          WRITE(*,*) '      times will be significantly increased).    '
          READ(*,18920) HOUR
             IF (HOUR.LT.0 .OR. HOUR.GT.2400) THEN
                CALL CLS   
                PRINT*,    '      Entry must be 0 to 2400 '
                GOTO 18912
             ENDIF
          ENDIF
      ENDIF  

18920 FORMAT(I4)
18921 FORMAT(A10)

      END SUBROUTINE PICK_FLT_OPTIONS
!_______________________________________________________________________        
!FUTURE WORK....
!     2. EVALUATE A SHORT LIST
!_______________________________________________________________________        
!     3. EVALUATE A WAYPOINT(DEG) FILE
!     USER CAN OVERRIDE INTERNAL DATE & MUST SELECT DOSE OUTPUT
      SUBROUTINE RUNDEG!(FILENAME,YMD1,HP1,PARTICLE,DOSEKIND)
!
      IMPLICIT NONE

      INTEGER(4)::PARTICLE,DOSEKIND,HOUR,BAD,WP,I
        
      REAL(8) :: FLIGHTDOSE
        
      CHARACTER(1)::C
      CHARACTER(10)::YMD1,PARTSTR
      CHARACTER(30)::PNAME,OUTFILE1NAME,OUTFILE2NAME
      CHARACTER(39)::DKSTR
      CHARACTER(30)::FILENAME,FLIGHTNAME
      CHARACTER(45)::BADSTR
      CHARACTER(6)::OPORT,DPORT
      LOGICAL::TF
! COMMON BLOCK VARIABLES

         CHARACTER(10)::INIVAR
         CHARACTER(12)::VIEWER,INIVAL 
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE,DIAGALL
         CHARACTER(1)::NVRC
         INTEGER(4)::SP,GCR 

      LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR
      COMMON /NVC_TF/TF,NVRC

! LOCALLY OVERRIDE DIAGNOSE
          DIAGALL=DIAGNOSE
!         DIAGNOSE='YES'
          DIAGNOSE='NO!'
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'MENUS,OS,DISPLAY,VIEWER,OUTPUT'
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) MENUS,OS,DISPLAY,VIEWER,OUTPUT
      IF (MENUS.EQ.'NO!') THEN
!          FILENAME='DEFAULT.DEG'
          OPEN(UNIT=45,FILE='DEFAULT.INP',STATUS='OLD')
            
            READ(45,'(A10)') YMD1
            READ(45,*) HOUR
            READ(45,*) PARTICLE
            READ(45,*) DOSEKIND
            READ(45,'(A30)') FILENAME

          CLOSE(45)
      ELSE
          CALL PICKDEG(FILENAME)
          CALL PICK_FLT_OPTIONS(PARTICLE,GCR,DOSEKIND,SP,YMD1,HOUR)
      ENDIF
        
      OPEN (UNIT=18,FILE=FILENAME,STATUS='OLD')
!
!     READ PROFILE FROM DEG FILE, PRINT DOSE TO OUTPUT
!
      WP=SCAN(FILENAME,'.',BACK=.TRUE.)
      OUTFILE1NAME(WP:WP+4)='.SUM'
      OUTFILE1NAME(1:WP-1)=FILENAME(1:WP-1)
      OUTFILE2NAME(WP:WP+4)='.DAT'
      OUTFILE2NAME(1:WP-1)=FILENAME(1:WP-1)
      DO I=12,(WP+5),-1
         !pad unused characters in OUTFILENAME with blanks to avoid  
         !character spillage from older runs, KC 20150217
         OUTFILE1NAME(I:I)=" " 
         OUTFILE2NAME(I:I)=" " 
      ENDDO
      OPEN (UNIT=17,FILE=OUTFILE1NAME,STATUS='UNKNOWN')
      OPEN (UNIT=19,FILE=OUTFILE2NAME,STATUS='UNKNOWN')


!      DO !for running multiple flights
       IF(DIAGNOSE.EQ.'YES') WRITE(40,*)'YMD1, HOUR, PARTICLE, DOSEKIND'     
       IF(DIAGNOSE.EQ.'YES') WRITE(40,*) YMD1,HOUR,PARTICLE,DOSEKIND
         WRITE(19,*)' LAT(N)   LON(E)   DEPTH   STEP  DOSE RATE  TOTAL' &
     &//' DOSE'
       CALL DEG_FLT_DOSE(FLIGHTNAME,YMD1,HOUR,FLIGHTDOSE,               &
     &                   PARTICLE,DOSEKIND,BAD) 
       IF (BAD.EQ.0) THEN 
         WRITE(17,*)'FLIGHT FROM '//FILENAME
         WRITE(17,19924)' GCR MODEL:',GCR,'  TRANSPORT AND CUTOFFS: ',  &
     &                  NVRC, '   SUPERPOSITION:',SP                                    
         WRITE(17,*) 'DATE: ',YMD1,' HOUR: ', HOUR 
         WRITE(19,19924)' GCR MODEL:',GCR,'  TRANSPORT AND CUTOFFS: ',   &
     &                  NVRC, '   SUPERPOSITION:',SP                                    
         WRITE(19,*) 'DATE: ',YMD1,' HOUR: ', HOUR 
         WRITE(17,19922) FLIGHTDOSE, PARTSTR(PARTICLE), DKSTR(DOSEKIND)
         WRITE(19,19922) FLIGHTDOSE, PARTSTR(PARTICLE), DKSTR(DOSEKIND)
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'FINISHED ', FLIGHTNAME    
         WRITE(*,*) 'FINISHED ', FILENAME    
       ELSE
         WRITE(19,19923) BADSTR(BAD)
         WRITE(17,19923) BADSTR(BAD)
         IF (BAD.EQ.3) THEN
            DIAGNOSE=DIAGALL ! AVOID SWITCH OF GLOBAL DIAGNOSTIC PRINTS
!             EXIT !QUIT RUN AT END OF FILE 
!                   not needed for single flight dose but may be useful 
!                   this function is expanded to run multiple DEG files
         ENDIF
       ENDIF         
!      ENDDO ! not needed for single flight dose
      CLOSE(17) ! summary report file
      CLOSE(18) ! flight data input file 
      CLOSE(19) ! step by step report file
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'rundeg finished'      
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) 'rundeg finished'
      DIAGNOSE=DIAGALL
19922 FORMAT(ES11.4,1X,A10,1X,A42)
19923 FORMAT(A45)
19924 FORMAT(A10,I4,A25,A1,A17,I4)
      END SUBROUTINE RUNDEG
!_______________________________________________________________________

      SUBROUTINE DEG_FLT_DOSE(FLTNAME,YMD2,HR,FLIGHTDOSE,RADIAT,DOSEKIND&
     &                        ,BAD)
      ! SUB GENERATES A SET OF WAYPOINTS BASED ON THE FLIGHT PROFILE
      ! THEN CALCULATES A FLIGHT DOSE
         IMPLICIT NONE
         REAL(8)::LAT(2000),LON(2000),DEPTH(2000),DR
         REAL(8)::WPLAT(2000),WPLON(2000),WPFT(2000)
         INTEGER(4)::WPMIN(2000),TMARK,LATD,LOND

         CHARACTER(60)::HEADERS
         CHARACTER(10)::YMD,FLTDATE
         CHARACTER(6)::ENDSTR
         CHARACTER(1)::LEW,LNS     

         CHARACTER(10)::YMD2 
         INTEGER(4), INTENT(IN)::HR, RADIAT, DOSEKIND

         INTEGER(4)::Y,M,D,H,HP,DK
         INTEGER(4)::ALLSTEPS,COMMA
         INTEGER(4)::CLIMBSTEPS, DESCSTEPS 
         INTEGER(4)::NOSTEPS, STEPFEET, SEGMIN, SEGMENTS  
         INTEGER(4)::TRIPMIN, I, J ,N
          
         CHARACTER(30), INTENT(OUT)::FLTNAME
         INTEGER(4), INTENT(OUT)::BAD
         REAL(8), INTENT(OUT)::FLIGHTDOSE
          
         REAL(8) :: NS, EW, NSM, EWM, FT, SPEED, LATM, LONM
         REAL(8) :: DIST,TRIPMILES,RHP,TOTALDOSE,ALT,STEPMIN
         REAL(8) :: OLAT, OLON, DLAT, DLON, miles, HOUR
         REAL(8) :: ROAC, FAZ, DOSE, T
          
         CHARACTER(3)::DIAGNOG
         CHARACTER(12)::VIEWER 
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
         INTEGER(4)::SP,GCR
         LOGICAL::DIAG,FRACK,TFNS,TFEW
         LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE

      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE                   
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR

      DIAGNOG='YES'
!      DIAGNOG='no!'
         BAD=0 
         IF ((DIAGNOSE=='YES').or.(DIAGNOG=='YES')) THEN
            DIAG=.TRUE.
         ELSE
            DIAG=.FALSE.
         ENDIF
!         PRINT*, DIAGNOG, DIAGNOSE, DIAG !diagnostic
         DK=RADIAT+38*(DOSEKIND-1) 
! fyi     
!                RADIAT=1, NEUTRONS  
!                RADIAT=2, PROTONS  
!                RADIAT=3, PHOTONS
!                RADIAT=4, POSITRONS  
!                RADIAT=5, ELECTRONS  
!                RADIAT=6, +MUONS  
!                RADIAT=7, -MUONS  
!                RADIAT=8, +PIONS  
!                RADIAT=9, -PIONS
!                ...  
!                RADIAT=37, Fe  
!                RADIAT=38, TOTAL  
! 
!
!                DOSEKIND=1, Secondary Particle Flux
!                DOSEKIND=2, ICRP103 EFFECTIVE DOSE
!                DOSEKIND=3, ICRP_60 EFFECTIVE DOSE
!                DOSEKIND=4, ICRU AMBIENT DOSE EQUIVALENT, H*(10)
!                DOSEKIND=5, ESTIMATED WHOLE-BODY ABSORBED DOSE 

! INITIALIZE ARRAYS
      DO I=1,2000
         DEPTH(I)=0        
         LAT(I)=0
         LON(I)=0
         WPLAT(I)=0
         WPLON(I)=0
         WPFT(I)=0
         WPMIN(I)=0
      ENDDO
! INITIALIZE COUNTERS
      J=1
      DR=0
      FLIGHTDOSE=0
      ALLSTEPS=0
      TRIPMILES=0
! BEGIN THE REAL WORK
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'READING DEG FILE '
      READ(18,'(A60)') HEADERS   ! expect '(A10,A30)' but can vary
      COMMA=INDEX(HEADERS,',')
      FLTDATE=HEADERS(1:COMMA-1)
      FLTNAME=TRIM(ADJUSTL(HEADERS(COMMA+1:60)))
!      PRINT*, FLTDATE,FLTNAME !diagnostic
      IF (YMD2.EQ.'0000/00/00') THEN !use date from file
         IF(FLTDATE(3:3).EQ.'/') THEN !OLD MM/YYYY FORMAT
            YMD(1:4)=FLTDATE(4:7)
            YMD(5:5)='/'
            YMD(8:8)='/'
            YMD(6:7)=FLTDATE(1:2)
            YMD(9:10)='00'
         ELSEIF ((FLTDATE(5:5).EQ.'/').AND.(FLTDATE(8:8).EQ.' ')) THEN !EUROPEAN FORMAT
            YMD(1:4)=FLTDATE(1:4)
            YMD(5:5)='/'
            YMD(6:7)=FLTDATE(6:7)
            YMD(8:8)='/'
            YMD(9:10)='00'
         ELSEIF ((FLTDATE(5:5).EQ.'/').AND.(FLTDATE(8:8).EQ.'/')) THEN !EUROPEAN FORMAT
            YMD=FLTDATE
         ELSE
            BAD = 1
         ENDIF
         IF (BAD.EQ.1) THEN
            IF (DIAG.EQ.'YES') WRITE(40,*) 'EXITING RDDEGFLT, BAD=',BAD 
            RETURN
         ENDIF
      ENDIF
! WHICH DATE TO USE? YMD2='0000/00/00' FOR USE PROFILE INFO
      IF (YMD2.NE.'0000/00/00') THEN
         YMD=YMD2   
         H=HR/100
         HOUR=REAL(HR-H*100)*1.9026E-06 
         !minutes to add to T at start, expressed as fraction of a year
      ELSE
         H=0
         HOUR=0.
      ENDIF
           YMD2=YMD !flip printed value for return to RUNDEG
      CALL DATE2YMD(YMD,Y,M,D)
      CALL YMDH2T(Y,M,D,H,T)
      CALL DATE2HP(Y,M,D,HP)

      IF(M.EQ.0) then 
         YEARLYAVE=.TRUE.
      ELSE
         YEARLYAVE=.FALSE.
      ENDIF
      IF(D.EQ.0) then 
         MONTHLYAVE=.TRUE.
      ELSE
         MONTHLYAVE=.FALSE.
      ENDIF
      IF((H.EQ.0).AND.(HOUR.EQ.0.0))then 
         DAILYAVE=.TRUE.
      ELSE
         DAILYAVE=.FALSE.
      ENDIF
!      PRINT*, YMD, FLTNAME                  ! '(A10,A30)'
      READ(18,'(A60)') HEADERS
!         FRACK=.TRUE.
      DO WHILE (.NOT.EOF(18))
         READ(18,*) LATD, LATM, LNS, LOND, LONM, LEW, FT, TMARK
            TFNS=(LNS.EQ.'N').OR.(LNS.EQ.'S')
            TFEW=(LEW.EQ.'E').OR.(LEW.EQ.'W')
         IF (TFNS.AND.TFEW) THEN 
               WPLAT(J) = (LATD+(LATM/60.))*NS(LNS)
               WPLON(J) = (LOND+(LONM/60.))*EW(LEW)
               WPFT(J) = FT
               WPMIN(J) = TMARK
               J=J+1
         ELSE
               BAD=9 !data in wrong format somewhere in the profile
               TOTALDOSE=-5.0
               RETURN
         ENDIF
      ENDDO !J-1 is last good data point before EOF
! 19930    FORMAT(2(I3,F5.2,A1),2I8)
      SEGMENTS = J-2
      TRIPMIN = WPMIN(J-1)          
      IF (TRIPMIN.GT.2000) THEN 
         BAD=5  !CAN'T USE MINUTES AS STEPS, TRIP TOO LONG
         TOTALDOSE= -5.0
         RETURN
      ENDIF 
      DO N = 1,SEGMENTS 
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'CALLING INVERSE'
         OLAT=WPLAT(N)!lat at start of segment
         OLON=WPLON(N)!lon at start of segment
         DLAT=WPLAT(N+1)!lat at end of segment
         DLON=WPLON(N+1)!lon at end of segment 
         SEGMIN=(WPMIN(N+1)-WPMIN(N)) !minutes in segment
         CALL USE_INVERSE(OLAT,OLON,DLAT,DLON,miles,FAZ) !FAZ is forward azimuth
         TRIPMILES=TRIPMILES+miles  
!        SPEED FOR ESTIMATING SEGMENT DISTANCES           
         SPEED = miles/REAL(SEGMIN)
         IF (DIAG) WRITE(40,*) 'CALCULATED SPEED miles/min ', SPEED
!        RATE OF ALTITUDE CHANGE
         ROAC = REAL(WPFT(N+1)-WPFT(N))/REAL(SEGMIN)
         IF (DIAG) WRITE(40,*) 'CALCULATED ALT CHANGE ft/min ', ROAC
! BUILD A MINUTE BY MINUTE PROFILE OF THIS SEGMENT OF THE FLIGHT
!                                        
!        USE ONE STEP PER MINUTE, CENTERED ON 1/2 STEP 
!
         IF (DIAG) WRITE(40,*) 'FLIGHT COORDS FOR SEGMENT', N
         IF (DIAG) WRITE(40,*)'I,LAT(I),LON(I),DEPTH(I),DIST,FAZ'
         DO I = WPMIN(N)+1,WPMIN(N+1) !STEPS FOR SEGMENT N
            STEPMIN=REAL(I-WPMIN(N))-0.50
            ALT=REAL(WPFT(N))+ROAC*STEPMIN
            DIST = SPEED*STEPMIN !AVE POSITION FOR EACH MIN IN THIS SEGMENT
            CALL FT2GPCMS(ALT,DEPTH(I))
            CALL USE_FORWARD(OLAT,OLON,LAT(I),LON(I),DIST,FAZ)
            IF (DIAG) WRITE(40,*)I,LAT(I),LON(I),DEPTH(I),DIST,FAZ
!               WRITE(*,*)I,LAT(I),LON(I),DEPTH(I),DIST,FAZ
         ENDDO
      ENDDO
             
      T=T+HOUR
      IF (DIAG) WRITE(40,*) 'CALCULATING FLIGHT DOSE ',T,DK,GCR,SP
      IF (DIAG) WRITE(40,*)'STEP, LAT, LON, ALT, DOSERATE, TIME, CUMDOS'
      IF (DIAG) WRITE(*,*)'STEP, LAT, LON, ALT, DOSERATE, TIME, CUMDOSE'
!        STOP !DIAGNOSTIC
      DO I = 1, TRIPMIN 
         IF (DIAGNOSE.EQ.'YES') THEN !GLOBAL DIAGNOSTIC OUTPUT IS ON
            DR=DOSE(LAT(I),LON(I),DEPTH(I),T,DK,GCR,1,SP)
         ELSE ! LOCAL ONLY DIAGNOSTICS ONLY
            DR=DOSE(LAT(I),LON(I),DEPTH(I),T,DK,GCR,0,SP)
         ENDIF
            ! add 1 minutes each time for accurate Forbush and kp
         T=T+1.9026E-06 
         FLIGHTDOSE=DR/60.+FLIGHTDOSE
         WRITE(19,19926) LAT(I),LON(I),DEPTH(I),I,DR,FLIGHTDOSE
         IF (DIAG) WRITE(40,*)I, LAT(I), LON(I),                        &
     &                           DEPTH(I), DR,  1, FLIGHTDOSE    
         IF (DIAG) WRITE(*,*) I, LAT(I), LON(I),                        &
     &                           DEPTH(I), DR,  1, FLIGHTDOSE    
      ENDDO
      IF (DIAG) WRITE(40,*) FLIGHTDOSE,RADIAT,DOSEKIND,SP
19926 FORMAT(2(F8.4,1X),F8.3,1X,I4,2(1x,ES11.4))
      END SUBROUTINE DEG_FLT_DOSE
! FUTURE WORK...
!    4. EVALUATE A LIST OF DEGFILES
!      
!     ...
! *****************************************************************    
! END OF SUBS AND FUNCTIONS TO MANIPULATE/EVALUATE FLIGHT PROFILE DATA
! *****************************************************************
!_______________________________________________________________________
!20000-29999
! BEGIN SUBS AND FUNCTIONS TO CREATE AND EVALUATE FILES OF SINGLE LOCATIONS
!     1. CREATE NEW FILE OF SINGLE LOCATIONS
!     2. ADD NEW LOCATIONS TO AN EXISTING FILE
!     3. EVALUATE A FILE OF SINGLE LOCATIONS
!        A. RESULTS TO ARCHIVE 
!        B. RESULTS TO SCREEN
!        C. RESULTS TO ARCHIVE AND SCREEN
!     4. DELETE A FILE OF SINGLE LOCATIONS
!     5. EVALUATE A SINGLE LOCATION
!        A. RESULTS TO ARCHIVE 
!        B. RESULTS TO SCREEN
!        C. RESULTS TO ARCHIVE AND SCREEN
! END SINGLE LOCATION SUBS AND FUNCTIONS
!-----------------------------------------------------------------------&
!20000
      SUBROUTINE ONESPOT
! CALCULATE A DOSE RATE AT A SINGLE LOCATION
!
      CHARACTER(4)::Y,H
      CHARACTER(1)::FG,NS,EW,C
      CHARACTER(2)::M,D,LAD,LAM,LAS,LOM,LOS
      CHARACTER(3)::LOD
      CHARACTER(5)::ALTF
      CHARACTER(10)::YMD
      CHARACTER(12)::LATIN,LONIN
      CHARACTER(10)::PARTSTR
      CHARACTER(32)::DSTR(6)

      CHARACTER(3)::DIAGNOG,DIAGNOSE,DIAG

      INTEGER(4)::YEAR,MONTH,DAY,HOUR,PARTICLE
      INTEGER(4)::GCR, SP, HP, DT, DK, GCRM
      INTEGER(4)::CHAR2INT
      REAL(8)::LAT,LON,GM,W,N,RM,RS,SIGMA
      REAL(8)::CHAR2REAL,DOSE,T,VC
      REAL(8)::DOSERATE,G,RALT,FEET
      LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /ISOSTATS/SIGMA 
      COMMON /VCUT/VC 

      YEAR=0;MONTH=0;DAY=0
      DIAGNOG='YES'
      IF ((DIAGNOG.EQ.'YES').OR.(DIAGNOSE.EQ.'YES')) DIAG='YES' 
        
      WRITE(40,*) 'CALLING ONESPOT'
! STEP 1, USER INPUTS LOCATION, PICKS DOSE RATE OUTPUT
! USE DATE

      YEARLYAVE=.FALSE. 
      MONTHLYAVE=.FALSE. 
      DAILYAVE=.FALSE. 
       
      WRITE(*,*)'     What year?'
      CALL READKB(Y,4);YEAR=CHAR2INT(Y,4)
      WRITE(*,*)'     What month (ENTER 0 FOR YEARLY AVERAGE)?'
      CALL READKB(M,2);MONTH=CHAR2INT(M,2)
      IF (MONTH.EQ.0) YEARLYAVE=.TRUE.
      CALL INT2CHAR(MONTH,M,2)
      IF (YEARLYAVE) THEN
         DAY=0
         D=M
         HOUR=0
         H=M
         GOTO 20000
      ENDIF
      WRITE(*,*)'     What day (ENTER 0 FOR MONTHLY AVERAGE)?'
      CALL READKB(D,2);DAY=CHAR2INT(D,2)
      IF (DAY.EQ.0) MONTHLYAVE=.TRUE. 
      CALL INT2CHAR(DAY,D,2)
      IF (MONTHLYAVE) THEN
         HOUR=0
         H=D
         GOTO 20000
      ENDIF
      WRITE(*,*)'     What hour (e.g. 14) in UT?'
      WRITE(*,*)'    (ENTER 0 FOR DAILY AVERAGE)?'
      CALL READKB(H,2);HOUR=CHAR2INT(H,2)
      IF (HOUR.EQ.24) DAILYAVE=.TRUE. 
      CALL INT2CHAR(HOUR,H,2)

20000 CONTINUE
      YMD(1:4)=Y
      YMD(6:7)=M
      YMD(9:10)=D
        
      WRITE(*,*)' '
      WRITE(*,*)'     LATITUDE'
      WRITE(*,*)'     NORTH OR SOUTH <N/S>?'
      CALL READKB(NS,1)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) NS 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) NS
      IF (NS.EQ.'N') THEN 
         N=1.
      ELSEIF(NS.EQ.'S')THEN
         N=-1.
      ELSE
         WRITE(*,*) '     ENTRY MUST BE ''N'' OR ''S''!'
         GOTO 20000
      ENDIF        
      WRITE(*,*)'      DEGREES LATITUDE (UP TO 12 CHARACTERS).'
      CALL READKB(LATIN,12)
      IF (SCAN(LATIN,'.').NE.0) THEN 
      IF (DIAG.EQ.'YES') WRITE(40,*) 'USING FRACTIONAL DEGREES' 
      IF (DIAG.EQ.'YES') WRITE(*,*)  'USING FRACTIONAL DEGREES'
      LAT=CHAR2REAL(LATIN,12)
      ELSE       
      LAD=TRIM(LATIN(1:2))
      IF (DIAG.EQ.'YES') WRITE(40,*) LAD 
      IF (DIAG.EQ.'YES') WRITE(*,*) LAD 
      WRITE(*,*)'      MINUTES (UP TO 2 CHARACTERS).'
      CALL READKB(LAM,2)
      IF (DIAG.EQ.'YES') WRITE(40,*) LAM 
      IF (DIAG.EQ.'YES') WRITE(*,*) LAM 
      WRITE(*,*)'      SECONDS (UP TO 2 CHARACTERS).'
      CALL READKB(LAS,2)
      IF (DIAG.EQ.'YES') WRITE(40,*) LAS 
      IF (DIAG.EQ.'YES') WRITE(*,*) LAS 
      RD=CHAR2INT(LAD,2)/1.     
      RM=CHAR2INT(LAM,2)/60.    
      RS=CHAR2INT(LAS,2)/3600.     
      LAT=N*(RD+RM+RS)
      ENDIF
20001 WRITE(*,*)'      LONGITUDE'
      WRITE(*,*)'      EAST OR WEST <E/W>.'
      CALL READKB(EW,1)
      IF (DIAG.EQ.'YES') WRITE(40,*) EW 
      IF (DIAG.EQ.'YES') WRITE(*,*) EW 
      IF (EW.EQ.'E') THEN 
           W=1.
      ELSEIF (EW.EQ.'W')THEN !KC 20160801 Fixed only E allowed bug
         W=-1.
      ELSE
         WRITE(*,*) '     ENTRY MUST BE ''E'' OR ''W''!'
         GOTO 20001
      ENDIF
      WRITE(*,*)'      DEGREES (UP TO 12 CHARACTERS).'
      CALL READKB(LONIN,12)
      IF (SCAN(LONIN,'.').NE.0) THEN 
      IF (DIAG.EQ.'YES') WRITE(40,*) 'USING FRACTIONAL DEGREES' 
      IF (DIAG.EQ.'YES') WRITE(*,*)  'USING FRACTIONAL DEGREES'
      LON=CHAR2REAL(LONIN,12)
      ELSE    
      LOD=TRIM(LONIN(1:3))
      IF (DIAG.EQ.'YES') WRITE(40,*) LOD 
      IF (DIAG.EQ.'YES') WRITE(*,*) LOD
      RD=CHAR2INT(LOD,3)/1.     
      WRITE(*,*)'      MINUTES (UP TO 2 CHARACTERS).'
      CALL READKB(LOM,2)
      IF (DIAG.EQ.'YES') WRITE(40,*) LOM
      IF (DIAG.EQ.'YES') WRITE(*,*) LOM
      RM=CHAR2INT(LOM,2)/60.    
      WRITE(*,*)'      SECONDS (UP TO 2 CHARACTERS).'
      CALL READKB(LOS,2)
      IF (DIAG.EQ.'YES') WRITE(40,*) LOS 
      IF (DIAG.EQ.'YES') WRITE(*,*) LOS
      RS=CHAR2INT(LOS,2)/3600.     
      LON=W*(RD+RM+RS)
      ENDIF
      PRINT*, ' '
20002 WRITE(*,*)'      ALTITUDE'
      WRITE(*,*)'      ARE UNITS FT, KM, OR G/CM^2 <F/K/G>.'
      CALL READKB(FG,1)
      WRITE(*,*)'      ENTER ALTITUDE.'
      READ*, RALT
      IF (DIAG.EQ.'YES') WRITE(40,*) RALT 
      IF (DIAG.EQ.'YES') WRITE(*,*) RALT 
           
      IF (FG.EQ.'F') THEN 
         CALL FT2GPCMS(RALT,G)
         IF (DIAG.EQ.'YES') PRINT*, 'CONVERTED ',RALT,' ft TO ',G
         IF (DIAG.EQ.'YES') WRITE(40,*)'CONVERTED ',RALT,' ft TO ',G
      ELSEIF (FG.EQ.'K') THEN
           G=KM2DEPTH(RALT)
         IF (DIAG.EQ.'YES') PRINT*, 'CONVERTED ',RALT,' km TO ',G
         IF (DIAG.EQ.'YES') WRITE(40,*)'CONVERTED ',RALT,' km TO ',G
      ELSEIF (FG.EQ.'KG') THEN
           G=RALT
      ELSE
         GOTO 20002
      ENDIF     
      PRINT*, ' '
      CALL DATE2HP(YEAR,MONTH,DAY,HP)
!      PRINT*, ' HELIOCENTRIC POTENTIAL USED WILL BE ',HP,' MV'  
      IF (DIAG.EQ.'YES') WRITE(40,*) ' HELIOCENTRIC POTENTIAL USED WILL &
     & BE ',HP,' MV'
!     TYPE, 38 (n,gamma,e-,e+,mu-,mu+,p, pi-,pi+,...,total)
! 1 = ICRP Pub. 103 Effective dose from neutrons
! 2 = ICRP Pub. 103 Effective dose from photons 
! 3 = ICRP Pub. 103 Effective dose from electrons 
! 4 = ICRP Pub. 103 Effective dose from positrons 
! 5 = ICRP Pub. 103 Effective dose from neg. muons 
! 6 = ICRP Pub. 103 Effective dose from pos. muons 
! 7 = ICRP Pub. 103 Effective dose from protons 
! 8 = ICRP Pub. 103 Effective dose from neg. pions 
! 9 = ICRP Pub. 103 Effective dose from pos. pions
! 10 = ICRP Pub. 103 Effective dose from deuterons
! 11 = ICRP Pub. 103 Effective dose from tritons
! 12 = ICRP Pub. 103 Effective dose from helions
! 13 = ICRP Pub. 103 Effective dose from alphas
! 14 = ICRP Pub. 103 Effective dose from Li ions
! 15 = ICRP Pub. 103 Effective dose from B ions
! 16 = ICRP Pub. 103 Effective dose from Be ions
! 17 = ICRP Pub. 103 Effective dose from C ions
! 18 = ICRP Pub. 103 Effective dose from N ions
! 19 = ICRP Pub. 103 Effective dose from O ions
! 20 = ICRP Pub. 103 Effective dose from F ions
! 21 = ICRP Pub. 103 Effective dose from Ne ions
! 22 = ICRP Pub. 103 Effective dose from Na ions
! 23 = ICRP Pub. 103 Effective dose from Mg ions
! 24 = ICRP Pub. 103 Effective dose from Al ions
! 25 = ICRP Pub. 103 Effective dose from Si ions
! 26 = ICRP Pub. 103 Effective dose from P ions
! 27 = ICRP Pub. 103 Effective dose from S ions
! 28 = ICRP Pub. 103 Effective dose from Cl ions
! 29 = ICRP Pub. 103 Effective dose from Ar ions
! 30 = ICRP Pub. 103 Effective dose from K ions
! 31 = ICRP Pub. 103 Effective dose from Ca ions
! 32 = ICRP Pub. 103 Effective dose from Sc ions
! 33 = ICRP Pub. 103 Effective dose from Ti ions
! 34 = ICRP Pub. 103 Effective dose from V ions
! 35 = ICRP Pub. 103 Effective dose from Cr ions
! 36 = ICRP Pub. 103 Effective dose from Mn ions
! 37 = ICRP Pub. 103 Effective dose from Fe ions  
20003 WRITE(*,*)'      SELECT RADIATION.'
      PRINT*,   '      <0> TOTAL       <10> DEUTERONS   <20>F    <30>K '
      PRINT*,   '      <1> NEUTRONS    <11> TRITONS     <21>Ne   <31>Ca'
      PRINT*,   '      <2> PHOTONS     <12> HELIONS     <22>Na   <32>Sc'
      PRINT*,   '      <3> ELECTRONS   <13> ALPHAS      <23>Mg   <33>Ti'
      PRINT*,   '      <4> POSITRONS   <14> Li          <24>Al   <34>V '
      PRINT*,   '      <5> NEG. MUONS  <15> Be          <25>Si   <35>Cr'
      PRINT*,   '      <6> POS. MUONS  <16> B           <26>P    <36>Mn'
      PRINT*,   '      <7> PROTONS     <17> C           <27>S    <37>Fe'
      PRINT*,   '      <8> POS. PIONS  <18> N           <28>Cl         '
      PRINT*,   '      <9> NEG. PIONS  <19> O           <29>Ar         '
      PRINT*,   '      Enter 0-37 and press <enter>.'
!      CALL READKB(C,1)
!         PARTICLE = CHAR2INT(C,1)
      READ*, PARTICLE 
      IF (DIAG.EQ.'YES') WRITE(40,*) PARTICLE 
      IF (DIAG.EQ.'YES') WRITE(*,*) PARTICLE
      IF (PARTICLE.LT.0 .OR. PARTICLE.GT.37) THEN
         CALL CLS   
         PRINT*,   '      Entry must be 0 to 37'
         GOTO 20003
      ENDIF
      IF (PARTICLE.EQ.0) PARTICLE=38
      IF (DIAG.EQ.'YES') WRITE(40,*) PARTICLE 
      IF (DIAG.EQ.'YES') WRITE(*,*) PARTICLE
      PRINT*, ' '
!      CALL CLS
20004 WRITE(*,*)'      SELECT DOSE TYPE'

      PRINT*,   '      <1> Secondary Particle Flux (Any rad but TOTAL)'
      PRINT*,   '      <2> ICRP PUB 103 EFFECTIVE DOSE'
      PRINT*,   '      <3> ICRP PUB 60 EFFECTIVE DOSE'
      PRINT*,   '      <4> ICRU H*(10) AMBIENT DOSE EQUIVALENT'
      PRINT*,   '      <5> WHOLE BODY ABSORBED DOSE'
      READ*,DT
      IF (DT.LT.1 .OR. DT.GT.5) THEN
         CALL CLS   
         PRINT*,   '      Entry must be 1 to 5 '
         GOTO 20004
      ENDIF
      IF (DIAG.EQ.'YES') WRITE(40,*) 'dose type is', DT 
      IF (DIAG.EQ.'YES') WRITE(*,*) 'dose type is', DT 
      DK=38*(DT-1)+PARTICLE

20006 PRINT*, ' '
!      WRITE(*,*)'      SELECT COSMIC RAY MODEL              '
!      PRINT*,   '      <1> GCR: ISO TS15390:2004-MSU-NYMMIK '
!      PRINT*,   '      <2> GCR: BADHWAR-ONEILL 2011         '
!      PRINT*,   '      <3> USE DEFAULT         '
!      PRINT*,   '      <4> GCR: HP MODULATED ISO          '
!      PRINT*,   '      <5> SPE: LaRC SEP 1989 EVENT TOTAL '
!      PRINT*,   '      <6> SPE: LaRC FEB 1956 EVENT TOTAL '
!      PRINT*,   '      <7> USER PROVIDED: GCR_MODELS\MY_MODEL.OUT'
!      READ(*,*) GCR
      GCR=4 !ONLY CHOICE IN CARI-7
!      IF (GCR.LT.1 .OR. GCR.GT.7) THEN
!         CALL CLS   
!           PRINT*,   '      Entry must be 1 to 7 '
!           GOTO 20006
!      ENDIF
!      WRITE(*,*)'      USE THE SUPERPOSITION APPROXIMATION <Y,N> ?'
!      PRINT*,   '         (PRIMARY NUCLEI ARE TREATED AS COLLECTIONS'
!      PRINT*,   '          OF N-Z FREE NEUTRONS AND Z FREE PROTONS.)'
!      WRITE(*,*)'                   (DEFAULT IS N)'
!      CALL READKB(C,1)
!      IF (C.EQ.'Y') then
!         SP=1
!      ELSE
         SP=0
!      ENDIF

!      IF (DIAG.EQ.'YES') THEN 
!         WRITE(40,*) 'gcr model is #',GCR,' superposition is', SP 
!      ENDIF
!
! STEP 2, CALCULATE SELECTED DOSE RATE
! YMD IS YYYY/MM/DD
      IF (DIAG.EQ.'YES') WRITE(40,*)'Converting ymdh data to T' 
      CALL YMDH2T(YEAR,MONTH,DAY,HOUR,T)
         IF(MONTH.EQ.0) then 
            YEARLYAVE=.TRUE.
         ELSE
            YEARLYAVE=.FALSE.
         ENDIF
         IF(DAY.EQ.0) then 
            MONTHLYAVE=.TRUE.
         ELSE
            MONTHLYAVE=.FALSE.
         ENDIF
         IF(HOUR.EQ.0) then 
            DAILYAVE=.TRUE.
         ELSE
            DAILYAVE=.FALSE.
         ENDIF

      IF (DIAG.EQ.'YES') THEN 
      WRITE(40,*)'SENDING data to DOSE', LAT, LON, G, T, DK, GCR, SP
         DOSERATE = DOSE(LAT,LON,G,T,DK,GCR,0,SP)
      ELSE
         DOSERATE = DOSE(LAT,LON,G,T,DK,GCR,0,SP)
      ENDIF  
! STEP 3, REPORT DOSE RATE
      DSTR(1) = 'ICRP PUB 103 EFFECTIVE DOSE     '
      DSTR(2) = 'ICRP PUB 60 EFFECTIVE DOSE      '
      DSTR(3) = 'AMBIENT DOSE EQUIVALENT (h*(10))'
      DSTR(4) = 'WHOLE BODY ABSORBED DOSE        '
      DSTR(5) = 'MICROSIEVERTS/HOUR              '
      DSTR(6) = 'MICROGRAY/HOUR                  '
      DSTR(7) = 'Particles per sq-cm per hour    '
      DSTR(8) = 'SECONDARY PARTICLE FLUX         '  
      CALL CLS
      WRITE(*,*) ' '
      WRITE(*,*) 'DATE                     ',T
      WRITE(*,*) 'HELIOCENTRIC POTENTIAL   ',HP, ' MV'
      WRITE(*,*) 'LATITUDE                 ',LAT, NS
      WRITE(*,*) 'LONGITUDE                ',LON, EW
      WRITE(*,*) 'ALTITUDE                 ',G,' g/cm^2'
      WRITE(*,*) 'VERTICAL CUTOFF RIGIDITY ', VC, ' GV'
      WRITE(*,*) 'PARTICLE                 ', PARTSTR(PARTICLE)

      IF (DT.EQ.1) WRITE(*,*) DSTR(8),DOSERATE,PARTSTR(PARTICLE),DSTR(7)   
      IF (DT.EQ.2) WRITE(*,*) DSTR(1),DOSERATE, DSTR(5)   
      IF (DT.EQ.3) WRITE(*,*) DSTR(2),DOSERATE, DSTR(5)   
      IF (DT.EQ.4) WRITE(*,*) DSTR(3),DOSERATE, DSTR(5)   
      IF (DT.EQ.5) WRITE(*,*) DSTR(4),DOSERATE, DSTR(6)   
      IF (DT.EQ.1) WRITE(40,*) DSTR(8),DOSERATE,PARTSTR(PARTICLE),      &
     &   DSTR(7)   
      IF (DT.EQ.2) WRITE(40,*) DSTR(1),DOSERATE, DSTR(5)   
      IF (DT.EQ.3) WRITE(40,*) DSTR(2),DOSERATE, DSTR(5)   
      IF (DT.EQ.4) WRITE(40,*) DSTR(3),DOSERATE, DSTR(5)   
      IF (DT.EQ.5) WRITE(40,*) DSTR(4),DOSERATE, DSTR(6)   

      CALL OOPS(' ',1)
      
      END SUBROUTINE ONESPOT
!                                                                      7
!----6-----------------------------------------------------------------2
        
21000 SUBROUTINE RUN_LOCATIONS
! RUN ALL LOCATIONS IN A DATABASE
      WRITE(40,*) 'CALLING READIT'
      CALL READIT
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE READIT 
      IMPLICIT NONE
      CHARACTER(72)::LOI
      CHARACTER(42)::DRKSTR
      CHARACTER(39)::DKSTR
      CHARACTER(10)::YMD,PARTSTR,YMD1
      CHARACTER(6)::CHRIN,ICAO
      CHARACTER(1)::HPSTR
      CHARACTER(1)::EW,NS,FG,C1
      CHARACTER(30)::PORTNAME,CITY,INFILENAME,OUTFILENAME
      REAL(8)::LAT,LON,ALTITUDE,DEPTH,PALT,DOSE,DOSERATE,VC
      REAL(8)::CHAR2REAL,NSM,EWM,T,SIGMA,KM2DEPTH
      REAL::STARTED,FINISHED,THISTIME,TODOSE,FROMDOSE,INDOSE
      REAL::GOTLAT, GOTLON,GOTALT,GOTDATE,GOTTIME
      INTEGER(4)::Y,M,D,CHAR2INT,HP,OLDY,OLDM,OLDD
      INTEGER(4)::I,HOUR,RAD,DT,SS(11),J,PARTICLE,DOSEKIND
      INTEGER(4)::TESTRUN,GCR,SP,DK,WP,HOURIN
      LOGICAL::NOGOOD
      CHARACTER(3)::DIAGNOG
            
         CHARACTER(12)::VIEWER
         CHARACTER(5)::OS
         CHARACTER(4)::OUTPUT
         CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
                   
      LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE,WITHVAR
      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /ISOSTATS/SIGMA
      COMMON /VCUT/VC
        
      DIAGNOG='no!'
      IF (DIAGNOSE.NE.'YES') TESTRUN=0 
      IF (MENUS.NE.'NO!') THEN
        CALL PICKLOC(INFILENAME)
!        INFILENAME='PLACES.LOC'
      ELSE
        OPEN(UNIT=45,FILE='DEFAULT.INP',STATUS='OLD')
        READ(45,'(A10)') YMD1
        READ(45,*) HOURIN
        READ(45,*) PARTICLE
        READ(45,*) DOSEKIND
        READ(45,'(A30)') INFILENAME
        CLOSE(45)
      ENDIF
      WP=SCAN(INFILENAME,'.',BACK=.TRUE.)
      OUTFILENAME(WP:WP+4)='.ANS'
      OUTFILENAME(1:WP-1)=INFILENAME(1:WP-1)
      DO I=12,(WP+5),-1
         !pad unused characters in OUTFILENAME with blanks to avoid  
         !character spillage from older runs, KC 20150217
         OUTFILENAME(I:I)=" " 
      ENDDO
      OPEN (UNIT=18, FILE=INFILENAME, STATUS='OLD')
      OPEN (UNIT=17, FILE=OUTFILENAME, STATUS='UNKNOWN')
        
        ! UNIT 18 IS EQUIVALENT TO OLD DATAIN, 17 OLD DATAOUT
      C1=','
      WRITE(17,21025)
      IF (DIAGNOSE.EQ.'YES')                                            & 
     &WRITE(40,*) 'CALLING READIT TO READ LIST OF LOCATIONS'  
      IF (DIAGNOG.EQ.'YES') WRITE(*,*) 'READING LIST OF LOCATIONS'  
      CHRIN='------' ! INITIALIZE CHRIN
      OLDY=1111
      OLDM=13
      OLDD=32
      I=0 
      DO 
         READ(18,21100,ERR=21200,END=21200,EOR=21010,ADVANCE='NO')LOI
21010    CHRIN=LOI(1:6) 
         I=I+1
         IF (DIAGNOG.EQ.'YES') WRITE(40,*) 'TESTING LINE',I
         IF (DIAGNOG.EQ.'YES') WRITE(40,*) LOI(1:72)
         IF (CHRIN.EQ.'START-') EXIT 
      ENDDO
      
      DO J = 1,9
         SS(J)=0
      ENDDO     

      DO
         CALL cpu_time(Started)!(Finished)
         I=I+1
         WRITE(*,*) " Analyzing data line: ", I
         READ(18,21100,ERR=21200,END=21200,EOR=21020,ADVANCE='NO')LOI
21020    CHRIN=LOI(1:6)
         IF (DIAGNOG.EQ.'YES') WRITE(40,*)CHRIN//LOI(7:72) 
         IF (CHRIN.EQ.'STOP--') EXIT 
         IF (LOI(1:1).EQ.'!') THEN
            IF (DIAGNOG.EQ.'YES') WRITE(40,*) 'SKIPPING COMMENT LINE'  
         ELSEIF ((LOI(1:1).EQ.'N') .OR. (LOI(1:1).EQ.'S')) THEN !FULL DATA
            IF (DIAGNOG.EQ.'YES') WRITE(40,*) 'READING '//CHRIN//LOI(   &
     &                                                             7:72) 
            CALL FINDCOMMAS(LOI,72,SS)
            HOUR=CHAR2INT(LOI(SCAN(LOI,'H')+1:SCAN(LOI,'H')+2),2)
!            HOUR=CHAR2INT(LOI(SS(7)+1:SS(8)-1),SS(8)-SS(7)-1)
            YMD= LOI(SCAN(LOI,'/')-4:SCAN(LOI,'/')+5)
            CALL DATE2YMD(YMD,Y,M,D)
            CALL YMDH2T(Y,M,D,HOUR,T)
            IF(M.EQ.0) then 
               YEARLYAVE=.TRUE.
            ELSE
               YEARLYAVE=.FALSE.
            ENDIF
            IF(D.EQ.0) then 
               MONTHLYAVE=.TRUE.
            ELSE
               MONTHLYAVE=.FALSE.
            ENDIF
            IF(HOUR.EQ.0) then 
               DAILYAVE=.TRUE.
            ELSE
               DAILYAVE=.FALSE.
            ENDIF
            CALL cpu_time(GOTTIME)!(Finished)
            NS=LOI(1:1)
            IF (NS.EQ.'S') THEN 
               NSM=-1.
            ELSE
               NSM=1.
            ENDIF       
            IF (SCAN(LOI(SS(1)+1:SS(2)-1),'.').NE.0) THEN
               LAT=CHAR2REAL(LOI(SS(1)+1:SS(2)-1),SS(2)-SS(1)-1)
            ELSE
               LAT=REAL(CHAR2INT(LOI(SS(1)+1:SS(2)-1),SS(2)-SS(1)-1)    &
     &                  ,KIND=8)
            ENDIF  
            CALL cpu_time(GOTLAT)!(Finished)
            EW=LOI(SS(3)-1:SS(3)-1)
            IF (SCAN(LOI(SS(3)+1:SS(4)-1),'.').NE.0) THEN
               LON=CHAR2REAL(LOI(SS(3)+1:SS(4)-1),SS(4)-SS(3)-1)
            ELSE
               LON=REAL(CHAR2INT(LOI(SS(3)+1:SS(4)-1),SS(4)-SS(3)-1)    &
     &                 ,KIND=8)
            ENDIF  
            FG=LOI(SS(5)-1:SS(5)-1)
            IF (SCAN(LOI(SS(5)+1:SS(6)-1),'.').NE.0) THEN
               ALTITUDE=CHAR2REAL(LOI(SS(5)+1:SS(6)-1),SS(6)-SS(5)-1)
            ELSE
               ALTITUDE=REAL(CHAR2INT(LOI(SS(5)+1:SS(6)-1),SS(6)-SS(5)  &
     &                    -1),KIND=8)
            ENDIF  
            IF (EW.EQ.'W') THEN 
               EWM=-1.
            ELSE
               EWM=1.
            ENDIF       
            CALL cpu_time(GOTLON)!(Finished)
            RAD = CHAR2INT(LOI(SCAN(LOI,'P')+1:SCAN(LOI,'P')+2),2) 
            IF (RAD.EQ.0) RAD=38
            DT = CHAR2INT(LOI(SCAN(LOI,'D')+1:SCAN(LOI,'D')+1),1)
            GCR = CHAR2INT(LOI(SCAN(LOI,'C')+1:SCAN(LOI,'C')+1),1)
            GCR = 4
            SP = CHAR2INT(LOI(SCAN(LOI,'S')+1:SCAN(LOI,'S')+1),1)
            SP = 0
            IF (DIAGNOSE.EQ.'YES')                                      & 
     &         WRITE(40,*) 'READING LAT AND LON ', NS,LAT,EW,LON
            IF (DIAGNOSE.EQ.'YES')                                      & 
     &         WRITE(40,*) 'READING ALTITUDE ', FG,ALTITUDE
            IF (DIAGNOSE.EQ.'YES')                                      & 
     &         WRITE(40,*) 'READING DATE & HOUR ', YMD, HOUR
            IF (DIAGNOSE.EQ.'YES')                                      & 
     &         WRITE(40,*) 'READING PARTICLE AND UNITS ', RAD, DT
            IF (FG.EQ.'F') THEN 
               IF (DIAGNOG.EQ.'YES') WRITE (40,*) 'CALLING FT2GPCM2'
               CALL FT2GPCMS(ALTITUDE,DEPTH)
            ELSEIF (FG.EQ.'K') THEN
               DEPTH=KM2DEPTH(ALTITUDE) 
            ELSE
               DEPTH=ALTITUDE
            ENDIF
            CALL cpu_time(GOTALT)!(Finished)
            DK=38*DT+RAD-38
!
!       IF (DIAGNOSE.EQ.'YES') THEN 
         WRITE(40,*) ' Sending DOSE latitude = ',NSM*LAT         
         WRITE(40,*) ' Sending DOSE longitude = ',EWM*LON
         WRITE(40,*) ' Sending DOSE altitude = ',DEPTH
         WRITE(40,*) ' Sending DOSE date and time = ',T
         WRITE(40,*) ' Sending DOSE dosekind = ',DK
         WRITE(40,*) ' Sending DOSE gcr = ',GCR
         WRITE(40,*) ' Sending DOSE test = ',TESTRUN
         WRITE(40,*) ' Sending DOSE superposition = ',SP
!       ENDIF
         CALL cpu_time(TODOSE)!(Finished)
         DOSERATE=DOSE(NSM*LAT,EWM*LON,DEPTH,T,DK,GCR,TESTRUN,SP)
         CALL cpu_time(FROMDOSE)!(Finished)
         INDOSE=FROMDOSE-TODOSE
         IF (DIAGNOSE.EQ.'YES')                                         & 
     &   WRITE(40,*)'DOSE RATE AT NORMAL LOCATION IS',DOSERATE,'+/-',   &
     &         SIGMA,DKSTR(DT)
         WRITE(17,21030)NSM*LAT,C1,EWM*LON,C1,REAL(ALTITUDE,KIND=8),C1, &
     &         FG,C1,YMD,C1,HOUR,C1,VC,C1,PARTSTR(RAD),C1,DOSERATE,C1,  &
     &         DRKSTR(DT)
!         WRITE(*,21030)NSM*LAT,C1,EWM*LON,C1,REAL(ALTITUDE,KIND=8),C1,  &
!     &         FG,C1,YMD,C1,HOUR,C1,VC,C1,PARTSTR(RAD),C1,DOSERATE,C1,  &
!     &         DRKSTR(DT)

       ELSEIF (LOI(1:1).EQ.'A' .OR. LOI(1:1).EQ.'a') THEN ! BY CODE
        IF (DIAGNOG.EQ.'YES') WRITE(40,*) 'READING '//CHRIN//LOI(7:66) 
           CALL FINDCOMMAS(LOI,72,SS)
           ICAO= ADJUSTL(TRIM(LOI(SS(1)+1:SS(2)-1)))
           YMD= LOI(SCAN(LOI,'/')-4:SCAN(LOI,'/')+5) 
           HOUR=CHAR2INT(LOI(SCAN(LOI,'H')+1:SCAN(LOI,'H')+2),2)
!          HOUR=CHAR2INT(LOI(SS(5)+1:SS(6)-1),SS(5)-SS(6)-1)

           CALL DATE2YMD(YMD,Y,M,D)
           CALL YMDH2T(Y,M,D,HOUR,T)
            IF(M.EQ.0) then 
               YEARLYAVE=.TRUE.
            ELSE
               YEARLYAVE=.FALSE.
            ENDIF
            IF(D.EQ.0) then 
               MONTHLYAVE=.TRUE.
            ELSE
               MONTHLYAVE=.FALSE.
            ENDIF
            IF(HOUR.EQ.0) then 
               DAILYAVE=.TRUE.
            ELSE
               DAILYAVE=.FALSE.
            ENDIF
           CALL cpu_time(GOTTIME)!(Finished)

!         SP = CHAR2INT(LOI(SCAN(LOI,'S',BACK=.TRUE.)+1:                 &
!     &                      SCAN(LOI,'S',BACK=.TRUE.)+1),1)
!         GCR = CHAR2INT(LOI(SCAN(LOI,'C',BACK=.TRUE.)+1:                &
!     &                      SCAN(LOI,'C',BACK=.TRUE.)+1),1)
          GCR=4
          SP=0
         RAD = CHAR2INT(LOI(SCAN(LOI,'P',BACK=.TRUE.)+1:                &
     &                      SCAN(LOI,'P',BACK=.TRUE.)+2),2)
         IF (RAD.EQ.0) RAD=38 
         DT = CHAR2INT(LOI(SCAN(LOI,'D',BACK=.TRUE.)+1:                 &
     &                     SCAN(LOI,'D',BACK=.TRUE.)+1),1)
         FG=LOI(SS(3)-1:SS(3)-1)
         IF (SCAN(LOI(SS(3)+1:SS(4)-1),'.').NE.0) THEN
            ALTITUDE=CHAR2REAL(LOI(SS(3)+1:SS(4)-1),SS(4)-SS(3)-1)
         ELSE
            ALTITUDE=REAL(CHAR2INT(LOI(SS(3)+1:SS(6)-1),SS(4)-SS(3)     &
     &                    -1),KIND=8)
            ENDIF  

         CALL PORT_INFO(ICAO,PORTNAME,CITY,NS,LAT,EW,LON,PALT,NOGOOD)
         IF (DIAGNOSE.EQ.'YES') THEN 
          WRITE(40,*) 'READING CODE ', ICAO
          WRITE(40,*) 'READING LAT AND LON ', NS,LAT,EW,LON
          WRITE(40,*) 'READING ALTITUDE ', FG,ALTITUDE
          WRITE(40,*) 'READING DATE & HOUR ', YMD, HOUR
          WRITE(40,*) 'READING PARTICLE AND UNITS ', RAD, DT
          WRITE(40,*) 'READING gcr = ',GCR
          WRITE(40,*) 'READING test = ',TESTRUN
          WRITE(40,*) 'READING superposition = ',SP
         ENDIF
         IF (NS.EQ.'S') THEN 
            NSM=-1.
         ELSE
            NSM=1.
         ENDIF       
         CALL cpu_time(GOTLAT)!(Finished)
         IF (EW.EQ.'W') THEN 
            EWM=-1.
         ELSE
            EWM=1.
         ENDIF       
         CALL cpu_time(GOTLON)!(Finished)
         IF (FG.EQ.'F') THEN 
            IF (DIAGNOG.EQ.'YES') WRITE (40,*) 'CALLING FT2GPCM2'
            CALL FT2GPCMS(ALTITUDE,DEPTH)
         ELSEIF (FG.EQ.'K') THEN
            DEPTH=KM2DEPTH(ALTITUDE) 
         ELSE
            DEPTH=ALTITUDE
         ENDIF
         CALL cpu_time(GOTALT)!(Finished)
         DK=38*DT+RAD-38

!         IF (DIAGNOSE.EQ.'YES') THEN 
          WRITE(40,*) ' Sending DOSE latitude = ',NSM*LAT         
          WRITE(40,*) ' Sending DOSE longitude = ',EWM*LON
          WRITE(40,*) ' Sending DOSE altitude = ',DEPTH
          WRITE(40,*) ' Sending DOSE date and time = ',T
          WRITE(40,*) ' Sending DOSE dosekind = ',DK
          WRITE(40,*) ' Sending DOSE gcr = ',GCR
          WRITE(40,*) ' Sending DOSE test = ',TESTRUN
          WRITE(40,*) ' Sending DOSE superposition = ',SP
!         ENDIF

         IF(M.EQ.0) then 
            YEARLYAVE=.TRUE.
         ELSE
            YEARLYAVE=.FALSE.
         ENDIF
         IF(D.EQ.0) then 
            MONTHLYAVE=.TRUE.
         ELSE
            MONTHLYAVE=.FALSE.
         ENDIF
         IF(HOUR.EQ.0) then 
            DAILYAVE=.TRUE.
         ELSE
            DAILYAVE=.FALSE.
         ENDIF

         CALL cpu_time(TODOSE)!(Finished)
         DOSERATE=DOSE(NSM*LAT,EWM*LON,DEPTH,T,DK,GCR,TESTRUN,SP)
         CALL cpu_time(FROMDOSE)!(Finished)
         INDOSE=FROMDOSE-TODOSE

         IF (DIAGNOSE.EQ.'YES') THEN 
            WRITE(40,*)'DOSE RATE AT', ICAO,' IS ', DOSERATE,'+/-',     &
     &      SIGMA, DKSTR(DT)
         ENDIF
         WRITE(17,21030)NSM*LAT,C1,EWM*LON,C1,REAL(ALTITUDE,KIND=8),C1, &
     &  FG,C1,YMD,C1,HOUR,C1,VC,C1,PARTSTR(RAD),C1,DOSERATE,C1,         &
     &  DRKSTR(DT)
       ELSE
         WRITE(40,*)'COMMENT OR IMPROPER FORMAT FOR DATA AT LINE: ', I, &
     &              'IN PLACES.DAT'
         WRITE(40,*) LOI(1:60)
         WRITE(17,*)'COMMENT OR IMPROPER FORMAT FOR DATA AT LINE: ', I, &
     &              'IN PLACES.DAT'
         WRITE(17,*) LOI(1:60)
       ENDIF 
       WRITE(40,*)'TIME BEFORE CALLING DOSE WAS ', TODOSE-STARTED
       WRITE(40,*)'  TIME FOR DATE WAS ',GOTTIME-STARTED
       WRITE(40,*)'  TIME FOR LAT WAS ',GOTLAT-GOTTIME
       WRITE(40,*)'  TIME FOR LON WAS ',GOTLON-GOTLAT
       WRITE(40,*)'  TIME FOR ALT WAS ',GOTALT-GOTLON
       WRITE(40,*)'TIME IN DOSE WAS ',INDOSE

       CALL CPU_TIME(FINISHED)
       THISTIME=FINISHED-STARTED
       WRITE(40,*)'TIME AFTER DOSE WAS ',FINISHED-FROMDOSE
       WRITE(40,*)'TIME ON THIS LINE DATA WAS ',THISTIME
       WRITE(40,*)' '
      ENDDO
      CLOSE(17)
      CLOSE(18)
      RETURN
        ! TABLE HEADER
21025 FORMAT('     LAT,       LON,     ALTITUDE,    DATE,    HR, ',     &
     &   'VCR(GV), PARTICLE,  DOSE RATE,        UNIT,      ',           &
     &   'QUANTITY ')   
        ! TABLE CONTENTS
21030 FORMAT(F10.5,A1,F10.5,A1,F11.4,A1,A1,A1,A10,A1,I4,A1,F6.2,A1,A10, &
     &       A1,ES11.4,A1,A45)        
21100 FORMAT(A66)
21200 WRITE(40,*) 'CANNOT READ LOCATIONS, FILE IS CORRUPT OR EMPTY'   
      CALL OOPS('CANNOT READ LOCATIONS, FILE IS CORRUPT OR EMPTY',49)
      WRITE(17,*) 'CANNOT READ LOCATIONS, FILE IS CORRUPT OR EMPTY'   
      CLOSE(17)
      CLOSE(18)
      RETURN
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2    
      SUBROUTINE FINDCOMMAS(SI,L,COMMAS)
         INTEGER(4)::L,COMMAS(11),RM,LM
         CHARACTER(L)::SI
         CHARACTER(3)::DIAGNOSE

         DIAGNOSE='NO!'
         
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'FINDING COMMAS IN'
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) SI
      RM=SCAN(SI,',',BACK=.TRUE.)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'RM ', RM
      I=0
      LM=1 
      DO
         I=I+1
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'SCANNING FOR COMMA', I
         COMMAS(I)= SCAN(SI(LM+1:RM+1),',')
         LM=COMMAS(I)+LM
         COMMAS(I)=LM
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'LM ', LM
         IF (LM.EQ.RM .OR. I.GT.9) EXIT
      ENDDO
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'COMMAS AT'
      DO J=1,I
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) COMMAS(J)
      ENDDO
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
        
22000 SUBROUTINE NEWLOCS
        INTEGER(4)::LOF
        CHARACTER(52)::FILENAME
! View or revise the locations database.    
        WRITE(40,*) 'CALLING NEWLOCS'
        CALL PICKLOC(FILENAME)
        LOF=LEN_TRIM(FILENAME)
!        CALL SHOWPICK('PLACES.LOC',10)
        CALL SHOWPICK(FILENAME,LOF)
        CALL LOCATIONS
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
!30000-39999
! BEGIN SUBS AND FUNCTIONS FOR AIRPORT DATA MANIPULATION
!  A. LOOK UP AN EXISTING AIRPORT GIVEN A CITY NAME
!  B. LOOK UP AN EXISTING AIRPORT GIVEN AN AIRPORT NAME        
!  C. LOOK UP AN EXISTING AIRPORT GIVEN AN ICAO CODE
30000 SUBROUTINE PORT_INFO(ICAO_IN,APORT,CITY,NS,LAT,EW,LON,ALT,FRACK)
       
         INTEGER(4)::I,N,K,CHAR2INT
         INTEGER(4)::ALTF,LOND,LONM,LONS,LATD,LATM,LATS
          REAL(8)::LON,LAT,ALT
         CHARACTER(6), INTENT(IN)::ICAO_IN
         CHARACTER(6)::ICAO,TESTCODE 
           CHARACTER(3)::ALTCODE
         CHARACTER(1)::LONEW,LATNS,NSS,EWS,NS,EW
           CHARACTER(30)::APORT,CITY
         CHARACTER(91)::PORTINFO
         LOGICAL::FRACK

           CHARACTER(3)::DIAGNOSE

      DIAGNOSE='NO!'
      FRACK=.FALSE. !assume record is in the DataBases
         TESTCODE=ADJUSTL(ICAO_IN)
          ! FIND THE RECORD
      OPEN (UNIT=98, FILE='AIRPORTS/PORT.NDX',STATUS='OLD')
      DO 
         READ(98,FMT=30008,ERR=30007,END=30007)PORTINFO
         ICAO=PORTINFO(31:36)
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) TESTCODE, ICAO
         IF (LEN_TRIM(ICAO).EQ.4) THEN
           IF (TESTCODE(1:4).EQ.ICAO(1:4)) THEN 
            IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'MATCHED ',ICAO_IN, ICAO
            EXIT
           ENDIF  
         ENDIF
         IF (LEN_TRIM(ICAO).EQ.5) THEN
           IF (TESTCODE(1:5).EQ.ICAO(1:5)) THEN
            IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'MATCHED ',ICAO_IN, ICAO
            EXIT
           ENDIF  
         ENDIF
         IF (LEN_TRIM(ICAO).EQ.6) THEN
           IF (TESTCODE(1:6).EQ.ICAO(1:6)) THEN
            IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'MATCHED ',ICAO_IN, ICAO
            EXIT
           ENDIF  
         ENDIF
      ENDDO
      APORT=PORTINFO(1:30)
      CITY=PORTINFO(37:66)
      LATNS=PORTINFO(74:74)
      LATD=CHAR2INT(PORTINFO(67:68),2) 
      LATM=CHAR2INT(PORTINFO(69:70),2) 
      LATS=CHAR2INT(PORTINFO(71:72),2) 
      LONEW=PORTINFO(83:83)
      LOND=CHAR2INT(PORTINFO(75:77),3)
      LONM=CHAR2INT(PORTINFO(78:79),2) 
      LONS=CHAR2INT(PORTINFO(80:81),2) 
      ALTF=CHAR2INT(PORTINFO(84:89),6) 
      LON=REAL(LOND,KIND=8)+(LONM/60.)+(LONS/3600.) 
      LAT=REAL(LATD,KIND=8)+(LATM/60.)+(LATS/3600.) 
      ALT=REAL(ALTF,KIND=8)/10.
      NS=NSS(LATNS)
      EW=EWS(LONEW)
      CLOSE(98)
      IF (DIAGNOSE.EQ.'YES') THEN 
         WRITE(40,30001) ICAO_IN 
         WRITE(40,30002) APORT
         WRITE(40,30003) CITY
         WRITE(40,30004) LATD, LATM, LATS, NS
         WRITE(40,30005) LOND, LONM, LONS, EW
         WRITE(40,30006) REAL(ALTF,KIND=8)/10.0
      ENDIF
30001 FORMAT(10X,'ICAO CODE: ',A6, '                OTHER CODE: ',A3)
30002 FORMAT(10X,'AIRPORT NAME: ',A30)
30003 FORMAT(10X,'CITY NAME: ',A30)
30004 FORMAT(10X,'LATITUDE: ',I2,' DEGS ',I2,' MINS ',I2,' SECS ',A1)
30005 FORMAT(10X,'LONGITUDE: ',I3,' DEGS ',I2,' MINS ',I2,' SECS ',A1)
30006 FORMAT(10X,'ALTITUDE: ',F7.1, ' FEET')
      RETURN   
30007 CLOSE(98)
      FRACK=.TRUE. !could not match the codes
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'PORT_INFO failed for ', ICAO_IN 
      WRITE(*,*)'PORT_INFO failed for ', ICAO_IN 
      CALL OOPS ('ICAO CODE NOT FOUND IN DATABASES',32)
30008 FORMAT(A91)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE MENU_MAIN
      CHARACTER(1)::CHOICE
70000 CALL CLS                                       ! BEGIN INPUT LOOP
      CALL MENUHEADER
      WRITE(*,70004)'                 MAIN MENU                        '
      PRINT*, ' '
      WRITE(*,70004)'<1>  HELP file (Read me).                         '
      PRINT*, ' '
      WRITE(*,70004)'<2>  Galactic cosmic radiation received on flights.&      
     &'
      PRINT*, ' '
      WRITE(*,70004)'<3>  View, add, or change airport information.    '
      PRINT*, ' '
      WRITE(*,70004)'<4>  Galactic cosmic radiation at user-specified  '
      WRITE(*,70004)'     altitude and geographic coordinates.         '
      PRINT*, ' '
      WRITE(*,70004)'<5>  View or update heliocentric potentials,      '
      WRITE(*,70004)'     Kp indices, or Forbush effect data.          '
      PRINT*, ' '
      WRITE(*,70004)'<6>  Change output settings. View old results.    '
      PRINT*, ' '
      WRITE(*,70004)'<7>  Exit program.                                '
      PRINT*, ' '
      WRITE(*,70001)'.'
      CALL READKB(CHOICE,1)
      IF (CHOICE.EQ.'q'.OR.CHOICE.EQ.'Q') CHOICE='7'
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'4'.OR.CHOICE.EQ.'5'.OR.CHOICE.EQ.'6'.OR.           &
     &    CHOICE.EQ.'7') THEN
         GOTO 70002
      ELSE 
         GOTO 70000 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
70002 CONTINUE
! DIAGNOSTIC PRINT*, CHOICE  
      IF (CHOICE.EQ.'1') CALL SHOWHELP
      IF (CHOICE.EQ.'2') CALL FLIGHTS
      IF (CHOICE.EQ.'3') CALL AIRPORTS
      IF (CHOICE.EQ.'4') CALL LOCATIONS
      IF (CHOICE.EQ.'5') CALL SOLARMOD
      IF (CHOICE.EQ.'6') CALL OUTPUTS
      IF (CHOICE.EQ.'7') STOP
      GOTO 70000 
70001 FORMAT(10X,'Type 1, 2, 3, 4, 5, 6, or 7 and press <ENTER> ',A1)
70004 FORMAT(10X,A50)
      END SUBROUTINE MENU_MAIN  
!                                                                      7
!----6-----------------------------------------------------------------2
! 71000
      SUBROUTINE SHOWHELP
      CHARACTER(8)::STARTHELP
      CALL CLS
      STARTHELP = 'HELP.TXT'
      CALL SHOWPICK(STARTHELP,8)
      END SUBROUTINE SHOWHELP
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE SHOWKP
      CHARACTER(21)::FILENAME
      CALL CLS
      CLOSE(21)
      FILENAME = 'KP_INDEX/KP_INDEX.TXT'
      CALL SHOWPICK(FILENAME,21)
      END SUBROUTINE SHOWKP
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE SHOWFORBUSH
      CHARACTER(19)::FILENAME
      CHARACTER(1)::CHOICE
      CALL CLS
      FILENAME = 'FORBUSH/FORBUSH.DAT'
      CLOSE(131)
      CALL SHOWPICK(FILENAME,19)
      PRINT*, 'Press any character and then <Enter>'
      CALL READKB(CHOICE,1)
      CALL LOADFORBUSH
      END SUBROUTINE SHOWFORBUSH
!                                                                      7
!----6-----------------------------------------------------------------2
! 72000
      SUBROUTINE FLIGHTS
      CALL MENU_FLIGHT
      END SUBROUTINE
!----6-----------------------------------------------------------------2
! 72000
      SUBROUTINE LOCATIONS
      CALL MENU_LOCATIONS
      END SUBROUTINE!                                                                      7
!----6-----------------------------------------------------------------2
! 72000
      SUBROUTINE MENU_FLIGHT
      CHARACTER(1)::CHOICE
72000 CALL CLS
      CALL MENUHEADER
      WRITE(*,72004)'                FLIGHT MENU                       '
      PRINT*, ' '
      WRITE(*,72004)'Indicate which kind of flight data you are using. '
      PRINT*, ' '
      WRITE(*,72004)'<1>  Flight path(s) nearly geodesic (i.e., great  '
      WRITE(*,72004)'     circle) route(s) between airports.           '
      WRITE(*,72004)'     (*.BIG files, same as CARI-6, -6P, and -6W). '
      PRINT*, ' '
      WRITE(*,72004)'<2>  Flight path(s) defined by many waypoints     '
      WRITE(*,72004)'     (*.DEG files, same as CARI-6M and -6PM).     '
      PRINT*, ' '
      WRITE(*,72004)'<3>  Return to Main Menu.                         '
      PRINT*, ' '
      WRITE(*,72004)'<4>  Exit program.                                '
      PRINT*,''
      PRINT*,''
      PRINT*,''
      WRITE(*,72001) '.' 
      CALL READKB(CHOICE,1)
72001 FORMAT(10X,'Type 1, 2, 3, or 4 and press <ENTER> ',A1)
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'4') THEN
         GOTO 72002
      ELSE 
         GOTO 72000 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
72002 CONTINUE
! DIAGNOSTIC PRINT   PRINT*, "SUCCESFUL DATA ENTRY ", CHOICE
      IF (CHOICE.EQ.'1') CALL RUNBIG                
      IF (CHOICE.EQ.'2') CALL RUNDEG            
      IF (CHOICE.EQ.'4') STOP 
      CALL MENU_MAIN
72003 FORMAT (A1)
72004 FORMAT (10X,A50)
72005 FORMAT (10X,A9)

      END SUBROUTINE MENU_FLIGHT
!                                                                      7
!----6-----------------------------------------------------------------2
! 73000
      SUBROUTINE AIRPORTS(OS)
      CHARACTER(5)::OS
      CALL MENU_AIRPORT(OS)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
!
      SUBROUTINE MENU_AIRPORT(OS)

      CHARACTER(1)::CHOICE
      CHARACTER(5)::OS

73000 CALL CLS      
      CALL MENUHEADER
      WRITE(*,73004)'                AIRPORT MENU                      '
      PRINT*, ' '
      WRITE(*,73004)'<1>  Find an airport by airport name.             '
      PRINT*, ' '
      WRITE(*,73004)'<2>  Find airport by city.                        '
      PRINT*, ' '
      WRITE(*,73004)'<3>  Find airport by code.                        '
      PRINT*, ' '
      WRITE(*,73004)'<4>  Add new airport or replace existing airport. '
      PRINT*, ' '
      WRITE(*,73004)'<5>  Associate a non-ICAO code with an airport.   '
      PRINT*, ' '
      WRITE(*,73004)'<6>  Return to Main Menu.                         '
      PRINT*, ' '
      WRITE(*,73004)'<7>  Exit program.                                '
      PRINT*,''
      PRINT*,''
      PRINT*,''
      WRITE(*,73001) '.' 
      CALL READKB(CHOICE,1)
73001 FORMAT(10X,'Type 1, 2, 3, 4, 5, 6, or 7 and press <ENTER> ',A1)
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'4'.OR.CHOICE.EQ.'5'.OR.CHOICE.EQ.'6'.OR.           &
     &    CHOICE.EQ.'7') THEN
         GOTO 73002
      ELSE 
         GOTO 73000 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
73002 CONTINUE
! DIAGNOSTIC PRINT   PRINT*, "SUCCESFUL DATA ENTRY ", CHOICE
      IF (CHOICE.EQ.'1') CALL FINDPORT(1)                !73100
      IF (CHOICE.EQ.'2') CALL FINDPORT(2)                !73100 
      IF (CHOICE.EQ.'3') CALL FINDPORT(3)                !73100 
      IF (CHOICE.EQ.'4') CALL ADDAPORT                   !73400
      IF (CHOICE.EQ.'5') CALL ADDACODE                   !73500
      IF (CHOICE.EQ.'7') STOP 
      CALL MENU_MAIN
73003 FORMAT (A1)
73004 FORMAT (10X,A50)
73005 FORMAT (10X,A9)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE MAKE_NDX(L,MXSIZ,OS)
! THIS SUB CREATES THE PORT.NDX FILE FROM THE AIRPORT DATABASES        
! THIS SUB CREATES THE CITY.NDX FILE FROM THE AIRPORT DATABASES        

         INTEGER(4)::I,J,K,L,IREC
         INTEGER(4)::MXSIZ
         CHARACTER(L)::PORTINFO(MXSIZ),PI2(MXSIZ),TESTPORT
         INTEGER(4)::PLATD, PLATM, PLATS, PALTF
         INTEGER(4)::PLOND, PLONM, PLONS
         INTEGER(4)::CHAR2INT
         CHARACTER(1)::PLONEW,PLATNS
         CHARACTER(3)::OS
         CHARACTER(6)::ICAO
         CHARACTER(6)::ICAO_IN 
         CHARACTER(30)::PNAME,CNAME
         LOGICAL::LEXIST
         CHARACTER(3)::DIAGNOSE

!           DIAGNOSE='YES'
            DIAGNOSE='NO!'           
!SAMPLE RECORD, 91 CHARACTERS LON 
!PORTNAME 1-30 
!CODE 31-36 
!CITY 37-66 
!LAT 67-73 DDMMSSBB
!N/S 74 ANY NON-NUMBER OR } INDICATES S
!LON 75-82
!E/W 83 ANY NON-NUMBER OR } INDICATES E
!ALT(FT X10)84-89 
!CARRAIGE RETURN AND LINEFEED CHARACTERS 90-91      
!     
      REWIND(27)
      REWIND(32)                   
      J=0
      PRINT*,'Reading primary airport database'
      DO 
         J=J+1
         READ(27,FMT=30103,ERR=30101) TESTPORT
!       IF (DIAGNOSE.EQ.'YES')
         PORTINFO(J)=TESTPORT
           
         WRITE(40,30104) J, PORTINFO(J)
         ICAO=TESTPORT(31:36)
         PALTF=CHAR2INT(TESTPORT(84:88),5)
         PNAME=TESTPORT(1:30)
         CNAME=TESTPORT(37:66)
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) ICAO, PNAME, CNAME, PALTF 
         IF (TRIM(ICAO).EQ.'----'.AND.PALTF.EQ.99999) THEN
            K=J-1
            EXIT
         ENDIF       
         IF (J.EQ.MXSIZ) THEN
            K=J
            EXIT
         ENDIF          
      ENDDO
        GOTO 30102
30101 CALL EPITATH('AIRPORTS.DAT IS CORRUPTED',25)
30102 IREC=K
30103 FORMAT (A91)
30104 FORMAT (I5,2X,A91)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'LAST GOOD RECORD WAS', IREC
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) ICAO, PNAME, CNAME, PALTF 

! INCLUDE USER ADDED AIRPORT LIST
      PRINT*, 'Adding any user entered airports'
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
        INQUIRE (FILE='AIRPORTS\NEWPORTS.DAT',EXIST=LEXIST)
      ELSE
        INQUIRE (FILE='AIRPORTS/NEWPORTS.DAT',EXIST=LEXIST)
      ENDIF
      IF (LEXIST) THEN
         J=0
      ELSE
         CALL EPITATH('NEWPORTS.DAT IS MISSING  ',25)
      ENDIF
      DO 
         J=J+1
         READ(32,FMT=30103,ERR=30115)PORTINFO(IREC+J)
         !IF (DIAGNOSE.EQ.'YES') 
         WRITE(40,30104) J, PORTINFO(IREC+J)
         K=IREC+J
         TESTPORT=PORTINFO(K)
         PI2(K)=TESTPORT
         ICAO=TESTPORT(31:36)
         PALTF=CHAR2INT(TESTPORT(84:88),5)
         PNAME=TESTPORT(1:30)
         CNAME=TESTPORT(37:66)
         IF (TRIM(ICAO).EQ.'----'.AND.PALTF.EQ.99999) THEN
            K=K-1
            EXIT
         ENDIF
      ENDDO
30115 CONTINUE
       
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'LAST GOOD RECORD WAS', K
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) ICAO, PNAME, CNAME, PALTF
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'SORTING ON AIRPORT NAME'
!      CALL SHELLSORT(PORTINFO,L,K+1)
      PRINT*, 'Generating sorted airport lists' 
      CALL BUBBLESORT(PORTINFO,L,K+1)
      OPEN (UNIT=98, FILE='AIRPORTS/PORT.NDX',STATUS='UNKNOWN')
      DO I=1,K
         WRITE(98,FMT=30103) PORTINFO(I)
      ENDDO
      CLOSE(98)
      IF (DIAGNOSE.EQ.'YES') THEN 
         WRITE(40,*) 'FINISHED WRITING PORT.NDX'
         OPEN (UNIT=98, FILE='AIRPORTS/PORT.NDX',STATUS='OLD')
         DO J=100,K-100
         IF (MOD(J,100) .EQ. 0) THEN
         READ(98,FMT=30103,ERR=30116)TESTPORT
         WRITE(40,30104) J, TESTPORT
         ICAO=TESTPORT(31:36)
         PNAME=TESTPORT(1:30)
         CNAME=TESTPORT(37:66)
         PLATNS=TESTPORT(74:74)
         PLATD=CHAR2INT(TESTPORT(67:68),2) 
         PLATM=CHAR2INT(TESTPORT(69:70),2) 
         PLATS=CHAR2INT(TESTPORT(71:72),2) 
         PLONEW=TESTPORT(83:83)
         PLOND=CHAR2INT(TESTPORT(75:77),3)
         PLONM=CHAR2INT(TESTPORT(78:79),2) 
         PLONS=CHAR2INT(TESTPORT(80:81),2) 
         PALTF=CHAR2INT(TESTPORT(84:88),5)
         WRITE(40,*) ICAO, PNAME, CNAME
         ENDIF
          ENDDO
30116    CLOSE(98) 
      ENDIF  
! NOW SWAP NAME AND CITY, THEN RE-SORT ON CITY NAME

      PI2=PORTINFO
      DO I=1,K
         
         PI2(I)(1:30)=PORTINFO(I)(37:66)
         PI2(I)(37:66)=PORTINFO(I)(1:30)
          IF (DIAGNOSE.EQ.'YES') WRITE(40,*) I, PI2(I)
      ENDDO 
      CALL BUBBLESORT(PI2,L,K+1)  
! WRITE PRELIM VERSION OF CITY.NDX
      OPEN (UNIT=98, FILE='AIRPORTS/CITY.NDX',STATUS='UNKNOWN')
      DO I=1,K
           WRITE(98,FMT=30103) PI2(I)
      ENDDO
      CLOSE(98)
      IF (DIAGNOSE.EQ.'YES') THEN 
         WRITE(40,*) 'FINISHED WRITING CITY.NDX'
         OPEN (UNIT=98, FILE='AIRPORTS/CITY.NDX',STATUS='OLD')
         DO J=100,K-100
            READ(98,FMT=30103,ERR=30118)TESTPORT
            IF (MOD(J,100).EQ. 0) THEN
               WRITE(40,30104) J, TESTPORT
               ICAO=TESTPORT(31:36)
               PNAME=TESTPORT(1:30)
               CNAME=TESTPORT(37:66)
               PLATNS=TESTPORT(74:74)
               PLATD=CHAR2INT(TESTPORT(67:68),2) 
               PLATM=CHAR2INT(TESTPORT(69:70),2) 
               PLATS=CHAR2INT(TESTPORT(71:72),2) 
               PLONEW=TESTPORT(83:83)
               PLOND=CHAR2INT(TESTPORT(75:77),3)
               PLONM=CHAR2INT(TESTPORT(78:79),2) 
               PLONS=CHAR2INT(TESTPORT(80:81),2) 
               PALTF=CHAR2INT(TESTPORT(84:88),5)
               WRITE(40,*) ICAO, PNAME, CNAME
            ENDIF
          ENDDO
30118    CLOSE(98) 
      ENDIF         
      END SUBROUTINE MAKE_NDX
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE FINDPORT(I)
 
         INTEGER(4)::I,J,K,IREC
         INTEGER(4), DIMENSION(11)::PLATD, PLATM, PLATS, PALTF
         INTEGER(4), DIMENSION(11)::PLOND, PLONM, PLONS
             INTEGER(4):: CHAR2INT
         CHARACTER(1), DIMENSION(11)::PLONEW,PLATNS
         CHARACTER(6), DIMENSION(11)::ICAO,SICAO
           CHARACTER(6)::ICAO_IN 
         CHARACTER(30), DIMENSION(11)::PNAME,CNAME,SORTED
         CHARACTER(30)::INPUT
         CHARACTER(1)::CHOICE
         CHARACTER(3)::DIAGNOSE
         CHARACTER(17)::FILENAME
         LOGICAL::FRACK
      DIAGNOSE='NO!'
! I=1, SEARCH FOR NEAREST PORTS BY AIRPORT NAME
! I=2, SEARCH FOR NEAREST PORTS BY CITY NAME
! I=3, 

! USER INPUT TARGET STRING

73112 CONTINUE
      CALL CLS
      CALL MENUHEADER
      SELECT CASE(I)
      CASE(1)
        FILENAME = 'AIRPORTS/PORT.NDX'
        CALL SHOWPICK(FILENAME,17)
      CASE(2)
        FILENAME = 'AIRPORTS/CITY.NDX'
        CALL SHOWPICK(FILENAME,17)
      CASE(3)
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,'         Please enter an ICAO code ' 
        CALL READKB(INPUT,6)
        IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'CALLING REPORT BY ICAO'
        CALL REP_BY_ICAO(INPUT,FRACK)
      END SELECT

      PRINT*, ' '
      PRINT*, ' '
      PRINT*, ' '
      PRINT*, ' '
      PRINT*, ' '
      WRITE(*,73114)'<1>  Enter an ICAO code to view a full report.    '
      PRINT*, ' '
      WRITE(*,73114)'<2>  Return to Airport Menu.                      '
      PRINT*, ' '
      WRITE(*,73114)'<3>  Return to Main Menu.                         '
      PRINT*, ' '
      WRITE(*,73114)'<4>  Exit the program.                            '
      PRINT*, ' '
      PRINT*, ' '
      PRINT*, ' '
      PRINT*,''
      WRITE(*,73111) '.' 
73110 CALL READKB(CHOICE,1)
73111 FORMAT(10X//'Type 1, 2, 3, or 4 and press <ENTER> ',A1)
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'4') THEN
         GOTO 73122
      ELSE 
         GOTO 73110 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
73122 CONTINUE
      SELECT CASE (CHAR2INT(CHOICE,1))
      CASE(1)
         CALL CLS
         CALL MENUHEADER
         PRINT*,' '
         PRINT*,'          Enter an ICAO code to view a full report.  '
         CALL READKB(ICAO_IN,6)     
         CALL REP_BY_ICAO(ICAO_IN,FRACK)
      CASE(2)
         CALL MENU_AIRPORT
      CASE(3)
         CALL MENU_MAIN 
      CASE(4)
         STOP
      END SELECT 

73103 FORMAT(5X,'CITY',13X,12X,'AIRPORT',21X,'CODE')
73104 FORMAT(5X,'AIRPORT',12X,13X,'CITY',21X,'CODE')
73105 FORMAT('CODE',13X,'CITY',13X,12X,'AIRPORT')
73106 FORMAT(A30,1X,A30,1X,A6)
73107 FORMAT(A6,1X,A30,1X,A30)
73113 FORMAT('-------------------------------------------------------'  &
     &'---------------')
73114 FORMAT(10X,A50)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! 73200 !A simple bubble sort, fast enough and doesn't corrupt the chars

      SUBROUTINE BUBBLESORT(SORTED, L, N) !-----needs numvals and custlist
       IMPLICIT NONE
       LOGICAL, EXTERNAL :: GTQ !funtion to tell which person goes first
       INTEGER(4), INTENT(IN) :: L 
       INTEGER(4), INTENT(IN) :: N  
       INTEGER(4) :: i,j
       CHARACTER(L):: TEMP1, TEMP2
       CHARACTER(L), INTENT(INOUT), DIMENSION(N) :: SORTED 

       CHARACTER(3)::DIAGNOSE

          DIAGNOSE='NO!'
       IF (DIAGNOSE.EQ.'YES') WRITE (*,*) 'Now sorting airports by name'
       IF (DIAGNOSE.EQ.'YES') THEN
           WRITE (40,*)'Now sorting airports by name'
           WRITE(40,*)'As passed to bubblesort'
        
           DO i = 1,N
             WRITE(40,*) SORTED(i)
           ENDDO
       ENDIF
       DO j=2, N-2
          DO i=1,N-j
              TEMP1=SORTED(i)
              TEMP2=SORTED(i+1)
              IF (GTQ(TEMP1,TEMP2,L)) THEN
                 SORTED(i)=TEMP2
                 SORTED(i+1)=TEMP1
              ENDIF
          ENDDO
       ENDDO    
          
         IF (DIAGNOSE.EQ.'YES') WRITE (*,*) 'End of sort'
         IF (DIAGNOSE.EQ.'YES') THEN 
           WRITE (40,*) 'End of sort'
           WRITE(40,*)'Checking after bubblesort'
           DO i = 1,N
             WRITE(40,*) SORTED(i)
           ENDDO
         ENDIF 
      END SUBROUTINE BUBBLESORT
!----------------------------------------------------------------------
      LOGICAL FUNCTION GTQ (a, b, c) 
! Greater Than Query
! This function takes airport names or codes as the arguments 
! and figures out which sorts out first.
    
         IMPLICIT NONE
         INTEGER(4) :: c
         CHARACTER(c), INTENT(IN)::a,b 
!grab the arguments and format them to make Fortran happy

           GTQ = .FALSE.
!if no other conditions are met then the 2nd AIRPORT comes first
         IF ((a==b) .OR. (LLT ( b , a))) THEN
             GTQ = .TRUE.
           END IF 
      END FUNCTION GTQ    
!                                                                      7
!----6-----------------------------------------------------------------2
! 73300
      SUBROUTINE REP_BY_ICAO(ICAO_IN,FRACK)
       
         INTEGER(4)::I,N,K,CHAR2INT
         INTEGER(4)::ALTF,LOND,LONM,LONS,LATD,LATM,LATS
         CHARACTER(6)::ICAO_IN,ICAO 
         CHARACTER(3)::ALTCODE
         CHARACTER(1)::LONEW,LATNS,NSS,EWS
         CHARACTER(30)::APORT,CITY
         CHARACTER(91)::PORTINFO
         CHARACTER(3)::DIAGNOSE
         LOGICAL::FRACK
         FRACK=.FALSE.
         DIAGNOSE='NO!'
! FIND THE RECORD
      OPEN (UNIT=98, FILE='AIRPORTS/PORT.NDX',STATUS='OLD')
      I=0 
      DO 
         I=I+1
         READ(98,FMT=73308,ERR=73307,END=73307)PORTINFO
         ICAO=PORTINFO(31:36)
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*) ICAO_IN, ICAO
         IF (LEN_TRIM(ICAO_IN).EQ.4) THEN
            IF (ICAO_IN(1:4).EQ.ICAO(1:4)) EXIT
         ENDIF
         IF (LEN_TRIM(ICAO_IN).EQ.5) THEN
            IF (ICAO_IN(1:5).EQ.ICAO(1:5)) EXIT
         ENDIF
         IF (LEN_TRIM(ICAO_IN).EQ.6) THEN
            IF (ICAO_IN(1:6).EQ.ICAO(1:6)) EXIT
         ENDIF
      ENDDO
      APORT=PORTINFO(1:30)
      CITY=PORTINFO(37:66)
      LATNS=PORTINFO(74:74)
      LATD=CHAR2INT(PORTINFO(67:68),2) 
      LATM=CHAR2INT(PORTINFO(69:70),2) 
      LATS=CHAR2INT(PORTINFO(71:72),2) 
      LONEW=PORTINFO(83:83)
      LOND=CHAR2INT(PORTINFO(75:77),3)
      LONM=CHAR2INT(PORTINFO(78:79),2) 
      LONS=CHAR2INT(PORTINFO(80:81),2) 
      ALTF=CHAR2INT(PORTINFO(84:89),6) 
      CLOSE(98)
      CALL CLS
      CALL MENUHEADER
      WRITE(*,73301) ICAO_IN, ALTCODE(ICAO_IN) 
      WRITE(*,73302) APORT
      WRITE(*,73303) CITY
      WRITE(*,73304) LATD, LATM, LATS, NSS(LATNS)
      WRITE(*,73305) LOND, LONM, LONS, EWS(LONEW)
      WRITE(*,73306) REAL(ALTF,KIND=8)/10.0
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
        PRINT*,' '
73301 FORMAT(10X,'ICAO CODE: ',A6, '                OTHER CODE: ',A3)
73302 FORMAT(10X,'AIRPORT NAME: ',A30)
73303 FORMAT(10X,'CITY NAME: ',A30)
73304 FORMAT(10X,'LATITUDE: ',I2,' DEGS ',I2,' MINS ',I2,' SECS ',A1)
73305 FORMAT(10X,'LONGITUDE: ',I3,' DEGS ',I2,' MINS ',I2,' SECS ',A1)
73306 FORMAT(10X,'ALTITUDE: ',F7.1, ' FEET')
      WRITE(*,*)' '
      CALL OOPS (' ',1)
      CALL AIRPORTS
      RETURN   
73307 CLOSE(98)
      IF (DIAGNOSE.EQ.'YES')WRITE(40,*)'REP_BY_ICAO FAILED FOR ',ICAO_IN
      WRITE(*,*)'REP_BY_ICAO FAILED FOR ',ICAO_IN
      FRACK=.TRUE. 
      CALL OOPS ('ICAO CODE NOT FOUND IN DATABASES',32)
      CALL CLS
      CALL MENUHEADER
73308 FORMAT(A91)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION NSS(A)
! INTERPRETS NORTH/SOUTH CODING IN AIRPORT DATABASES          
      CHARACTER(1)::A,B,NSS

      IF (A.EQ.'1'.OR.A.EQ.'2'.OR.A.EQ.'3'.OR.A.EQ.'4'.OR.A.EQ.'5'.OR.  &
     & A.EQ.'6'.OR.A.EQ.'7'.OR.A.EQ.'8'.OR.A.EQ.'9'.OR.A.EQ.'0') THEN
        B='N'
      ELSE
        B='S'
      ENDIF
      NSS=B 

      END FUNCTION
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION NS(A)
! CONVERTS NORTH/SOUTH CODING TO +/- 
        CHARACTER(1)::A
        REAL(8)::NS
      IF (A.EQ.'N') THEN
        NS=1.
      ELSE
        NS=-1.
      ENDIF
!      Print*, A,' latitude is ',NS  !diagnostic

      END FUNCTION
!                                                                      7
!----6-----------------------------------------------------------------2

      FUNCTION EW(A)
! CONVERTS EAST/WEST CODING TO +/- 
        CHARACTER(1)::A
        REAL(8)::EW
        
      IF (A.EQ.'E') THEN
        EW=1.
      ELSE
        EW=-1.
      ENDIF
!      Print*, A,' latitude is ',EW !diagnostic 

      END FUNCTION
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION EWS(A)
! INTERPRETS EAST/WEST CODING IN AIRPORT DATABASES
        CHARACTER(1)::A,B,EWS
        
      IF (A.EQ.'1'.OR.A.EQ.'2'.OR.A.EQ.'3'.OR.A.EQ.'4'.OR.A.EQ.'5'.OR.  &
     & A.EQ.'6'.OR.A.EQ.'7'.OR.A.EQ.'8'.OR.A.EQ.'9'.OR.A.EQ.'0') THEN
        B='W'
      ELSE
        B='E'
      ENDIF
      EWS=B 

      END FUNCTION
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION ALTCODE(A)
! FINDS ALTERNATES TO ICAO CODES FROM FILE 'CODES'   
         CHARACTER(1) :: JUNK
         CHARACTER(6) :: A,AA
         CHARACTER(4) :: AAA, B
         CHARACTER(3) :: C, ALTCODE
         INTEGER(4) :: I
         CHARACTER(3)::DIAGNOSE
         DIAGNOSE='NO!'
         AA=ADJUSTL(A); AAA=AA(1:4)
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'TESTING ALTCODE'  
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'CODE IN:', AAA   
         DO 
            READ(30,73330,ERR=73331,END=73331) C, JUNK, B
            IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'VS:', C, B
            IF(AAA.EQ.B) THEN
               REWIND(30) 
               ALTCODE=C
               RETURN
            ENDIF
         ENDDO 
73331 CONTINUE
      REWIND(30) 
      ALTCODE='---'
73330 FORMAT(A3,A1,A4) 
      END FUNCTION  
!                                                                      7
!----6-----------------------------------------------------------------2
! 73500 
      SUBROUTINE ADDAPORT
      IMPLICIT NONE
! Add an airport to NEWPORTS.DAT   
      INTEGER(4)::CHAR2INT,I,OLDNUM
      CHARACTER(30)::CITY,PORT
      CHARACTER(21)::NEWPORTPATH
      CHARACTER(14)::CODEPATH
      CHARACTER(6)::ICAO,ALTF
      CHARACTER(5)::OS
      CHARACTER(3)::IATA,LOD
      CHARACTER(2)::LOM,LOS,LAD,LAM,LAS
      CHARACTER(1)::NSS,EWS,YN
      CHARACTER(89)::NEWPORT,LASTPORT,TESTPORT
      CHARACTER(89), DIMENSION(200)::OLDPORTS
      CHARACTER(4)::OUTPUT
      CHARACTER(12)::VIEWER
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT  
        
         DIAGNOSE='NO!'
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
         CODEPATH='AIRPORTS\CODES'
         NEWPORTPATH='AIRPORTS\NEWPORTS.DAT'
      ELSE
         CODEPATH='AIRPORTS/CODES'
         NEWPORTPATH='AIRPORTS/NEWPORTS.DAT'
      ENDIF
      DO I = 1,89
         NEWPORT(I:I)='0'
      ENDDO

         
      CALL CLS
      CALL MENUHEADER
        PRINT*,' ' 
      WRITE(*,*)'          AIRPORT ENTRY FORM                      '
      PRINT*, ' '
      WRITE(*,*)'          <1>  Airport name (UP TO 30 CHARACTERS).'
      CALL READKB(PORT,30)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,73501) PORT 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,73501) PORT 
      WRITE(*,*)'          <2>  City name (UP TO 30 CHARACTERS).'
      CALL READKB(CITY,30)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,73501) CITY 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,73501) CITY 
      WRITE(*,*)'          <3>  ICAO code (UP TO 6 CHARACTERS).'
      CALL READKB(ICAO,6)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,73501) ICAO 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,73501) ICAO

73502 WRITE(*,*)'          <4>  Add a 3 letter code? (YN)'
      CALL READKB(YN,1)
      IF (YN.EQ.'Y') THEN
      WRITE(*,*)' -----------------INSTRUCTIONS--------------------' 
      WRITE(*,*)' The format of the 3-Character code match file is:'
      WRITE(*,*)' XXX,ICAO_CODE. (i.e., a 3Char_code followed by a '
      WRITE(*,*)' comma followed by the ICAO code'
      WRITE(*,*)' Add data at any place in the file. It need not be'
      WRITE(*,*)' at the end. But the search uses the first match'
      WRITE(*,*)' found with any particulat ICAO code.' 
         CLOSE (30)
          
         CALL SHOWPICK(CODEPATH,14)
         OPEN (UNIT=30,FILE=CODEPATH,STATUS='OLD')      
      ELSEIF (YN.EQ.'N') THEN
           CONTINUE
      ELSE
           PRINT*,'          Entry invalid, must be Y or N.'
         PRINT*,''
         GOTO 73502
      ENDIF
      WRITE(*,*)'          LATITUDE'
      WRITE(*,*)'          <6a>  NORTH OR SOUTH <N/S>'
      CALL READKB(NSS,1)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) NSS 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) NSS
      IF (NSS.EQ.'S') THEN 
           NSS='}'
      ELSE
         NSS='0'
      ENDIF        
      WRITE(*,*)'          <6b>  DEGREES (UP TO 2 CHARACTERS).'
      CALL READKB(LAD,2)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) LAD 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) LAD 
      WRITE(*,*)'          <6c>  MINUTES (UP TO 2 CHARACTERS).'
      CALL READKB(LAM,2)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) LAM 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) LAM 
      WRITE(*,*)'          <6d>  SECONDS (UP TO 2 CHARACTERS).'
      CALL READKB(LAS,2)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) LAS 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) LAS 
      WRITE(*,*)'          LONGITUDE'
      WRITE(*,*)'          <7a>  EAST OR WEST <E/W>.'
      CALL READKB(EWS,1)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) EWS 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) EWS 
      IF (EWS.EQ.'E') THEN 
           EWS='}'
      ELSE
         EWS='0'
      ENDIF        
      WRITE(*,*)'          <7b>  DEGREES (UP TO 3 CHARACTERS).'
      CALL READKB(LOD,3)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) LOD 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) LOD 
      WRITE(*,*)'          <7c>  MINUTES (UP TO 2 CHARACTERS).'
      CALL READKB(LOM,2)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) LOM
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) LOM 
      WRITE(*,*)'          <7d>  SECONDS (UP TO 2 CHARACTERS).'
      CALL READKB(LOS,2)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) LOS 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) LOS 
      PRINT*, ' '
      WRITE(*,*)'          <8>  ALTITUDE IN FEET (UP TO 5 CHARACTERS).'
      CALL READKB(ALTF,5)
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) ALTF 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) ALTF 
      PRINT*, ' '

      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'Defining NEWPORT' 
      IF (DIAGNOSE.EQ.'YES') WRITE(*,*) 'Defining NEWPORT' 
        
      NEWPORT(1:30)=PORT
      NEWPORT(31:36)=ICAO
      NEWPORT(37:66)=CITY
      IF(CHAR2INT(LAD,2).GT.9) THEN
         NEWPORT(67:68)=LAD
      ELSE
           NEWPORT(67:67)='0'
           NEWPORT(68:68)=LAD
      ENDIF  
      IF(CHAR2INT(LAM,2).GT.9) THEN
         NEWPORT(69:70)=LAM
      ELSE
           NEWPORT(69:69)='0'
           NEWPORT(70:70)=LAM
      ENDIF  
      IF(CHAR2INT(LAS,2).GT.9) THEN
         NEWPORT(71:72)=LAS
      ELSE
           NEWPORT(71:71)='0'
           NEWPORT(72:72)=LAS
      ENDIF
      NEWPORT(73:73)='0'
      NEWPORT(74:74)=NSS
      IF(CHAR2INT(LOD,3).GT.9 .AND. CHAR2INT(LOD,3).LT.100) THEN
           NEWPORT(75:75)='0'
         NEWPORT(76:77)=LOD
      ELSEIF(CHAR2INT(LOD,3).GT.99) THEN
           NEWPORT(75:77)=LOD
      ELSE
           NEWPORT(75:76)='00'
         NEWPORT(77:77)=LOD
      ENDIF  
      IF(CHAR2INT(LOM,2).GT.9) THEN
         NEWPORT(78:79)=LOM
      ELSE
         NEWPORT(78:78)='0'
         NEWPORT(79:79)=LOM
      ENDIF  
      IF(CHAR2INT(LOS,2).GT.9) THEN
         NEWPORT(80:81)=LOS
      ELSE
         NEWPORT(80:80)='0'
         NEWPORT(81:81)=LOS
      ENDIF
      NEWPORT(82:82)='0'
      NEWPORT(83:83)=EWS
      IF(CHAR2INT(ALTF,5).GT.9 .AND. CHAR2INT(ALTF,5).LT.100) THEN
         NEWPORT(84:86)='000'
         NEWPORT(87:88)=ALTF
      ELSEIF(CHAR2INT(ALTF,5).GT.99 .AND. CHAR2INT(ALTF,5).LT.1000) THEN
         NEWPORT(84:85)='00'
         NEWPORT(86:88)=ALTF
      ELSEIF(CHAR2INT(ALTF,5).GT.999 .AND. CHAR2INT(ALTF,5).LT.10000)   &
     & THEN
         NEWPORT(84:84)='0'
         NEWPORT(85:88)=ALTF
      ELSEIF(CHAR2INT(ALTF,5).LT.10) THEN
         NEWPORT(84:87)='0000'
         NEWPORT(88:88)=ALTF
      ELSE
         NEWPORT(84:88)=ALTF
      ENDIF  
        
      NEWPORT(89:89)='0'
!      NEWPORT(90:90)=CHAR(13)!CR
!      NEWPORT(91:91)=CHAR(10)!LF
        
      IF (DIAGNOSE.EQ.'YES') THEN
      WRITE(40,*) 'NEWPORT DEFINED'
      WRITE(*,*) 'CITY NAME: ', NEWPORT(37:66)
      WRITE(*,*) 'AIRPORT NAME: ', NEWPORT(1:30)
      WRITE(*,*) 'ICAO CODE: ', NEWPORT(31:36)!,' ALTERNATE CODE: ',IATA
      WRITE(40,*) 'CITY NAME: ', NEWPORT(37:66)
      WRITE(40,*) 'AIRPORT NAME: ', NEWPORT(1:30)
      WRITE(40,*)'ICAO CODE: ', NEWPORT(31:36)!,' ALTERNATE CODE: ',IATA
      CALL OOPS('_',1)
      ENDIF
        
      CALL CLS
        CALL MENUHEADER
      WRITE(*,*) 'CITY NAME: ', NEWPORT(37:66)
      WRITE(*,*) 'AIRPORT NAME: ', NEWPORT(1:30)
      WRITE(*,*) 'ICAO CODE: ', NEWPORT(31:36)!,' ALTERNATE CODE: ',IATA
      WRITE(*,*) 'COORDINATES:' 
      IF (NSS.NE.'0') THEN
         WRITE(*,*) CHAR2INT(NEWPORT(67:68),2),'DEGS',                  &
     &                CHAR2INT(NEWPORT(69:70),2),'MINS',                &
     &                CHAR2INT(NEWPORT(71:72),2),'SECS SOUTH LAT'
      ELSE
         WRITE(*,*) CHAR2INT(NEWPORT(67:68),2),'DEGS',                  &
     &                CHAR2INT(NEWPORT(69:70),2),'MINS',                &
     &                CHAR2INT(NEWPORT(71:72),2),'SECS NORTH LAT'
      ENDIF
      IF (EWS.NE.'0') THEN
         WRITE(*,*) CHAR2INT(NEWPORT(75:77),3),'DEGS',                  &
     &                CHAR2INT(NEWPORT(78:79),2),'MINS',                &
     &                CHAR2INT(NEWPORT(80:81),2),'SECS EAST LON'
      ELSE
         WRITE(*,*) CHAR2INT(NEWPORT(75:77),3),'DEGS',                  &
     &                CHAR2INT(NEWPORT(78:79),2),'MINS',                &
     &                CHAR2INT(NEWPORT(80:81),2),'SECS WEST LON'
      ENDIF
      WRITE(*,*) 'ALTITUDE (FT): ',CHAR2INT(NEWPORT(84:88),5)
      PRINT*,' '
      PRINT*,'          IS THIS INFORMATION CORRECT <Y/N>?' 
      LASTPORT(1:45) ='END OF FILE                   ----  END OF FI'   
      LASTPORT(46:89)='LE                   99999990999999999999990'
!      LASTPORT(90:90)=CHAR(13)!CR
!      LASTPORT(91:91)=CHAR(10)!LF
      CALL READKB(YN,1)
        IF (YN.EQ.'Y') THEN
           REWIND(32)
           DO I=1,200
              READ(32,73504) OLDPORTS(I)
              TESTPORT=OLDPORTS(I)
              IF (TESTPORT(1:11).EQ.'END OF FILE')THEN
                 OLDNUM=I-1
                 EXIT
              ENDIF
           ENDDO
           CLOSE(32)
           OPEN(UNIT=32,FILE=NEWPORTPATH,STATUS='UNKNOWN')  
           WRITE(32,73504) NEWPORT    
           DO I = 1,OLDNUM
              WRITE(32,73504) OLDPORTS(I)    
           ENDDO
           WRITE(32,73504) LASTPORT    
           CLOSE(32)
           OPEN(UNIT=32,FILE=NEWPORTPATH,STATUS='OLD')
           CALL MAKE_NDX(91,7000,OS) 
           CALL MENU_MAIN
        ELSE
           CALL MENU_AIRPORT
        ENDIF      
73501 FORMAT('ENTERED',A30)
73503 FORMAT(A91)
73504 FORMAT(A89)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! 73600 
      SUBROUTINE ADDACODE
      
      CHARACTER(14)::CODEPATH
!      CHARACTER(5)::OS

      CALL CLS
      CALL MENUHEADER
      WRITE(*,*)' -----------------INSTRUCTIONS--------------------' 
      WRITE(*,*)''
      WRITE(*,*)' The format of the 3-Character code match file is:'
      WRITE(*,*)' XXX,ICAO_CODE. (i.e., a 3Char_code followed by a '
      WRITE(*,*)' comma followed by the ICAO code.'
      WRITE(*,*)''
      WRITE(*,*)' Add data at any place in the file. It need not be'
      WRITE(*,*)' at the end. But the search uses the first match'
      WRITE(*,*)' found with any particular ICAO code.' 
      WRITE(*,*)''
      WRITE(*,*)''
      WRITE(*,*)''
      WRITE(*,*)''
      WRITE(*,*)''

      CLOSE(UNIT=30)
            
         CODEPATH='AIRPORTS/CODES'

      CALL SHOWPICK(CODEPATH,LEN(CODEPATH))
      CALL OOPS('         When finished... ',26) 
      OPEN(UNIT=30,FILE=CODEPATH,STATUS='OLD')

      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! 74000 
      SUBROUTINE MENU_LOCATIONS
!     menu for single locations
      CHARACTER(1)::CHOICE

74000 CALL CLS      
      CALL MENUHEADER
      WRITE(*,74004)'               LOCATION MENU                      '
      PRINT*, ' '
      WRITE(*,74004)'<1>  Calculate the dose rate at a single location.'
      PRINT*, ' '
      WRITE(*,74004)'<2>  Calculate dose rates for locations in a      '
      WRITE(*,74004)'     locations file (*.LOC).                      '
      PRINT*, ' '
      WRITE(*,74004)'<3>  Open a location file (*.LOC).                '
      WRITE(*,74004)'     (uses the default text editor)               '
      PRINT*, ' '
      WRITE(*,74004)'<4>  Open a dose rate archive file (*.ANS).       '
      PRINT*, ' '
      WRITE(*,74004)'<5>  Open the HELP file.                          '
      PRINT*, ' '
      WRITE(*,74004)'<6>  Return to Main Menu.                         '
      PRINT*, ' '
      WRITE(*,74004)'<7>  Exit program.                                '
      PRINT*,''
      WRITE(*,74001) '.' 
      CALL READKB(CHOICE,1)
74001 FORMAT(10X,'Type 1, 2, 3, 4, 5, 6, or 7 and press <ENTER> ',A1)
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'4'.OR.CHOICE.EQ.'5'.OR.CHOICE.EQ.'6'.OR.           &
     &    CHOICE.EQ.'7') THEN
         GOTO 74002
      ELSE 
         GOTO 74000 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
74002 CONTINUE
! DIAGNOSTIC PRINT   PRINT*, "SUCCESFUL DATA ENTRY ", CHOICE
      IF (CHOICE.EQ.'1') CALL ONESPOT                    !20000
      IF (CHOICE.EQ.'2') CALL RUN_LOCATIONS              !21000 
      IF (CHOICE.EQ.'3') CALL NEWLOCS                    !22000
      IF (CHOICE.EQ.'4') CALL SHOWDOSERATES              !
      IF (CHOICE.EQ.'5') CALL SHOWHELP                   !
      IF (CHOICE.EQ.'7') STOP 
      CALL MENU_MAIN
74003 FORMAT (A1)
74004 FORMAT (10X,A50)
74005 FORMAT (10X,A9)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! 75000 
      SUBROUTINE SOLARMOD
      CALL MENU_SOLARMOD   
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! 75000 
      SUBROUTINE MENU_SOLARMOD

      CHARACTER(1)::CHOICE

75000 CALL CLS      
      CALL MENUHEADER
      WRITE(*,75004)'           SOLAR DATA MENU                        '
      PRINT*, ' '
      WRITE(*,75004)'<1>  View/modify heliocentric potential data.     '
      PRINT*, ' '
      WRITE(*,75004)'<2>  View/modify Forbush effect data.             '
      PRINT*, ' '
      WRITE(*,75004)'<3>  View/modify Kp index data.                   '
      PRINT*, ' '
      WRITE(*,75004)'<4>  Return to Main Menu.                         '
      PRINT*, ' '
      WRITE(*,75004)'<5>  Exit program.                                '
      PRINT*,''
      PRINT*,''
      PRINT*,''
      PRINT*,''
      PRINT*,''
      WRITE(*,75001) '.' 
      CALL READKB(CHOICE,1)
75001 FORMAT(10X,'Type 1, 2, 3, 4, or 5 and press <ENTER> ',A1)
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'4'.OR.CHOICE.EQ.'5'.OR.                            &
     &    CHOICE.EQ.'q'.OR.CHOICE.EQ.'Q') THEN
         GOTO 75002
      ELSE 
         GOTO 75000 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
75002 CONTINUE
      IF (CHOICE.EQ.'Q'.OR.CHOICE.EQ.'q') CHOICE='6' 
! DIAGNOSTIC PRINT   PRINT*, "SUCCESFUL DATA ENTRY ", CHOICE
      IF (CHOICE.EQ.'1') CALL POTENTIALS                 
      IF (CHOICE.EQ.'3') CALL SHOWKP                    !75600     
      IF (CHOICE.EQ.'2') CALL SHOWFORBUSH               !75700
      IF (CHOICE.EQ.'4') CALL MENU_MAIN                  
      IF (CHOICE.EQ.'5') STOP 
      CALL MENU_MAIN
75003 FORMAT (A1)
75004 FORMAT (10X,A50)
75005 FORMAT (10X,A9)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE POTENTIALS
      CALL MENU_POTENTIALS   
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! 75000 
      SUBROUTINE MENU_POTENTIALS

      CHARACTER(1)::CHOICE

75000 CALL CLS      
      CALL MENUHEADER
      WRITE(*,75004)'     HELIOCENTRIC POTENTIAL MENU                  '
      PRINT*, ' '
      WRITE(*,75004)'<1>  Look up the potential for a date.            '
      PRINT*, ' '
      WRITE(*,75004)'<2>  Add a potential to the monthly database.     '
      PRINT*, ' '
      WRITE(*,75004)'<3>  Add a potential for a specific day or a      '
      WRITE(*,75004)'     temporary value for a past potential.        '
      PRINT*, ' '
      WRITE(*,75004)'<4>  Return to Main Menu.                         '
      PRINT*, ' '
      WRITE(*,75004)'<5>  Exit program.                                '
      PRINT*,''
      PRINT*, ' '
      PRINT*,''
      PRINT*,''
      WRITE(*,75001) '.' 
      CALL READKB(CHOICE,1)
75001 FORMAT(10X,'Type 1, 2, 3, 4, 5, or 6 and press <ENTER> ',A1)
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'4'.OR.CHOICE.EQ.'5'.OR.                            &
     &    CHOICE.EQ.'q'.OR.CHOICE.EQ.'Q') THEN
         GOTO 75002
      ELSE 
         GOTO 75000 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
75002 CONTINUE
      IF (CHOICE.EQ.'Q'.OR.CHOICE.EQ.'q') CHOICE='6' 
! DIAGNOSTIC PRINT   PRINT*, "SUCCESFUL DATA ENTRY ", CHOICE
      IF (CHOICE.EQ.'1') CALL VIEWHP                     !75100
      IF (CHOICE.EQ.'2') CALL ADD2L99                    !75200 
      IF (CHOICE.EQ.'3') CALL ADDADAILY                  !75300 
      IF (CHOICE.EQ.'5') CALL MENU_MAIN                  !75500
      IF (CHOICE.EQ.'6') STOP 
      CALL MENU_MAIN
75003 FORMAT (A1)
75004 FORMAT (10X,A50)
75005 FORMAT (10X,A9)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! 75100
      SUBROUTINE VIEWHP

      CHARACTER(10)::DATESTR
      CHARACTER(1)::CHOICE
      INTEGER(4)::Y,M,D,HP

      CALL CLS
75100 CALL MENUHEADER 
      WRITE(*,75104)'     VIEW HELIOCENTRIC POTENTIALS                 '
      PRINT*, ' '
      WRITE(*,75104)'     Enter a date as YYYY/MM/DD                   '
      WRITE(*,75104)'     (Enter YYYY/00 for a yearly average)         ' 
      WRITE(*,75104)'     (Enter YYYY/MM/00 for a monthly average)     ' 
      PRINT*, ' '
      PRINT*, '                   -- OR --                             '
      PRINT*,''
      WRITE(*,75104)'<1>  Return to Heliocentric Potential MENU.       '
      PRINT*, ' '
      WRITE(*,75104)'<2>  Return to Main Menu.                         '
      PRINT*, ' '
      WRITE(*,75104)'<3>  Exit program.                                '
      PRINT*,''
      PRINT*,''
      WRITE(*,75101) '.' 
      READ(*,75106) DATESTR
      IF (LEN_TRIM(DATESTR).GT.1) THEN 
         CALL DATE2YMD(DATESTR,Y,M,D)
         CALL DATE2HP(Y,M,D,HP)
         WRITE (*,*) '           For YEAR:',Y,' MONTH:',M ,' and DAY:',D
         WRITE (*,75105) HP
         CALL OOPS(' ',1)
         CHOICE='1'
      ELSE
         CHOICE=DATESTR(1:1)
      ENDIF
75101 FORMAT(10X,'Type date, 1, 2, or 3 and press <ENTER> ',A1)
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'q'.OR.CHOICE.EQ.'Q') THEN
         GOTO 75102
      ELSE 
         GOTO 75100 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
75102 CONTINUE
      IF (CHOICE.EQ.'Q'.OR.CHOICE.EQ.'q') CHOICE='3' 
      IF (CHOICE.EQ.'1') CALL MENU_POTENTIALS
      IF (CHOICE.EQ.'2') CALL MENU_MAIN
      IF (CHOICE.EQ.'3') STOP 
75104 FORMAT (10X,A50)
75105 FORMAT (10X,'Heliocentric potential was: ',I5,' MV')
75106 FORMAT (A10)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! 75200 ADD to monthly database of heliocentric potentials
      SUBROUTINE ADD2L99
       PRINT*,'ADD NEW DATE AND POTENTIAL TO END OF FILE MV-DATES'
       PRINT*,'E.G. 12/2020,   1345'
         CLOSE(28)
! EDIT MV-DATES WITH DEFAULT TEXT EDITOR
         CALL SHOWPICK('SOLARMOD/MV-DATES.L99',21)
!         OPEN(UNIT=28,FILE='SOLARMOD/MV-DATES.L99',STATUS='OLD')
      END SUBROUTINE    
!                                                                      7
!----6-----------------------------------------------------------------2
! 75300 ADD to user database of heliocentric potentials
      SUBROUTINE ADDADAILY
       PRINT*,'ADD NEW DATE AND POTENTIAL TO END OF FILE MORDATES.2K'
       PRINT*,'E.G. 2020/12/25,   1345'
         CLOSE(31)
! EDIT MORDATES WITH DEFAULT TEXT EDITOR
         CALL SHOWPICK('SOLARMOD/MORDATES.2K',20)
!         OPEN(UNIT=31,FILE='SOLARMOD/MORDATES.2K',STATUS='OLD')
      END SUBROUTINE    
!                                                                      7
!----6-----------------------------------------------------------------2
! 75400 ADD to user database of heliocentric potentials
      SUBROUTINE ADDAMONTH
         PRINT*,'ADD NEW DATE AND POTENTIAL TO END OF FILE MORDATES.2K'
         PRINT*,'E.G. 12/2020,   1345'
         CLOSE(31)
! EDIT MORDATES WITH DEFAULT TEXT EDITOR
         CALL SHOWPICK('SOLARMOD/MORDATES.2K',20)
!         OPEN(UNIT=31,FILE='SOLARMOD/MORDATES.2K',STATUS='OLD')
      END SUBROUTINE    
!                                                                      7
!----6-----------------------------------------------------------------2
! 76000
      SUBROUTINE OUTPUTS
      CHARACTER(1)::CHOICE
      CALL MENU_OUTPUT(CHOICE)
!     DEFAULT IS SET OUTPUT TO FILE ONLY        
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE MENU_OUTPUT(CHOICE)

      CHARACTER(1)::CHOICE, OUTPUT, JUNK
      CHARACTER(12)::VARIABLE, SETTING
      CHARACTER(9)::CS(3)
        CHARACTER(50)::MYTEXT
      INTEGER(4)::I
      CHARACTER(3)::DIAGNOSE
           DIAGNOSE='NO!'

76000 CALL CLS
      CALL MENUHEADER
      WRITE(*,76005)'                OUTPUT MENU                       '
      PRINT*, ' '
      WRITE(*,76005)'<1>  Review/Edit the default settings in CARI.INI.'
      PRINT*, ' '
      WRITE(*,76005)'<2>  Review/Edit defaults specified for shell     '
      WRITE(*,76005)'     script use (i.e., batch mode) in DEFAULT.DAT.'
      PRINT*, ' '
      WRITE(*,76005)'<3>  Review a city-pair flight dose (OUT) archive.'
      PRINT*, ' '
      WRITE(*,76005)'<4>  Review a waypoint flight dose (SUM) archive. '
      PRINT*, ' '
      WRITE(*,76005)'<5>  Review a location dose rate (ANS) archive.   '
      PRINT*, ' '
      WRITE(*,76005)'<6>  Return to Main Menu.                         '
      PRINT*, ' '
      WRITE(*,76005)'<7>  Exit program.                                '
      PRINT*,''
      WRITE(*,76001) '.' 
      CALL READKB(CHOICE,1)
76001 FORMAT('     Type 1, 2, 3, 4, 5, 6, or 7 and press <ENTER> ',A1)
      IF (CHOICE.EQ.'1'.OR.CHOICE.EQ.'2'.OR.CHOICE.EQ.'3'.OR.           &
     &    CHOICE.EQ.'4'.OR.CHOICE.EQ.'5'.OR.CHOICE.EQ.'6'.OR.           &
     &    CHOICE.EQ.'7') THEN
         GOTO 76002
      ELSE 
         GOTO 76000 ! TRY AGAIN, NOT A VALID CHOICE
      END IF
76002 CONTINUE
! DIAGNOSTIC PRINT 
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*) 'SUCCESFUL DATA ENTRY ', CHOICE
      IF (CHOICE.EQ.'1') THEN 
         CALL SHOWPICK('CARI.INI',8)
         CALL GETINI !must reload after changes
      ELSEIF(CHOICE.EQ.'2') THEN 
         CALL SHOWPICK('DEFAULT.INP',11)
      ELSEIF(CHOICE.EQ.'3') THEN 
         CALL SHOWDOSES
      ELSEIF(CHOICE.EQ.'4') THEN 
         CALL SHOWDOSES2
      ELSEIF(CHOICE.EQ.'5') THEN 
         CALL SHOWDOSERATES
      ELSEIF(CHOICE.EQ.'7') THEN 
         STOP 
      ELSE
         CALL MENU_MAIN
      ENDIF
76003 FORMAT (A1)
76005 FORMAT (10X,A50)

      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE SHOWDOSES

!  OPEN A FLIGHT DOSE ARCHIVE

      CALL PICKFROMLIST('OUT') 

      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE SHOWDOSES2

!  OPEN A FLIGHT DOSE ARCHIVE

      CALL PICKFROMLIST('SUM') 

      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE SHOWDOSERATES

!  OPEN A DOSE RATE ARCHIVE (RESULTS AT SINGLE LOCATIONS)

      CALL PICKFROMLIST('ANS') 

      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE PICKFROMLIST(FILE_EXT)

      CHARACTER(2)::SELECTION,GOON
      CHARACTER(3)::FILE_EXT
      CHARACTER(60)::MESSAGE, FILE_LIST(5,16)
      CHARACTER(52)::FILENAME 

      INTEGER(4)::FILENUM,SELECTED,PAGE,I,J,M,N,FNUM,CHAR2INT
! .INI read block
      CHARACTER(10)::INIVAR
      CHARACTER(12)::INIVAL
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
      INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR
!      
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN ! Win/DOS syntax
         CALL SYSTEM('DIR /B *.'//FILE_EXT//' > FILELIST.TXT')
      ELSE 
         CALL SYSTEM('ls -1 *.'//FILE_EXT//' > FILELIST.TXT')
      ENDIF

      OPEN(UNIT=60,FILE='FILELIST.TXT',STATUS='OLD')
     
      N=-1
      DO I=1,5
         DO J = 1,16
            N=N+1
            READ(60,*,ERR=76014,END=76014) FILE_LIST(I,J)
         ENDDO
      ENDDO
76014 CONTINUE
      CLOSE(60)  
      SELECT CASE (N)
      CASE (0)        ! No such files to view
         MESSAGE = ' No such files exist in this directory!'
         CALL OOPS(MESSAGE,39)
         RETURN 
      CASE (1:16)     ! Show what there is  
         MESSAGE = ' '
      CASE (17:80)    ! Start with the first page 
         MESSAGE = ' '
      CASE DEFAULT    ! Start with the first page, but there are more 
         MESSAGE = ' '  !    than 80 files to choose from 
      END SELECT 

      M=16 !number of files to show on a page
      PAGE=1
76009 CALL CLS
         WRITE(*,*) ' '
         PRINT*,' PAGE ',PAGE,' OF ',(N+16)/16  
         WRITE(*,*) ' '
         WRITE(*,*)'    Showing available files 16 at a time.'
         WRITE(*,*) ' '
! SHOW FILES 1 PAGE AT  A TIME
         IF (M>N-(PAGE-1)*16) THEN ! DO NOT HAVE A FULL PAGE 
           M=N
         ELSE
           M=16
         ENDIF
         DO I = 1, M ! 
           WRITE(*,76010) I+(PAGE-1)*16, FILE_LIST(PAGE,I)
         ENDDO
         PRINT *,' '
         SELECT CASE(PAGE)
             CASE(1)
                IF (N.LT.17) THEN 
                    WRITE(*,76011)
                ELSE 
                    WRITE(*,76012) 
                ENDIF            
              CASE(2:4)
                    WRITE(*,76013) 
              CASE(5)
                    WRITE(*,76015) 
           END SELECT    
         CALL READKB(SELECTION,2)
      IF (SELECTION=='n' .OR. SELECTION=='N') THEN
         PAGE = PAGE + 1
         GOTO 76009
      ENDIF
      IF (SELECTION=='P' .OR. SELECTION=='p') THEN
         PAGE = PAGE - 1
         GOTO 76009
      ENDIF
      IF (SELECTION=='q' .OR. SELECTION=='Q') CALL MENU_MAIN
      
      FILENUM=CHAR2INT(SELECTION,2)
      FILENAME = FILE_LIST(PAGE,FILENUM)
      CALL SHOWPICK(FILENAME,LEN(FILENAME))
76010 FORMAT (10X,'<',I2,'> ',A60)
76011 FORMAT (10X,'Enter a file number or quit <Q,q>: ' )
76012 FORMAT (10X,'Enter a file number, see next page <N,n>'            &
     &  ,', or quit <Q,q>: ') 
76013 FORMAT (10X,'Enter a file number, see next page <N,n>'            &
     &  ,', see preivous page <P,p>, or quit <Q,q>: ')  
76015 FORMAT (10X,'Enter a file number',                                &
     &   ', see preivous page <P,p>, or quit <Q,q>: ')  
      END SUBROUTINE PICKFROMLIST
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE PICKBIG(FILENAME)

      CHARACTER(2)::SELECTION,GOON
      CHARACTER(3)::FILE_EXT
      CHARACTER(60)::MESSAGE, FILE_LIST(5,16)
      CHARACTER(52)::FILENAME 

      INTEGER(4)::FILENUM,SELECTED,PAGE,I,J,M,N,FNUM,CHAR2INT
! .INI read block
      CHARACTER(10)::INIVAR
      CHARACTER(12)::INIVAL
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
      INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR
!      
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN ! Win/DOS syntax
         CALL SYSTEM('DIR /B *.BIG > FILELIST.TXT')
      ELSE 
         CALL SYSTEM('ls -1 *.BIG > FILELIST.TXT')
      ENDIF
      
      OPEN(UNIT=60,FILE='FILELIST.TXT',STATUS='OLD')
     
      N=-1
      DO I=1,5
         DO J = 1,16
            N=N+1
            READ(60,*,ERR=76014,END=76014) FILE_LIST(I,J)
         ENDDO
      ENDDO
76014 CONTINUE
      CLOSE(60)  
      SELECT CASE (N)
      CASE (0)        ! No such files to view
         MESSAGE = ' No such files exist in this directory!'
         CALL OOPS(MESSAGE,39)
         RETURN 
      CASE (1:16)     ! Show what there is  
         MESSAGE = ' '
      CASE (17:80)    ! Start with the first page 
         MESSAGE = ' '
      CASE DEFAULT    ! Start with the first page, but there are more 
         MESSAGE = ' '  !    than 80 files to choose from 
      END SELECT 

      M=16 !number of files to show on a page
      PAGE=1
76009 CALL CLS
         WRITE(*,*) ' '
         PRINT*,' PAGE ',PAGE,' OF ',(N+16)/16  
         WRITE(*,*) ' '
         WRITE(*,*)'    Showing available files 16 at a time.'
         WRITE(*,*) ' '
! SHOW FILES 1 PAGE AT  A TIME
         IF (M>N-(PAGE-1)*16) THEN ! DO NOT HAVE A FULL PAGE 
           M=N
         ELSE
           M=16
         ENDIF
         DO I = 1, M ! 
           WRITE(*,76010) I+(PAGE-1)*16, FILE_LIST(PAGE,I)
         ENDDO
         PRINT *,' '
         SELECT CASE(PAGE)
             CASE(1)
                IF (N.LT.17) THEN 
                    WRITE(*,76011)
                 ELSE 
                  WRITE(*,76012) 
               ENDIF            
              CASE(2:4)
                    WRITE(*,76013) 
              CASE(5)
                    WRITE(*,76015) 
           END SELECT    
         CALL READKB(SELECTION,2)
      IF (SELECTION=='n' .OR. SELECTION=='N') THEN
         PAGE = PAGE + 1
         GOTO 76009
      ENDIF
      IF (SELECTION=='P' .OR. SELECTION=='p') THEN
         PAGE = PAGE - 1
         GOTO 76009
      ENDIF
      IF (SELECTION=='q' .OR. SELECTION=='Q') CALL MENU_MAIN
      
      FILENUM=CHAR2INT(SELECTION,2)
      FILENAME = FILE_LIST(PAGE,FILENUM)
!      CALL SHOWPICK(FILENAME,LEN(FILENAME))
76010 FORMAT (10X,'<',I2,'> ',A60)
76011 FORMAT (10X,'Enter a file number or quit <Q,q>: ' )
76012 FORMAT (10X,'Enter a file number, see next page <N,n>'            &
     &  ,', or quit <Q,q>: ') 
76013 FORMAT (10X,'Enter a file number, see next page <N,n>'            &
     &  ,', see preivous page <P,p>, or quit <Q,q>: ')  
76015 FORMAT (10X,'Enter a file number',                                &
     &   ', see preivous page <P,p>, or quit <Q,q>: ')  
      END SUBROUTINE PICKBIG
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE PICKDEG(FILENAME)

      CHARACTER(2)::SELECTION,GOON
      CHARACTER(3)::FILE_EXT
      CHARACTER(60)::MESSAGE, FILE_LIST(5,16)
      CHARACTER(52)::FILENAME 

      INTEGER(4)::FILENUM,SELECTED,PAGE,I,J,M,N,FNUM,CHAR2INT
! .INI read block
      CHARACTER(10)::INIVAR
      CHARACTER(12)::INIVAL
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
      INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR
!      
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN ! Win/DOS syntax
         CALL SYSTEM('DIR /B *.DEG > FILELIST.TXT')
      ELSE 
         CALL SYSTEM('ls -1 *.DEG > FILELIST.TXT')
      ENDIF
      
      OPEN(UNIT=60,FILE='FILELIST.TXT',STATUS='OLD')
     
      N=-1
      DO I=1,5
         DO J = 1,16
            N=N+1
            READ(60,*,ERR=77014,END=77014) FILE_LIST(I,J)
         ENDDO
      ENDDO
77014 CONTINUE
      CLOSE(60)  
      SELECT CASE (N)
      CASE (0)        ! No such files to view
         MESSAGE = ' No such files exist in this directory!'
         CALL OOPS(MESSAGE,39)
         RETURN 
      CASE (1:16)     ! Show what there is  
         MESSAGE = ' '
      CASE (17:80)    ! Start with the first page 
         MESSAGE = ' '
      CASE DEFAULT    ! Start with the first page, but there are more 
         MESSAGE = ' '  !    than 80 files to choose from 
      END SELECT 

      M=16 !number of files to show on a page
      PAGE=1
77009 CALL CLS
         WRITE(*,*) ' '
         PRINT*,' PAGE ',PAGE,' OF ',(N+16)/16  
         WRITE(*,*) ' '
         WRITE(*,*)'    Showing available files 16 at a time.'
         WRITE(*,*) ' '
! SHOW FILES 1 PAGE AT  A TIME
         IF (M>N-(PAGE-1)*16) THEN ! DO NOT HAVE A FULL PAGE 
           M=N
         ELSE
           M=16
         ENDIF
         DO I = 1, M ! 
           WRITE(*,77010) I+(PAGE-1)*16, FILE_LIST(PAGE,I)
         ENDDO
         PRINT *,' '
         SELECT CASE(PAGE)
             CASE(1)
                IF (N.LT.17) THEN 
                    WRITE(*,77011)
                 ELSE 
                    WRITE(*,77012) 
               ENDIF            
              CASE(2:4)
                    WRITE(*,77013) 
              CASE(5)
                    WRITE(*,77015) 
           END SELECT    
         CALL READKB(SELECTION,2)
      IF (SELECTION=='n' .OR. SELECTION=='N') THEN
         PAGE = PAGE + 1
         GOTO 77009
      ENDIF
      IF (SELECTION=='P' .OR. SELECTION=='p') THEN
         PAGE = PAGE - 1
         GOTO 77009
      ENDIF
      IF (SELECTION=='q' .OR. SELECTION=='Q') CALL MENU_MAIN
      
      FILENUM=CHAR2INT(SELECTION,2)
      FILENAME = FILE_LIST(PAGE,FILENUM)
!      CALL SHOWPICK(FILENAME,LEN(FILENAME))
77010 FORMAT (10X,'<',I2,'> ',A60)
77011 FORMAT (10X,'Enter a file number or quit <Q,q>: ' )
77012 FORMAT (10X,'Enter a file number, see next page <N,n>'            &
     &  ,', or quit <Q,q>: ') 
77013 FORMAT (10X,'Enter a file number, see next page <N,n>'            &
     &  ,', see preivous page <P,p>, or quit <Q,q>: ')  
77015 FORMAT (10X,'Enter a file number',                                &
     &   ', see preivous page <P,p>, or quit <Q,q>: ')  
      END SUBROUTINE PICKDEG
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE PICKLOC(FILENAME)

      CHARACTER(2)::SELECTION,GOON
      CHARACTER(3)::FILE_EXT
      CHARACTER(60)::MESSAGE, FILE_LIST(5,16)
      CHARACTER(52)::FILENAME 

      INTEGER(4)::FILENUM,SELECTED,PAGE,I,J,M,N,FNUM,CHAR2INT
! .INI read block
      CHARACTER(10)::INIVAR
      CHARACTER(12)::INIVAL
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
      INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR
!      
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN ! Win/DOS syntax
         CALL SYSTEM('DIR /B *.LOC > FILELIST.TXT')
      ELSE 
         CALL SYSTEM('ls -1 *.LOC > FILELIST.TXT')
      ENDIF
      
      OPEN(UNIT=60,FILE='FILELIST.TXT',STATUS='OLD')
     
      N=-1
      DO I=1,5
         DO J = 1,16
            N=N+1
            READ(60,*,ERR=78014,END=78014) FILE_LIST(I,J)
         ENDDO
      ENDDO
78014 CONTINUE
      CLOSE(60)  
      SELECT CASE (N)
      CASE (0)        ! No such files to view
         MESSAGE = ' No such files exist in this directory!'
         CALL OOPS(MESSAGE,39)
         RETURN 
      CASE (1:16)     ! Show what there is  
         MESSAGE = ' '
      CASE (17:80)    ! Start with the first page 
         MESSAGE = ' '
      CASE DEFAULT    ! Start with the first page, but there are more 
         MESSAGE = ' '  !    than 80 files to choose from 
      END SELECT 

      M=16 !number of files to show on a page
      PAGE=1
78009 CALL CLS
         WRITE(*,*) ' '
         PRINT*,' PAGE ',PAGE,' OF ',(N+16)/16  
         WRITE(*,*) ' '
         WRITE(*,*)'    Showing available files 16 at a time.'
         WRITE(*,*) ' '
! SHOW FILES 1 PAGE AT  A TIME
         IF (M>N-(PAGE-1)*16) THEN ! DO NOT HAVE A FULL PAGE 
           M=N
         ELSE
           M=16
         ENDIF
         DO I = 1, M ! 
           WRITE(*,78010) I+(PAGE-1)*16, FILE_LIST(PAGE,I)
         ENDDO
         PRINT *,' '
         SELECT CASE(PAGE)
             CASE(1)
                IF (N.LT.17) THEN 
                    WRITE(*,78011)
                 ELSE 
                    WRITE(*,78012) 
               ENDIF            
              CASE(2:4)
                    WRITE(*,78013) 
              CASE(5)
                    WRITE(*,78015) 
           END SELECT    
         CALL READKB(SELECTION,2)
      IF (SELECTION=='n' .OR. SELECTION=='N') THEN
         PAGE = PAGE + 1
         GOTO 78009
      ENDIF
      IF (SELECTION=='P' .OR. SELECTION=='p') THEN
         PAGE = PAGE - 1
         GOTO 78009
      ENDIF
      IF (SELECTION=='q' .OR. SELECTION=='Q') CALL MENU_MAIN
      
      FILENUM=CHAR2INT(SELECTION,2)
      FILENAME = FILE_LIST(PAGE,FILENUM)
!      CALL SHOWPICK(FILENAME,LEN(FILENAME))
78010 FORMAT (10X,'<',I2,'> ',A60)
78011 FORMAT (10X,'Enter a file number or quit <Q,q>: ' )
78012 FORMAT (10X,'Enter a file number, see next page <N,n>'            &
     &  ,', or quit <Q,q>: ') 
78013 FORMAT (10X,'Enter a file number, see next page <N,n>'            &
     &  ,', see preivous page <P,p>, or quit <Q,q>: ')  
78015 FORMAT (10X,'Enter a file number',                                &
     &   ', see preivous page <P,p>, or quit <Q,q>: ')  
      END SUBROUTINE PICKLOC
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE SHOWPICK(FILENAME,LOFN)
! OPENS 'FILENAME' WITH DEFAULT VIEWER
     
      CHARACTER(80)::COMMAND
      CHARACTER(LOFN)::FILENAME
      CHARACTER(12)::VIEWER,VARNAME
      CHARACTER(1)::TRASH

      INTEGER(4)::I,LOFN

      I=0
! GET DEFAULT VIEWER FROM FILE
      OPEN (UNIT=99,FILE='CARI.INI',STATUS='OLD')
      DO WHILE (I.EQ.0)             
         READ(99,*) VARNAME, TRASH, VIEWER
         IF (VARNAME=='VIEWER') THEN 
            I=1
         ELSEIF (VARNAME=='END') THEN 
            I=2
         ELSE
            I=0
         ENDIF  
      ENDDO     
76100 CLOSE(99)
      IF (I==1) THEN
         COMMAND = VIEWER//' '//FILENAME
      ELSE
         COMMAND = 'NOTEPAD '//FILENAME
      ENDIF 
      CALL SYSTEM(COMMAND)
      END SUBROUTINE  
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE MENUHEADER
      PRINT*,' ' 
      PRINT*,                                                           &
     &'  CARI-7                       Civil Aerospace Medical Institute'
      PRINT*,                                                           &
     &'  April 13, 2017                 Federal Aviation Administration'
      PRINT*,''
      PRINT*,''

      END SUBROUTINE      
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE READKB(STROUT,L)
! SPECIAL READ FUNCTION FOR USER INPUT FROM KEYBOARD
! KC, 7 JUNE 2012
      INTEGER(4)::L
      CHARACTER(L)::STRIN,STROUT
      CHARACTER(3)::DIAGNOSE

      DIAGNOSE='NO!'
      IF (L.GT.1) PRINT*,'(If entering more than a word, enclose answer'&
     &,' in quotes)' 
      READ*,STRIN
      CALL LC2UC(STRIN,STROUT,L)    

      END SUBROUTINE      
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE LC2UC(STRIN,STROUT,L)
! CONVERTS LOWERCASE LETTERS IN STRINGS TO UPPERCASE
! KC, 7 JUNE 2012
      INTEGER(4)::L,I
      CHARACTER(L)::STRIN,STROUT
      CHARACTER(1)::A
      CHARACTER(3)::DIAGNOSE
      DIAGNOSE='NO!'

      DO I = 1,L
         A=STRIN(I:I)
         IF (IACHAR(A) .GT. 96 .AND. IACHAR(A).LT. 123) THEN
            ! SHIFT IN ASCII DOWN 32 TO GET UPPERCASE FROM LOWERCASE
            STROUT(I:I)=CHAR(IACHAR(A)-32)
         ELSE  
            STROUT(I:I)=STRIN(I:I)
         ENDIF
      ENDDO
      IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'CONVERTED ',STRIN,' TO ',STROUT
      END SUBROUTINE      
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION DKSTR(I)
            INTEGER(4)::I
        CHARACTER(39)::DKSTR
          SELECT CASE (I)    
             CASE(1)
                 DKSTR=' p/sq.cm, secondary particle fluence   '
             CASE(2)
                 DKSTR=' microSv, ICRP Pub. 103 EFFECTIVE DOSE '
             CASE(3)
                 DKSTR=' microSv, ICRP Pub. 60 EFFECTIVE DOSE  '
             CASE(4)
                 DKSTR=' microSv, ICRU AMBIENT DOSE EQ. H(*10) '
             CASE(5)
                 DKSTR=' microGy, AVE WHOLE-BODY ABSORBED DOSE '
          END SELECT
      END FUNCTION DKSTR
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION DRKSTR(I)
            INTEGER(4)::I
        CHARACTER(42)::DRKSTR
          SELECT CASE (I)    
             CASE(1)
                 DRKSTR=' particles/sq.cm/sec, SECONDARY FLUX      '
             CASE(2)
                 DRKSTR=' microSv/hr, ICRP Pub. 103 EFFECTIVE DOSE '
             CASE(3)
                 DRKSTR=' microSv/hr, ICRP Pub. 60 EFFECTIVE DOSE  '
             CASE(4)
                 DRKSTR=' microSv/hr, ICRU AMBIENT DOSE EQ. H(*10) '
             CASE(5)
                 DRKSTR=' microGy/hr, AVE WHOLE-BODY ABSORBED DOSE '
          END SELECT
      END FUNCTION DRKSTR
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION PARTSTR(I)
            INTEGER(4)::I
        CHARACTER(10)::PARTSTR
          SELECT CASE (I)    
             CASE(0)
                PARTSTR='TOTAL     '
             CASE(1)
                PARTSTR='NEUTRONS  '
             CASE(2)
                PARTSTR='PHOTONS   '
             CASE(3)
                PARTSTR='ELECTRONS '
             CASE(4)
                PARTSTR='POSITRONS '
             CASE(5)
                PARTSTR='NEG MUONS '
             CASE(6)
                PARTSTR='POS MUONS '
             CASE(7)
                PARTSTR='PROTONS   '
             CASE(8)
                PARTSTR='NEG PIONS '
             CASE(9)
                PARTSTR='POS PIONS '
             CASE(10)
                PARTSTR='DEUTERONS '
             CASE(11)
                PARTSTR='TRITONS   '
             CASE(12)
                PARTSTR='HELIONS   '
             CASE(13)
                PARTSTR='ALPHAS    '
             CASE(14)
                PARTSTR='LITHIUM   '
             CASE(15)
                PARTSTR='BERYLLIUM '
             CASE(16)
                PARTSTR='BORON     '
             CASE(17)
                PARTSTR='CARBON    '
             CASE(18)
                PARTSTR='NITROGEN  '
             CASE(19)
                PARTSTR='OXYGEN    '
             CASE(20)
                PARTSTR='FLUORINE  '
             CASE(21)
                PARTSTR='NEON      '
             CASE(22)
                PARTSTR='SODIUM    '
             CASE(23)
                PARTSTR='MAGNESIUM '
             CASE(24)
                PARTSTR='ALUMINUM  '
             CASE(25)
                PARTSTR='SILICON   '
             CASE(26)
                PARTSTR='PHOSPHORUS'
             CASE(27)
                PARTSTR='SULPHUR   '
             CASE(28)
                PARTSTR='CHLORINE  '
             CASE(29)
                PARTSTR='ARGON     '
             CASE(30)
                PARTSTR='POTASSIUM '
             CASE(31)
                PARTSTR='CALCIUM   '
             CASE(32)
                PARTSTR='SCANDIUM  '
             CASE(33)
                PARTSTR='TITANIUM  '
             CASE(34)
                PARTSTR='VANADIUM  '
             CASE(35)
                PARTSTR='CHROMIUM  '
             CASE(36)
                PARTSTR='MANGANESE '
             CASE(37)
                PARTSTR='IRON      '
             CASE(38)
                PARTSTR='TOTAL     '
          END SELECT
      END FUNCTION PARTSTR
!----------------------------------------------------------------------

