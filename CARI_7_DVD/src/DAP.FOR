!      program DAP2
!     FUNCTION DAP(lat,lon,alt,t,testrun,gcrmodel,superpos,dosekind)
!DAP.FOR by Kyle Copeland
!  Patched for database misreading and top altitude fluence error in databases 
!  by KC 30OCT2014, see LOADFTDCCS in UTILITY.for for details
!Calculates Dose At Point (lat,lon,alt,t) of kind d using gcrmodel n 
! uses superposition if superpos=1
! Unit 130 for diagnostic output if testrun=1
!      IMPLICIT NONE
!      REAL(8)::DOSE,d
!      INTEGER(4)::testrun,gcrmodel,superpos,dosekind,dkind,yr,mo,da,hr
!      REAL(8)::lat,lon,alt,t
!      CHARACTER*3::diagnose
! toggle diagnostic print and write statements
!      diagnose="YES" 
! function inputs from CARI-7
!      lat=0. !North is positive
!      lon=0. !East is positive
!      alt=100. !Depth in g/cm2
!      t=1996.380 !a Solar min, May 1996, 372 MV
!      gcrmodel=1 !1 or >3 is ISO/Nymmik 2013, 2 is Badhwar and O'Neill 2O11, 3 is O'Brien 2000
!      testrun=0  !TF on testing, 0 for production use in CARI   
!      superpos=0 ! TF on superposition approximation. 
                 ! 0 = use HZEs transport in atmosphere, 1 = no HZE transport in atmosphere, uses
                 ! only p and n nucleons for primaries after geomagnetic and solar filters 
    
!      DO dkind=76,190,38
!         dkind=40,80,120,160,200 ! which dose to calculate?
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
! 38 = ICRP Pub 103 total effective dose
! 38 INACTIVE = ICRP Pub. 103 Effective dose from Co ions 
! 39 INACTIVE = ICRP Pub. 103 Effective dose from Ni ions
!
!But Since 38 and 39 are inactive
! 1-38 = secondary particle fluxes (38 is TOTAL ION FLUX)
! 39-76 = ICRP Pub. 103 Total Effective dose from all available ions
! 77-114 =  Same but  for ICRP Pub 60 Effective dose
! 115-152 = "   "    "  ICRU ambient dose equivalent H*(10)
! 153-190 = "   "    "  whole body average absorbed dose
   
! Get result         
!      d=DOSE(lat,lon,alt,t,dkind,gcrmodel,testrun,superpos)
!      IF (diagnose.eq.'YES') then
!         Write(*,9999) lat,lon,alt,t,d
!         Write(*,*) ''
!         call dates(t,yr,mo,da,hr)
!         WRITE(*,9998)'Date is:',yr,mo,da,' at ',hr,' h' 
!      ENDIF
!      STOP !diagnostic
!      ENDDO
!9998  FORMAT(A8,3I6,A4,I4,A2)
!9999  FORMAT('The dose at ',f6.2,' N,',f6.2,' E,',f9.4,' g/cm2 for date'&
!     & ,f11.5,' is ',es10.3)
!      end program DAP2
!     DAP=d 
!     END FUNCTION DAP
!----------------------------------------------------------------------&

      function DOSE(la,lo,al,t,dosekind,gcrmodel,testrun,superpos)
      IMPLICIT NONE
      INTEGER, PARAMETER::NZ = 1, NA = 1
      REAL(8)::t,la,lo,al,d,dd,DOSE,NVALT,SLANT,S,NMF
! Variance variables and functions
      REAL(8)::LVAR,PVAR,NVAR,VS,stddev
      REAL(8)::SHOWS,nSHOWS,pSHOWS,PHIVAR,FTDVAR,PFVAR,NFVAR
! Altitude grid
      REAL(8),DIMENSION(19)::depth
      Real(8)::VC, norm, FTD, ftdcc, pftdcc, nftdcc, FHF, hiflux
      Real(8)::fluxmod
! Vertical cutoff rigidity in GV
! normalization constants
      REAL(8)::j1,j2,j3,j4,j5,j7,j8
! dummy reads
!      REAL(8), DIMENSION(28,19)::flux,R,RU,RL,KE,KEPN
! There are 29 primary particles:n+ion species 1-28
!   (Ni and Co, while counted, are currently not in use)
      REAL(8), DIMENSION(27,100)::gcrflux,gcrKE,gcrKEPN,gcrR
      REAL(8), DIMENSION(27,100)::flux,R,UR,LR,KE,KEPN,GCRVAR
! There are 29 primary particles:ion species 1-28
!   (Ni and Co, while counted, are currently not in use)
! There are 19 primary kinetic energy bins 1-19 
! flux(z,e) = differential flux in particles/(m**2-sr-sec-GeV)
! RL,RU, and R(z,e) = upper lower and characteristic rigidity of 
!                     primary type z in energy bin e
! KE = Characteristic energy of primary type z in energy bin e
! KEPN = Characteristic energy/nucleon of primary type z in energy bin e
!      REAL(8), DIMENSION(27,19)::pflux,nflux
! pflux, nflux = proton and neutron fluxes in superposition approximation
      REAL(8), DIMENSION(100)::egrid
! energy grid values in GeV
      REAL(8), DIMENSION(19)::deltaKE,GeVgrid
! UE, LE = upper and lower energies of each E band, in GeV
      REAL(8)::Sky, LL, Kinetic, UE, LE, km, fbf, Zetastep, zstep
      REAL(8), DIMENSION(NZ,NA)::Skycut,Skywt,Skypass
! The sky around the aircraft is divided into NA x NZ 
!     sectors for evaluation of approaching cosmic rays primaries 
! Skycut = rigidity cutoff for each sector of the sky in GV
! Skywt = fractional area of that sector of the sky relative to 
!         the celestial sphere around the craft 
! Skypass = how much view of this sector is blocked or 'occulted' by the 
!         Earth (assumes no significant contribution from primaries 
!         approaching from below horizon)   
      INTEGER(4)::dosekind,testrun,gcrmodel,superpos,e,es,z
!     t = date and time converted to digital yyyy.XXXXX
!     la = E latitude
!     lo= N longitude
!     al = altitude in g/cm**2
!     dk = dose kind (2ndary particle to calculate dose or flux from)
!     test = 1 or 0; write selected diagnostic prints
!     sp = 1 or 0; use superposition for protons primaries
!     e = 19 primary energies
!     es = 100 pt spectral energy grids
!     z = 1-28, nuclear charge, a.k.a. atomic number
!         NOTE: currently doses for 1-26 are available
      CHARACTER*72::header,HEADER2
      CHARACTER*48::HEADER1
      CHARACTER*8::SUBDIR
      CHARACTER*12, DIMENSION(29)::FTDCCF
      INTEGER(4), DIMENSION(28)::N
!     N(z)=nucleon number
      Real(8)::A, m, Te, pi
      DIMENSION::A(28), m(28)
      INTEGER(4)::IZ,IA,IE,IS,i,ii,nothere
      REAL(8)::pdose,ndose,afluence,pfluence,nfluence,pftd,nftd,Rigidity
      INTEGER(4), DIMENSION(28)::zfilter
! zfilter = pass band matrix for which ions of the primary spectrum are allowed         
      REAL(8), DIMENSION(27,19,19,38)::PHI,E103,E60,H10,GRAY
      REAL(8), DIMENSION(27,19,19,38)::PHIS,E103S,E60S,H10S,GRAYS
      REAL(8), DIMENSION(100)::X,Y
!     ENERGY grid for MCNPX calculated coefficients conversion
      REAL(8), DIMENSION(19)::CEGRID  
!
      real :: start, calcstart, endgcr, endfb, finish, rj1
!
      LOGICAL::WV
!
      COMMON /energy/egrid,CEGRID
      COMMON /kepn/gcrKEPN
      COMMON /phys/depth 
      COMMON /atomic/N,A,m 
      COMMON /ftdcc/PHI,E103,E60,H10,GRAY
      COMMON /ftdccs/PHIS,E103S,E60S,H10S,GRAYS
      COMMON /ISOSTATS/stddev
      COMMON /VCUT/VC
! CARI.INI VARIABLES
      CHARACTER(10)::INIVAR
      CHARACTER(12)::INIVAL
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
! END CARI.INI BLOCK        
!
! Average atomic masses/nucleon numbers to nearest integer for each primary 
! ion type    
      DATA N/1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,            &
     & 40,39,40,45,48,51,52,55,56,59,59/
!THIS zfilter allows only gcr protons to enter the atmosphere
!      DATA zfilter/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,                   &
!     & 0,0,0,0,0,0,0,0,0,0,0/
!THIS zfilter allows only gcr alphas to enter the atmosphere
!      DATA zfilter/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,                   &
!     & 0,0,0,0,0,0,0,0,0,0,0/
!THIS zfilter allows all gcr ions to enter the atmosphere
      DATA zfilter/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,                   &
     & 1,1,1,1,1,1,1,1,1,1,1/
!THIS zfilter allows only HZE gcr ions to enter the atmosphere
!      DATA zfilter/0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,                   &
!     & 1,1,1,1,1,1,1,1,1,1,1/
! energy bin data
      DATA (deltaKE(e),e=1,19)/0.001,0.001,0.005,0.005,0.015,0.045,0.055&
     & ,0.145,0.455,0.545,1.455,4.545,5.455,14.545,45.455,54.545,145.455&
     & ,454.545,545.455/ !GeV
      DATA CEGRID/1.,2.,5.,10.,20.,50.,100.,200.,500.,1000.,2000.,5000.,&
     & 10000.,20000.,50000.,100000.,200000.,500000.,1000000./ !MeV
! GCR spectrum energy/nucleon grid
      EGRID(1)=1.000E-02 ! 1.929E+02
      EGRID(2)=1.150E-02 ! 2.203E+02
      EGRID(3)=1.322E-02 ! 2.513E+02
      EGRID(4)=1.520E-02 ! 2.863E+02
      EGRID(5)=1.748E-02 ! 3.259E+02
      EGRID(6)=2.009E-02 ! 3.703E+02
      EGRID(7)=2.310E-02 ! 4.201E+02
      EGRID(8)=2.656E-02 ! 4.757E+02
      EGRID(9)=3.054E-02 ! 5.375E+02
      EGRID(10)=3.511E-02!  6.058E+02
      EGRID(11)=4.037E-02!  6.809E+02
      EGRID(12)=4.642E-02!  7.630E+02
      EGRID(13)=5.337E-02!  8.518E+02
      EGRID(14)=6.136E-02!  9.472E+02
      EGRID(15)=7.055E-02!  1.049E+03
      EGRID(16)=8.111E-02!  1.155E+03
      EGRID(17)=9.326E-02!  1.265E+03
      EGRID(18)=1.072E-01!  1.377E+03
      EGRID(19)=1.233E-01!  1.489E+03
      EGRID(20)=1.417E-01!  1.597E+03
      EGRID(21)=1.630E-01!  1.700E+03
      EGRID(22)=1.874E-01!  1.792E+03
      EGRID(23)=2.154E-01!  1.871E+03
      EGRID(24)=2.477E-01!  1.933E+03
      EGRID(25)=2.848E-01!  1.974E+03
      EGRID(26)=3.275E-01!  1.992E+03
      EGRID(27)=3.765E-01!  1.983E+03
      EGRID(28)=4.329E-01!  1.948E+03
      EGRID(29)=4.977E-01!  1.886E+03
      EGRID(30)=5.722E-01!  1.798E+03
      EGRID(31)=6.579E-01!  1.687E+03
      EGRID(32)=7.565E-01!  1.558E+03
      EGRID(33)=8.697E-01!  1.415E+03
      EGRID(34)=1.000E+00!  1.263E+03
      EGRID(35)=1.150E+00!  1.108E+03
      EGRID(36)=1.322E+00!  9.556E+02
      EGRID(37)=1.520E+00!  8.099E+02
      EGRID(38)=1.748E+00!  6.745E+02
      EGRID(39)=2.009E+00!  5.529E+02
      EGRID(40)=2.310E+00!  4.455E+02
      EGRID(41)=2.656E+00!  3.533E+02
      EGRID(42)=3.054E+00!  2.759E+02
      EGRID(43)=3.511E+00!  2.124E+02
      EGRID(44)=4.037E+00!  1.613E+02
      EGRID(45)=4.642E+00!  1.209E+02
      EGRID(46)=5.337E+00!  8.964E+01
      EGRID(47)=6.136E+00!  6.576E+01
      EGRID(48)=7.055E+00!  4.778E+01
      EGRID(49)=8.111E+00!  3.443E+01
      EGRID(50)=9.326E+00!  2.464E+01
      EGRID(51)=1.072E+01!  2.045E+01
      EGRID(52)=1.233E+01!  1.577E+01
      EGRID(53)=1.417E+01!  1.210E+01
      EGRID(54)=1.630E+01!  9.202E+00
      EGRID(55)=1.874E+01!  6.968E+00
      EGRID(56)=2.154E+01!  5.250E+00
      EGRID(57)=2.477E+01!  3.934E+00
      EGRID(58)=2.848E+01!  2.935E+00
      EGRID(59)=3.275E+01!  2.181E+00
      EGRID(60)=3.765E+01!  1.615E+00
      EGRID(61)=4.329E+01!  1.192E+00
      EGRID(62)=4.977E+01!  8.771E-01
      EGRID(63)=5.722E+01!  6.436E-01
      EGRID(64)=6.579E+01!  4.710E-01
      EGRID(65)=7.565E+01!  3.438E-01
      EGRID(66)=8.697E+01!  2.505E-01
      EGRID(67)=1.000E+02!  1.821E-01
      EGRID(68)=1.150E+02!  1.321E-01
      EGRID(69)=1.322E+02!  9.574E-02
      EGRID(70)=1.520E+02!  6.925E-02
      EGRID(71)=1.748E+02!  4.998E-02
      EGRID(72)=2.009E+02!  3.607E-02
      EGRID(73)=2.310E+02!  2.597E-02
      EGRID(74)=2.656E+02!  1.868E-02
      EGRID(75)=3.054E+02!  1.341E-02
      EGRID(76)=3.511E+02!  9.625E-03
      EGRID(77)=4.037E+02!  6.896E-03
      EGRID(78)=4.642E+02!  4.935E-03
      EGRID(79)=5.337E+02!  3.529E-03
      EGRID(80)=6.136E+02!  2.521E-03
      EGRID(81)=7.055E+02!  1.798E-03
      EGRID(82)=8.111E+02!  1.282E-03
      EGRID(83)=9.326E+02!  9.129E-04
      EGRID(84)=1.072E+03!  6.498E-04
      EGRID(85)=1.233E+03!  4.614E-04
      EGRID(86)=1.417E+03!  3.280E-04
      EGRID(87)=1.630E+03!  2.324E-04
      EGRID(88)=1.874E+03!  1.647E-04
      EGRID(89)=2.154E+03!  1.167E-04
      EGRID(90)=2.477E+03!  8.254E-05
      EGRID(91)=2.848E+03!  5.834E-05
      EGRID(92)=3.275E+03!  4.118E-05
      EGRID(93)=3.765E+03!  2.906E-05
      EGRID(94)=4.329E+03!  2.049E-05
      EGRID(95)=4.977E+03!  1.443E-05
      EGRID(96)=5.722E+03!  1.016E-05
      EGRID(97)=6.579E+03!  7.140E-06
      EGRID(98)=7.565E+03!  5.014E-06
      EGRID(99)=8.697E+03!  3.520E-06
      EGRID(100)=1.000E+04!  2.468E-06
! Allows testing of effect of small shifts in EGRID    
!     EGRID=EGRID*0.95 
! Altitude data: MCNPX tally surface depths in g/cm**2
      DATA depth/1035.0837,918.3537,716.8487,483.0577,315.3317,198.7657,&
     & 124.2157,77.6617,48.5928,30.5677,19.3575,12.3383,7.9173,3.3999,  &
     & 0.7309,0.1508,0.02462,0.002941,0.0/ 
! Note: Since 0 depth data does not include the all GCR flux, only backscatter  
! and shallow incidence angle outbound crossers, a correction is applied at
! at 0.
 
! GeV ATOMIC MASSES PER NUCLEON 
      Data m/0.938,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939&
     &,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939&
     &,0.939,0.939,0.939,0.939,0.939,0.939,0.939/

! ATOMIC WEIGHTS
      Data A/1.,4.,6.9,9.,10.8,12.,14.,16.,19.,20.2,23.,24.3,27.,28.1,  &
     & 31.,32.1,35.4,39.9,39.1,40.1,44.9,47.9,50.9,52.,54.9,55.8,58.9,  &
     & 58.7/

! GENERAL NORMALIZATION
! norm converts isotropic fluence to fluence entering the atmosphere and 
!      changes units: 
!        0.5 since only the downward fluence matters and MCNPX sim is normalized 
!            to per incident particle,
!        A 2/pi could be considered to approximately undo MCNPX surface flux 
!            normalization of |sec(x)|, however, exmination of point detectector 
!            fluxes vs shells indicates this is not required. Fluences are 
!            consistently within 2x of reported std errors.  
!        4pi steradians in unit sphere, 3600 sec/hour, 10000 cm**2/m**2
      pi=3.1415927
      norm =0.5*(4*pi)*3600./10000. 
! NMF (neutron flux indicated magnetic field induced spreading) corrects approximately 
! for the lack of magnetic field effect presence in MCNPX on the transport path 
! length by assigning a slightly longer path length to shower depths.
! While path lengths really remain unchanged in a magnetic field, the field broadens 
! charged particle paths. MCNP6 results for GCR protons are equivalent to increasing the 
! proton->neutron cross section by a few percent - Paul Goldhagen conversation. 
! Since the MCNPX transport did not include magnetic fields, charge particle
! showers are overly at some angles of incidence penetrating. This is a simple 
! approximation of the resultant influence, since path length bending is based on 
! rigidity and orientation relative to the local field. 
      NMF= 1.0! 1.00 to turn off, ! Estimated 3-5% for protons.

! Set up GeV gridpoints
      DO e=1,19 
         GeVgrid(e)=cegrid(e)/1000. 
      ENDDO

      call cpu_time(start)
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
         OPEN(unit=130,file='DIAGNOSE\DIAGDAP.DAT',status='UNKNOWN')
      ELSE
         OPEN(unit=130,file='DIAGNOSE/DIAGDAP.DAT',status='UNKNOWN')
      ENDIF 
!      WRITE(*,*)' Starting ', la,lo,al,t,dosekind,gcrmodel,testrun,     &
!     &superpos
      
      if (testrun.eq.1)then
         WRITE(130,*)' Starting ', la,lo,al,t,dosekind,gcrmodel,testrun, &
     &superpos
         write(*,*) 'Loading fluence conversion coefficients' 
         write(130,*) 'Loading fluence conversion coefficients' 
      endif

!      CALL LOADFTDCCS ! this call moved to CARI to speed up calculations
!                        as it only needs to be done once.

!      Bring converted West longitudes from CARI-7 into the 
!        expected 0-360 East range used in the VCR tables
       IF (lo.GT.360.0) THEN
          DO WHILE (lo.GT.360.0)
             lo=lo-360.0
          ENDDO
       ELSEIF (lo.LT.0.0) THEN
          DO WHILE (lo.LT.0.0)
             lo=lo+360.
          ENDDO
       ELSE
          lo=lo
       ENDIF

! CALCULATE CR particle spectra for date
! All spectra *.out files have the same format, 'i4,2es11.3',
! units are E in GeV, phi in #/(cm**2 s sr GeV/n) 
! BUT different numbers of header rows
      CALL cpu_time(calcstart)
      WRITE(40,*)' PreGCR time in DOSE in seconds. ',calcstart-start
! unit 40 is diagnostic for CARI-7
      CALL GETGCR(4,t,OS)
! initialize arrays
      gcrflux=0.0
      gcrKE=0.0
      GCRVAR=0.0
!
      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
      ! USE WINDOWS/DOS PATHNAMES
         OPEN(unit=50,file='GCR_MODELS\ISO_GCR2.OUT',status='old')
         READ(50,*) j1
         READ(50,50011)header
         DO z = 1,26
            DO e=1,100
               READ(50,50012) j1,gcrKE(z,e),gcrflux(z,e)
               gcrKEPN(z,e)=gcrKE(z,e)/N(z)
               if (testrun.EQ.1) then
                  WRITE(130,*) z, gcrKE(z,e), gcrflux(z,e)
               endif
            ENDDO
         ENDDO 
         CLOSE(50)
      ELSE 
! REPEAT BUT USE POSIX/LINUX FILESYSTEM PATH NAMES
        OPEN(unit=50,file='GCR_MODELS/ISO_GCR2.OUT',status='old')
        READ(50,*) j1
        READ(50,50011)header
        DO z = 1,26
           DO e=1,100
              READ(50,50012) j1,gcrKE(z,e),gcrflux(z,e)
              gcrKEPN(z,e)=gcrKE(z,e)/N(z)
              if (testrun.EQ.1) then
                 WRITE(130,*) z, gcrKE(z,e), gcrflux(z,e)
              endif
           ENDDO
        ENDDO 
        CLOSE(50)
      ENDIF !END CALLING FOR SPECTRA
! KC 20161205
! Added the ISO GCR model driven by heliocentric potential instead of susnspots
! 
      CALL cpu_time(endgcr)
      WRITE(40,*)' GCR time in seconds. ',endgcr-calcstart

      if (testrun.eq.1)then
         WRITE(130,*)' GCR time in seconds. ',endgcr-calcstart
      endif
50001 FORMAT(I5,ES12.3,6ES11.3)
50011 FORMAT(A72)
50012 FORMAT(I4,2ES11.3)
50014 FORMAT(I4,4ES11.3)
50013 FORMAT(I4,ES11.3)
      CALL findcut(NA,NZ,t,la,lo,km(al),VC,Skycut,Skypass,Skywt)
        If (testrun.eq.1) then
          write(130,*) 'Geophysical and magnetic filters built'
          WRITE(130,*) ' Vertical cutoff rigidity = ', VC
          WRITE(130,*) ' SkyCut Skypass Skyweight'
          DO i=1,NZ
            DO ii=1,NA
            WRITE(130,*) SkyCut(i,ii), Skypass(i,ii), Skywt(i,ii)
            enddo
          enddo
!          Write(*,*)'Geophysical and magnetic filters built'
        ENDIF
! Find the zenith angle step size used for non-vertical related adjustments
      Zetastep = zstep(NZ,km(al))
! Forbush effects
      CALL FORBUSH(t,gcrflux,fbf)      
          if (testrun.eq.1) then
             write(130,*)' GCR Adjusted by ', fbf, 'for Forbush effects'
!
!             write(*,*)' For GCR at ',t,' Forbush Adjustment is', fbf
          ENDIF
! initialize DOSES AND VARIANCES
      d=0.0
      dd=0.0
      PVAR=0.0
      NVAR=0.0
      LVAR=0.0
      VS=0.0
! Finally, the Main Loops
      CALL cpu_time(endfb)
      WRITE(40,*)' FORBUSH time in seconds. ',endfb-endgcr
! run spectrum through the sky filter
      IF (testrun.eq.1) then
         write(130,*)'Calculating Dose Rate'
         write(*,*)'Calculating Dose Rate'
      ENDIF
      DO z=1,26 
        DO e=1,100
          PHIVAR=GCRVAR(z,e)
          LR(z,e)=Rigidity(LE(z,e)/N(z),z)
          UR(z,e)=Rigidity(UE(z,e)/N(z),z)
          hiflux=FHF(gcrflux,gcrke,z,e)
          DO i=1,NZ !zenith dependence
! function SKYLIBS/nvalt(Zetastep,alt,i,z) and SLANT(Zetastep,ALT,i,z) 
! control angular depth dependence 
          ftdcc=FTD(z,e,dosekind,nvalt(Zetastep,al*NMF,i,z))
!          nftdcc=nFTD(z,e,dosekind,nvalt(Zetastep,al*NMF,i,z))! KC trial proton replacement
          FTDVAR=SHOWS(z,e,dosekind,nvalt(Zetastep,al*NMF,i,z))**2 !recall VAR=SIGMA**2
          S=SLANT(Zetastep,al*NMF,i,z)
            DO ii=1,NA !azimuthal dependence + numerical integration
              Sky=Skypass(i,ii)*Skywt(i,ii)!horizon cut and geometric weight 
              IF (Sky.NE.0.0) then
                if (LR(z,e).GE.Skycut(i,ii)) then !everything gets through geomag cutoff
                  afluence=hiflux*gcrflux(z,e)*(UE(z,e)-LE(z,e))*norm*S
! KC proton trial cross section fix,convert ~3.5% of protons to neutrons right away
!                IF (z.EQ.1) Then
!                   nfluence=0.0 !afluence*(0.05)*norm*S
!                   ndose=nftdcc*nfluence*Sky*zfilter(z)
!                   afluence=afluence!*0.96 
!                ELSE
!                   ndose=0.0
 !               ENDIF
! end of trial proton cross section fix
                  dd=ftdcc*afluence*Sky*zfilter(z)!+ndose 
                  d=d+dd
                  IF ((gcrflux(z,e).eq.0.0).or.(ftdcc.eq.0.0)) then
                    LVAR=0
                  ELSE
                    LVAR=PHIVAR*(dd/gcrflux(z,e))**2 +                  &
     &                         FTDVAR*(dd/ftdcc)**2
                  ENDIF
                  VS=VS+LVAR
                else !could be partial passage or no passage
                  if (UR(z,e).GT.Skycut(i,ii)) then !some will get through geomag cutoff
                    LL=Kinetic(z,Skycut(i,ii))
                    afluence=hiflux*gcrflux(z,e)*(UE(z,e)-LL)*norm*S  
                    dd=ftdcc*afluence*Sky*zfilter(z)
                    d=d+dd
                    IF ((gcrflux(z,e).eq.0.0).or.(ftdcc.eq.0.0)) then
                      LVAR=0
                    ELSE
                    LVAR=PHIVAR*(dd/gcrflux(z,e))**2 +                  &
     &                   FTDVAR*(dd/ftdcc)**2
                    ENDIF
                    VS=VS+LVAR
                  else
                    dd=0.
                    d=d ! no change
                  endif
                endif
              ELSE
                dd=0.
                d=d ! no change
              ENDIF
            IF (testrun.eq.1) then
!                IF (testrun.NE.1) then
                  write(130,*)'Rigidities', LR(z,e), UR(z,e)
                  write(130,*)'UE:', UE(z,e),' LE:',LE(z,e),' LL:',LL
                  write(130,*)'Skycut',Skycut(i,ii), 'Vcut',VC
                  write(130,*)'Sky ',Sky,' afluence ',afluence
                  write(130,*)'flux',gcrflux(z,e),' norm',norm
                  write(130,*)'dosekind, alt ', dosekind, al
                  write(130,*)'ftdcc ',ftdcc,' dd ',dd,' VS ',VS 
                  write(130,*)'z,e,i,ii',z,e,i,ii
                  write(130,*)'dose=',d
            ENDIF
            ENDDO 
          ENDDO
      if (testrun.eq.1)then
        WRITE(130,13099)'z,e,ftdcc,dose,flux',z,e,ftdcc,d,gcrflux(z,e)
      endif 
        ENDDO 
      ENDDO
13099 FORMAT(A20,2I4,5ES12.4)
      stddev=SQRT(VS)
      IF (testrun.eq.1) THEN
         write(130,*)' Adjusted Dose of type ',dosekind,' is ',d
         write(*,*)' Adjusted Dose of type ',dosekind,' is ',d
      ENDIF
      CALL cpu_time(finish)
      WRITE(40,653)' Finished: ', la,lo,al,t,dosekind,gcrmodel,testrun,    &
     &superpos
      WRITE(40,*)' Adjusted Dose of type ',dosekind,' is ',d
      WRITE(40,*)' Doserate calculation time in seconds.',finish-start
!      WRITE(*,*)' Doserate calculation time in seconds.',finish-start
      CLOSE(130)
      DOSE=d
653   FORMAT(A11,4ES12.4,4I4)
!      STOP ! diagnostic
      end function dose
!                                                                      7
!----6-----------------------------------------------------------------2
!
      SUBROUTINE GETGCR(i,t,OS)
! GETGCR chooses which GCR model spectrum to use
! BO11 is default
      IMPLICIT NONE
      INTEGER(4)::i
      REAL(8)::t,oldt,newt 
      CHARACTER*3::DIAGNOSE
      CHARACTER*5::OS
      LOGICAL::NEWMO
!      DIAGNOSE='YES'
      DIAGNOSE='NO!'

      IF (DIAGNOSE.eq.'YES') THEN
         WRITE(130,*)'BUILDING PRIMARY GCR SPECTRUM FOR DATE ',t  
         WRITE(*,*)'BUILDING PRIMARY GCR SPECTRUM FOR DATE ',t
      ENDIF   

      OPEN(99,file='GCR_MODELS/last_gcr_time.txt',status='UNKNOWN') 
      WRITE(99,'(f11.6)') t
      REWIND(99)
      READ(99,'(f11.6)') newt
      CLOSE(99)     

! SINCE MODELS ARE MONTHLY AVERAGES AND HOURLY RESOLUTION IS ADDED ELSEWHERE, 
! T WITHIN 15 DAYS IS CLOSE ENOUGH TO NOT RECALCULATE GCR   
!     15/365 = 0.041

      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
      !USE DOS/WINDOWS PATHS 
         open(50,file='GCR_MODELS\ISO_GCR2.OUT',status='OLD')
            READ (50,'(f11.6)') oldt
         close(50)
         NEWMO=((oldt-0.041.GT.newt).OR.(oldt+0.041.LT.newt))
         IF (NEWMO) CALL GET_HPD_ISO(t)
      ELSE
!     USE LINUX/POSIX/UNIX PATHS
         open(50,file='GCR_MODELS/ISO_GCR2.OUT',status='OLD')
            READ (50,'(f11.6)') oldt
         close(50)
         NEWMO=((oldt-0.041.GT.newt).OR.(oldt+0.041.LT.newt))
         IF (NEWMO) CALL GET_HPD_ISO(t)
      ENDIF
        
      IF (DIAGNOSE.eq.'YES') THEN
         WRITE(130,*)'PRIMARY GCR SPECTRUM FOR DATE ',t,' FINISHED'  
         WRITE(*,*)'PRIMARY GCR SPECTRUM FOR DATE ',t,' FINISHED'
      ENDIF     
      END SUBROUTINE GETGCR
!                                                                      7
!----6-----------------------------------------------------------------2
! Kinetic Energy from Rigidty
      function Kinetic(Z,R)
      Real(8)::KE,KEP,Kinetic,R
      INTEGER(4)::I,J,Z
      INTEGER(4)::N
      Real(8)::A, m 
      DIMENSION::A(28), m(28), N(28)
!
      COMMON /atomic/N,A,m
!
      KEP=-m(Z)+SQRT(m(Z)**2+(Z*R/A(Z))**2)
      KE=KEP*INT(A(Z)+.2)
      Kinetic=KE
      end function 
!                                                                      7
!----6-----------------------------------------------------------------2
! FTD finds ftdcc at alt for specified z, e, and dosekind
      FUNCTION FTD(Z,E,DCODE,ALT)
      IMPLICIT NONE
      INTEGER(4)::I,J,Z,E,DCODE,DKIND,S
      REAL(8)::ALT,FTD,XA,DD
      REAL(8)::YOUT,XE,RENORM
      CHARACTER*3::DIAGNOSE
      REAL(8),DIMENSION(19)::X,Y,DY,DYY,DEPTH
      REAL(8),DIMENSION(19,19)::DY2
      REAL(8),DIMENSION(27,19)::KEPN      
      REAL(8),DIMENSION(19)::CEGRID
      REAL(8),DIMENSION(100)::EGRID
      REAL(8),DIMENSION(27,19,19,38)::FLUENCE,E103,E60,H10,GRAY
      REAL(8),DIMENSION(27,19,19,38)::FLUENCES,E103S,E60S,H10S,GRAYS
      INTEGER(4)::N
      Real(8)::A, m 
      DIMENSION::A(28), m(28), N(28)
!
      COMMON /atomic/N,A,m
      COMMON /energy/EGRID,CEGRID      
      COMMON /kepn/KEPN
      COMMON /phys/DEPTH 
      COMMON /ftdcc/FLUENCE,E103,E60,H10,GRAY
      COMMON /ftdccs/FLUENCES,E103S,E60S,H10S,GRAYS

!      DIAGNOSE='YES'
      DIAGNOSE='NO!'

      XE=EGRID(E)*N(Z)
! XE is lookup energy in GeV/nucleus

!get ftdccs at all altitudes for specific E from data file and interpolate      
      CALL DECODER(DCODE,DKIND,S)
      IF (DKIND.EQ.1) THEN 
         DO J=1,19
            DO I=1,19
               DY2(I,J)=FLUENCE(Z+1,J,I,S)!recall 1st var is Z+1, protons=2
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.3) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E60(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.4) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=H10(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.5) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=GRAY(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSE ! DEFAULT
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E103(Z+1,J,I,S)
            ENDDO
         ENDDO
      ENDIF

      X=CEGRID*0.001 !*0.001 converts to GeV

      DO I=1,19
         DO J =1,19
            DYY(J)=DY2(I,J)
!            IF (DIAGNOSE.eq.'YES') THEN
!               WRITE(130,*) 'DYY(I)=',DYY(J), 'DY2(I,J)=', DY2(I,J)
!            ENDIF     
         ENDDO
         CALL KCSPLINE(19,X,DYY,1.0e+30,1.0e+30,XE,YOUT)
         IF (DIAGNOSE.eq.'YES') THEN
            DO J=1,19
               WRITE(130,*) 'E(I)=',X(J), ' FTD(I)=', DYY(J)
            ENDDO
            WRITE(130,*)'FTD for E',XE,' at alt ',I,' is ',YOUT  
         ENDIF 
!         IF (S.EQ.1) YOUT=YOUT*1.5 !fixes low neutron fluxes
         DY(I)= YOUT 
      ENDDO
      X=DEPTH
!      WRITE(130,*) 'QUERY ALT: ', ALT
      IF (ALT.LE.1030.) then 
      ! interpolate
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,ALT,YOUT) 
!        DO I=1,19
!           WRITE(130,*) 'ALT(I)=',DEPTH(I), ' FTD(I)=', DY(I)
!        ENDDO
!        WRITE(130,*)'FTD calculated for ALT ',ALT,' is ',YOUT
      ELSE
      XA=1000.
      DD=132. 
! KC 20160918 DD is doubling depth estimate based on LANL report 
!    LA-UR-14-21331 (2014)
!Extrapolation using power law, assume doubling depth of DD g/cm2
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,XA,YOUT)
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' FTD(I)=', DY(I)
          ENDDO
          WRITE(130,*)'FTD for ALT=1000 gpcs is ',YOUT  
        ENDIF
!        WRITE(130,*)'FTD for ', XA, ' gpcs is ',YOUT  
        YOUT=YOUT*0.5**((ALT-XA)/DD)
 !       WRITE(130,*)'FTD extrapolated to ALT ',ALT,' is ',YOUT  
      ENDIF
        
      FTD=YOUT*RENORM(DKIND)
      END FUNCTION FTD 
!                                                                      7
!----6-----------------------------------------------------------------2
! PFTD finds ftdcc for protons at alt for specified z, e, and dosekind
      FUNCTION pFTD(Z,E,DCODE,ALT)
      IMPLICIT NONE
      INTEGER(4)::Z,E,DCODE,DKIND,S
      INTEGER(4)::I,J
      REAL(8)::ALT,pFTD,XA
      REAL(8)::YOUT,RENORM
      CHARACTER*3::DIAGNOSE
      REAL(8),DIMENSION(27,100)::KEPN
      REAL(8),DIMENSION(19)::X,DY,DYY
      REAL(8),DIMENSION(19,19)::DY2
      REAL(8),DIMENSION(19)::DEPTH
      REAL(8),DIMENSION(19)::CEGRID
      REAL(8),DIMENSION(100)::EGRID
      REAL(8),DIMENSION(27,19,19,38)::FLUENCE,E103,E60,H10,GRAY
      REAL(8),DIMENSION(27,19,19,38)::FLUENCES,E103S,E60S,H10S,GRAYS
      INTEGER(4)::N
      Real(8)::A, m 
      DIMENSION::A(28), m(28), N(28)
!
      COMMON /atomic/N,A,m
      COMMON /energy/EGRID,CEGRID      
      COMMON /kepn/KEPN
      COMMON /phys/DEPTH 
      COMMON /ftdcc/FLUENCE,E103,E60,H10,GRAY
      COMMON /ftdccs/FLUENCES,E103S,E60S,H10S,GRAYS

!      DIAGNOSE='YES'
      DIAGNOSE='NO!'

      !get ftdccs at all altitudes from data file and interpolate      
      CALL DECODER(DCODE,DKIND,S)
      IF (DKIND.EQ.1) THEN 
         DO J=1,19
            DO I=1,19
               DY2(I,J)=FLUENCE(2,J,I,S)!recall 1st var is Z+1, protons=2
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.3) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E60(2,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.4) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=H10(2,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.5) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=GRAY(2,J,I,S)
            ENDDO
         ENDDO
      ELSE ! DEFAULT
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E103(2,J,I,S)
            ENDDO
         ENDDO
      ENDIF
      ! create a grid for this KEPN for all altitudes
      X=CEGRID*0.001
      DO I=1,19
         DO J =1,19
            DYY(J)=DY2(I,J)
!            IF (DIAGNOSE.eq.'YES') THEN
!               WRITE(130,*) 'DYY(I)=',DYY(J), 'DY2(I,J)=', DY2(I,J)
!            ENDIF     
         ENDDO
         CALL KCSPLINE(19,X,DYY,1.0e+30,1.0e+30,KEPN(Z,E),YOUT)
         IF (DIAGNOSE.eq.'YES') THEN
            DO J=1,19
               WRITE(130,*) 'E(I)=',X(J), ' pFTD(I)=', DYY(J)
            ENDDO
            WRITE(130,*)'pFTD for KEPN',KEPN(Z,E),' is ',YOUT  
         ENDIF     
         DY(I)=YOUT 
      ENDDO
      X=DEPTH
      IF (ALT.LT.1000) then 
      ! interpolate
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,ALT,YOUT) 
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' pFTD(I)=', DY(I)
          ENDDO
          WRITE(130,*)'pFTD for ALT',ALT,' is ',YOUT  
        ENDIF
      ELSE
! KC 20160918 DD=132 is doubling depth estimate based on LANL report 
!    LA-UR-14-21331 (2014)
!Extrapolation using power law, assume doubling depth of DD g/cm2
      XA=1000.
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,XA,YOUT)
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' pFTD(I)=', DY(I)
          ENDDO
          WRITE(130,*)'pFTD for 1000 gpcs ',ALT,' is ',YOUT  
        ENDIF
        YOUT=YOUT*0.5**((ALT-1000.)/132.)
        WRITE(130,*)'pFTD extrapolated to ALT ',ALT,' is ',YOUT  
      ENDIF   
      pftd=YOUT*RENORM(DKIND)
      END FUNCTION pFTD 
!                                                                      7
!----6-----------------------------------------------------------------2
! NFTD finds ftdcc for neutrons at alt for specified z, e, and dosekind
      FUNCTION nFTD(Z,E,DCODE,ALT)
      IMPLICIT NONE
      INTEGER(4)::Z,E,DCODE,DKIND,S
      INTEGER(4)::I,J
      REAL(8)::ALT,nFTD,XA
      REAL(8)::YOUT,RENORM
      CHARACTER*3::DIAGNOSE
      REAL(8),DIMENSION(27,100)::KEPN
      REAL(8),DIMENSION(19)::X,DY,DYY
      REAL(8),DIMENSION(19,19)::DY2
      REAL(8),DIMENSION(19)::DEPTH
      REAL(8),DIMENSION(19)::CEGRID
      REAL(8),DIMENSION(100)::EGRID
      REAL(8),DIMENSION(27,19,19,38)::FLUENCE,E103,E60,H10,GRAY
      REAL(8),DIMENSION(27,19,19,38)::FLUENCES,E103S,E60S,H10S,GRAYS
      INTEGER(4)::N
      Real(8)::A, m 
      DIMENSION::A(28), m(28), N(28)
!
      COMMON /atomic/N,A,m
      COMMON /energy/EGRID,CEGRID
      COMMON /kepn/KEPN
      COMMON /phys/DEPTH 
      COMMON /ftdcc/FLUENCE,E103,E60,H10,GRAY
      COMMON /ftdccs/FLUENCES,E103S,E60S,H10S,GRAYS

 !    DIAGNOSE='YES'
      DIAGNOSE='NO!'

      !get ftdccs at all altitudes from data file and interpolate      
      CALL DECODER(DCODE,DKIND,S)
      IF (DKIND.EQ.1) THEN 
         DO J=1,19
            DO I=1,19
               DY2(I,J)=FLUENCE(1,J,I,S)!recall 1st var is Z+1, protons=2
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.3) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E60(1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.4) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=H10(1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.5) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=GRAY(1,J,I,S)
            ENDDO
         ENDDO
      ELSE ! DEFAULT
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E103(1,J,I,S)
            ENDDO
         ENDDO
      ENDIF
      ! create a grid for this KEPN for all altitudes
      X=N(z)*CEGRID*0.001
      DO I=1,19
         DO J =1,19
            DYY(J)=DY2(I,J)
         ENDDO
         CALL KCSPLINE(19,X,DYY,1.0e+30,1.0e+30,KEPN(Z,E),YOUT)
         IF (DIAGNOSE.eq.'YES') THEN
            DO J=1,19
               WRITE(130,*) 'E(I)=',X(J), ' nFTD(I)=', DYY(J)
            ENDDO
            WRITE(130,*)'nFTD for KEPN',KEPN(Z,E),' is ',YOUT  
         ENDIF     
         DY(I)=YOUT 
      ENDDO
      X=DEPTH
      IF (ALT.LT.1000) then 
      ! interpolate
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,ALT,YOUT) 
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' nFTD(I)=', DY(I)
          ENDDO
          WRITE(130,*)'nFTD for ALT',ALT,' is ',YOUT  
        ENDIF
      ELSE
       XA=1000.
! KC 20160918 DD=132 is doubling depth estimate based on LANL report 
!    LA-UR-14-21331 (2014)
!Extrapolation using power law, assume doubling depth of DD g/cm2
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,XA,YOUT)
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' nFTD(I)=', DY(I)
          ENDDO
          WRITE(130,*)'nFTD for 1000 gpcs is ',YOUT  
        ENDIF
        YOUT=YOUT*0.5**((ALT-1000.)/132.)
        WRITE(130,*)'nFTD extrapolated to ALT ',ALT,' is ',YOUT  
      ENDIF   
      nFTD=YOUT*RENORM(DKIND)
      END FUNCTION nFTD 
!                                                                      7
!----6-----------------------------------------------------------------2
! renormalization factor for flux instead of dose calculations, 
! DOSE is normalized to microSv/hr in M_READER.FOR
! but flux in shower files has not yet been properly normalized       
      FUNCTION RENORM(I)
      IMPLICIT NONE
      REAL(8)::RENORM,SCALE
      INTEGER(4)::I
!
! KC 20161030 SCALE value is based on reported inbound ion fluences at 1 TeV 3.70E-19  
!             before other normalizations    
      SCALE=1.22 !geometric factor needed to make mcnpx reported particle fluence 
!                 equal expected fluence at the top of the atmosphere

      IF (I.EQ.1) THEN !particle flux must be renormalized 
         RENORM=5.262E+18/(SCALE*3600)! = Sky area in cm^2 /(3600 sec/hour)
         ! 1. sky area is not included in fluence tables from m_reader.for   
         ! 2. output spectra are in units of /s, but doses are per hour  
      ELSE !dose rates are handled by norm
         RENORM=1.0/SCALE  
      ENDIF       

      END FUNCTION 
!                                                                      7
!----6-----------------------------------------------------------------2
! finds shower+ftdcc based uncertainty in the ftdcc data at alt 
! for specified primary z, e, and dosekind
      FUNCTION SHOWS(Z,E,DCODE,ALT)
      IMPLICIT NONE
      INTEGER(4)::I,J,Z,E,DCODE,DKIND,S
      REAL(8)::ALT,SHOWS,XA
      REAL(8)::YOUT,XE,RENORM
      CHARACTER*3::DIAGNOSE
      REAL(8),DIMENSION(19)::X,Y,DY,DYY,DEPTH
      REAL(8),DIMENSION(19,19)::DY2
      REAL(8),DIMENSION(27,19)::KEPN      
      REAL(8),DIMENSION(19)::CEGRID
      REAL(8),DIMENSION(100)::EGRID
      REAL(8),DIMENSION(27,19,19,38)::FLUENCE,E103,E60,H10,GRAY
      REAL(8),DIMENSION(27,19,19,38)::FLUENCES,E103S,E60S,H10S,GRAYS
      INTEGER(4)::N
      Real(8)::A, m 
      DIMENSION::A(28), m(28), N(28)
!
      COMMON /atomic/N,A,m
      COMMON /energy/EGRID,CEGRID      
      COMMON /kepn/KEPN
      COMMON /phys/DEPTH 
      COMMON /ftdcc/FLUENCE,E103,E60,H10,GRAY
      COMMON /ftdccs/FLUENCES,E103S,E60S,H10S,GRAYS

!      DIAGNOSE='YES'
      DIAGNOSE='NO!'

      XE=EGRID(E)*N(Z)
      CALL DECODER(DCODE,DKIND,S)
      IF (DKIND.EQ.1) THEN 
         DO J=1,19
            DO I=1,19
               DY2(I,J)=FLUENCES(Z+1,J,I,S)!recall 1st var is Z+1, protons=2
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.3) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E60S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.4) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=H10S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.5) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=GRAYS(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSE ! DEFAULT
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E103S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ENDIF

      X=CEGRID*0.001

      DO I=1,19
         DO J =1,19
            DYY(J)=DY2(I,J)
!            IF (DIAGNOSE.eq.'YES') THEN
!               WRITE(130,*) 'DYY(I)=',DYY(J), 'DY2(I,J)=', DY2(I,J)
!            ENDIF     
         ENDDO
         CALL KCSPLINE(19,X,DYY,1.0e+30,1.0e+30,XE,YOUT)
         IF (DIAGNOSE.eq.'YES') THEN
            DO J=1,19
               WRITE(130,*) 'E(I)=',X(J), ' UNCERTAINTY(I)=', DYY(J)
            ENDDO
            WRITE(130,*)'UNCERTAINTY for E',XE,' at alt ',I,' is ',YOUT  
         ENDIF     
         DY(I)=YOUT 
      ENDDO
      X=DEPTH
      XA=1000.
      IF (ALT.LT.XA) then 
      ! interpolate
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,ALT,YOUT) 
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' UNCERTAINTY(I)=', DY(I)
          ENDDO
          WRITE(130,*)'UNCERTAINTY for ALT',ALT,' is ',YOUT  
        ENDIF
      ELSE
! KC 20160918 DD is doubling depth estimate based on LANL report 
!    LA-UR-14-21331 (2014)
!Extrapolation using power law, assume doubling depth of DD g/cm2
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,XA,YOUT)
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' UNCERTAINTY(I)=', DY(I)
          ENDDO
          WRITE(130,*)'UNCERTAINTY for ', XA ,' gpcs is ',YOUT  
        ENDIF
        YOUT=YOUT*0.5**((ALT-XA)/132.)
!        WRITE(130,*)'UNCERTAINTY extrapolated to ALT ',ALT,' is ',YOUT  
      ENDIF   
     
      SHOWS=YOUT*RENORM(DKIND)
      END FUNCTION SHOWS 
!                                                                      7
!----6-----------------------------------------------------------------2
! Finds proton shower+ftdcc based uncertainty in the ftdcc data at alt 
! for specified z, e, and dosekind
      FUNCTION pSHOWS(Z,E,DCODE,ALT)
      IMPLICIT NONE
      INTEGER(4)::I,J,Z,E,DCODE,DKIND,S
      REAL(8)::ALT,pSHOWS,XA
      REAL(8)::YOUT,RENORM
      CHARACTER*3::DIAGNOSE
      REAL(8),DIMENSION(19)::X,Y,DY,DYY,DEPTH
      REAL(8),DIMENSION(19,19)::DY2
      REAL(8),DIMENSION(27,19)::KEPN      
      REAL(8),DIMENSION(19)::CEGRID
      REAL(8),DIMENSION(100)::EGRID
      REAL(8),DIMENSION(27,19,19,38)::FLUENCE,E103,E60,H10,GRAY
      REAL(8),DIMENSION(27,19,19,38)::FLUENCES,E103S,E60S,H10S,GRAYS
      INTEGER(4)::N
      Real(8)::A, m 
      DIMENSION::A(28), m(28), N(28)
!
      COMMON /atomic/N,A,m
      COMMON /energy/EGRID,CEGRID      
      COMMON /kepn/KEPN
      COMMON /phys/DEPTH 
      COMMON /ftdcc/FLUENCE,E103,E60,H10,GRAY
      COMMON /ftdccs/FLUENCES,E103S,E60S,H10S,GRAYS

!      DIAGNOSE='YES'
      DIAGNOSE='NO!'

      CALL DECODER(DCODE,DKIND,S)
      IF (DKIND.EQ.1) THEN 
         DO J=1,19
            DO I=1,19
               DY2(I,J)=FLUENCES(Z+1,J,I,S)!recall 1st var is Z+1, protons=2
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.3) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E60S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.4) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=H10S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.5) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=GRAYS(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSE ! DEFAULT
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E103S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ENDIF
      ! create a grid for this KEPN for all altitudes
      X=CEGRID*0.001
      DO I=1,19
         DO J =1,19
            DYY(J)=DY2(I,J)
!            IF (DIAGNOSE.eq.'YES') THEN
!               WRITE(130,*) 'DYY(I)=',DYY(J), 'DY2(I,J)=', DY2(I,J)
!            ENDIF     
         ENDDO
         CALL KCSPLINE(19,X,DYY,1.0e+30,1.0e+30,KEPN(Z,E),YOUT)
         IF (DIAGNOSE.eq.'YES') THEN
            DO J=1,19
               WRITE(130,*) 'E(I)=',X(J), ' pSHOWS(I)=', DYY(J)
            ENDDO
            WRITE(130,*)'pSHOWS for KEPN',KEPN(Z,E),' is ',YOUT  
         ENDIF     
         DY(I)=YOUT 
      ENDDO
      X=DEPTH
      XA=1000.
      IF (ALT.LT.XA) then 
      ! interpolate
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,ALT,YOUT) 
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' PSHOWS(I)=', DY(I)
          ENDDO
          WRITE(130,*)'pSHOWS for ALT',ALT,' is ',YOUT  
        ENDIF
      ELSE
! KC 20160918 DD is doubling depth estimate based on LANL report 
!    LA-UR-14-21331 (2014)
!Extrapolation using power law, assume doubling depth of DD g/cm2
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,XA,YOUT)
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' pSHOWS(I)=', DY(I)
          ENDDO
          WRITE(130,*)'pSHOWS for ', XA,' gpcs is ',YOUT  
        ENDIF
        YOUT=YOUT*0.5**((ALT-XA)/132.)
!        WRITE(130,*)'PSHOWS extrapolated to ALT ',ALT,' is ',YOUT  
      ENDIF   
      pSHOWS=YOUT*RENORM(DKIND)
      END FUNCTION pSHOWS 
!                                                                      7
!----6-----------------------------------------------------------------2
! Finds neutron shower+ftdcc based uncertainty in the ftdcc data at alt 
! for specified z, e, and dosekind
      FUNCTION nSHOWS(Z,E,DCODE,ALT)
      IMPLICIT NONE
      INTEGER(4)::I,J,Z,E,DCODE,DKIND,S
      REAL(8)::ALT,nSHOWS,XA
      REAL(8)::YOUT,RENORM
      CHARACTER*3::DIAGNOSE
      REAL(8),DIMENSION(19)::X,Y,DY,DYY,DEPTH
      REAL(8),DIMENSION(19,19)::DY2
      REAL(8),DIMENSION(27,19)::KEPN      
      REAL(8),DIMENSION(19)::CEGRID
      REAL(8),DIMENSION(100)::EGRID
      REAL(8),DIMENSION(27,19,19,38)::FLUENCE,E103,E60,H10,GRAY
      REAL(8),DIMENSION(27,19,19,38)::FLUENCES,E103S,E60S,H10S,GRAYS
      INTEGER(4)::N
      Real(8)::A, m 
      DIMENSION::A(28), m(28), N(28)
!
      COMMON /atomic/N,A,m
      COMMON /energy/EGRID,CEGRID      
      COMMON /kepn/KEPN
      COMMON /phys/DEPTH 
      COMMON /ftdcc/FLUENCE,E103,E60,H10,GRAY
      COMMON /ftdccs/FLUENCES,E103S,E60S,H10S,GRAYS

!      DIAGNOSE='YES'
      DIAGNOSE='NO!'

      CALL DECODER(DCODE,DKIND,S)
      IF (DKIND.EQ.1) THEN 
         DO J=1,19
            DO I=1,19
               DY2(I,J)=FLUENCES(Z+1,J,I,S)!recall 1st var is Z+1, protons=2
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.3) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E60S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.4) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=H10S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSEIF (DKIND.EQ.5) THEN
         DO J=1,19
            DO I=1,19
               DY2(I,J)=GRAYS(Z+1,J,I,S)
            ENDDO
         ENDDO
      ELSE ! DEFAULT
         DO J=1,19
            DO I=1,19
               DY2(I,J)=E103S(Z+1,J,I,S)
            ENDDO
         ENDDO
      ENDIF

      ! create a grid for this KEPN for all altitudes
      X=N(z)*CEGRID*0.001
      DO I=1,19
         DO J =1,19
            DYY(J)=DY2(I,J)
         ENDDO
         CALL KCSPLINE(19,X,DYY,1.0e+30,1.0e+30,KEPN(Z,E),YOUT)
         IF (DIAGNOSE.eq.'YES') THEN
            DO J=1,19
               WRITE(130,*) 'E(I)=',X(J), ' nSHOWS(I)=', DYY(J)
            ENDDO
            WRITE(130,*)'nFTD for KEPN',KEPN(Z,E),' is ',YOUT  
         ENDIF     
         DY(I)=YOUT 
      ENDDO
      X=DEPTH
      XA=1000.
      IF (ALT.LT.XA) then 
      ! interpolate
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,ALT,YOUT) 
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' nSHOWS(I)=', DY(I)
          ENDDO
          WRITE(130,*)'nSHOWS for ALT',ALT,' is ',YOUT  
        ENDIF
      ELSE
! KC 20160918 DD is doubling depth estimate based on LANL report 
!    LA-UR-14-21331 (2014)
!Extrapolation using power law, assume doubling depth of DD g/cm2
        CALL KCSPLINE(19,X,DY,1.0e+30,1.0e+30,XA,YOUT)
        IF (DIAGNOSE.eq.'YES') THEN
          DO I=1,19
            WRITE(130,*) 'ALT(I)=',DEPTH(I), ' nSHOWS(I)=', DY(I)
          ENDDO
          WRITE(130,*)'nSHOWS for ',XA ,' gpcs is ',YOUT  
        ENDIF
        YOUT=YOUT*0.5**((ALT-XA)/132.)
!        WRITE(130,*)'nSHOWS extrapolated to ALT ',ALT,' is ',YOUT  
      ENDIF   
      nSHOWS=YOUT*RENORM(DKIND)
      END FUNCTION nSHOWS 
!                                                                      7
!----6-----------------------------------------------------------------2
      SUBROUTINE DECODER(DCODE,DKIND,S)
      IMPLICIT NONE
      INTEGER(4)::DCODE,DKIND,S

!     S = SECONDARY PARTICLE INDICATOR IN *SHOWER.DAT, 
!         1 = neutrons
!         2 = photons
!         3 = e-
!          ...
!        38 = Co !!!INACTIVE
!        39 = Ni !!!INACTIVE
!        40 = total
!     DKIND = SECONDARY OUTPUT, 1-5: flux, E103, E60, H10, GRAY 
!
!     Since 38 and 39 are inactive numbers, total is at 38
!     For dkind=1 S=38 is ions/cm2
!      
!      IF (MOD(DCODE,40).EQ.0) THEN
!         DCODE=38*DCODE/40
!      ENDIF
!     Assign output kind and secondary designators
      IF (DCODE.LE.38) THEN
         DKIND=1
         S=DCODE
      ELSEIF((DCODE.GT.38).AND.(DCODE.LE.76)) THEN
         DKIND=2
         S=DCODE-38
      ELSEIF((DCODE.GT.76).AND.(DCODE.LE.114)) THEN
         DKIND=3
         S=DCODE-76
      ELSEIF((DCODE.GT.114).AND.(DCODE.LE.152)) THEN
         DKIND=4
         S=DCODE-114
      ELSEIF((DCODE.GT.152).AND.(DCODE.LE.190)) THEN
         DKIND=5
         S=DCODE-152
      ENDIF
      END SUBROUTINE 
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION LE(Z,E)
! Lower energy of nucleus for dose rate calculations
      IMPLICIT NONE
      INTEGER(4)::Z,E
      REAL(8),DIMENSION(19)::CEGRID
      REAL(8),DIMENSION(100)::EGRID
      INTEGER(4)::N
      Real(8)::A, m, LE, UE 
      DIMENSION::A(28), m(28), N(28)
      CHARACTER*3::DIAGNOSE
!
      COMMON /atomic/N,A,m
      COMMON /energy/EGRID,CEGRID
!
      DIAGNOSE='No!'

      IF (E.EQ.1) THEN
         LE=0.
         IF (DIAGNOSE.eq.'YES') WRITE(130,*) 'lower E is', LE
         RETURN
      ELSE
         LE=UE(Z,(E-1))
         IF (DIAGNOSE.eq.'YES') WRITE(130,*) 'lower E is', LE
         RETURN
      ENDIF
      END FUNCTION LE
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION UE(Z,E)
! Upper energy of nucleus for dose rate calculations
      INTEGER(4)::Z,E
      REAL(8),DIMENSION(19)::CEGRID
      REAL(8),DIMENSION(100)::EGRID
      INTEGER(4)::N
      Real(8)::A, m, UE 
      DIMENSION::A(28), m(28), N(28)
      CHARACTER*3::DIAGNOSE
!
      COMMON /atomic/N,A,m
      COMMON /energy/EGRID,CEGRID

      DIAGNOSE='No!'
!
      IF (E.NE.100) THEN  
         UE=N(Z)*(EGRID(E)+((EGRID(E+1)-EGRID(E))/2))
         IF (DIAGNOSE.eq.'YES') WRITE(130,*) 'Upper E is', UE
      ELSE 
         UE=N(Z)*(EGRID(E)+((EGRID(E)-EGRID(E-1))/2)) 
         IF (DIAGNOSE.eq.'YES') WRITE(130,*) 'Upper E is', UE
      ENDIF
      END FUNCTION UE
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION FHF(gcrflux,gcrKE,Z,E) 
!Puts all remaining estimated integral flux in energy flux bin 100
      IMPLICIT NONE
      INTEGER(4)::Z,E
      REAL(8)::alpha,f0,x1,x2,E1,E2
      REAL(8),DIMENSION(27,100)::gcrflux,gcrKE
      REAL(8)::ANSWER,FHF
!
      IF (E.eq.100) then
         x1=gcrflux(Z,100)
         x2=gcrflux(Z,99)
         E1=gcrKE(Z,100)
         E2=gcrKE(Z,99)
         alpha = -1*log(x2/x1)/log(E2/E1)  
         f0=x1/(E1**-alpha)
         ANSWER=f0*(alpha-1)*E1**(1-alpha)/gcrflux(Z,E)
      ELSE
         ANSWER=1.0
      ENDIF 
      FHF=ANSWER
      END FUNCTION FHF
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION km(G)
            IMPLICIT NONE
!     convert gpcms to km for use by SKY_LIBS
      REAL(8)::gpcms,km
         INTEGER(4) :: I,J
         REAL(8):: F,G
         REAL(8), DIMENSION(553) :: KFT, GRAMS
         CHARACTER(10)::HEADER1,HEADER2
         CHARACTER(3)::DIAGNOSE='no ' !'YES'

         OPEN(UNIT=29,FILE='ATMOSPHERE/FT-GM.DAT',STATUS='OLD')          
         READ(29,9100) HEADER1, HEADER2
         DO I = 1,553
            READ(29,*) KFT(I), GRAMS(I)
            KFT(I)=KFT(I)
!      KFT(I) TABLE IS IN 1000'S OF FEET
            IF (DIAGNOSE.EQ.'YES') WRITE(40,*) KFT(I), GRAMS(I)
         ENDDO
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'grams IN:', G
         IF ((G .LE. 1035.0837).AND.(G.GE.0.0001)) THEN
            CALL KCSPLINE(553,GRAMS, KFT,1.0E+30,1.0E+30,G,F)
         ELSE !ALTITUDE IS OUT OF RANGE
            CLOSE(29)
            CALL EPITATH(' ALTITUDE CONVERSION OUT OF RANGE ',35)
         ENDIF
9100  FORMAT(A10,A10)
      CLOSE(29)
      IF (DIAGNOSE.EQ.'YES') THEN
         WRITE(40,*)'G IN: ', G, ' kft OUT:', F,' km out: ',F*0.3048
      ENDIF 
      km=F*0.3048       

      END FUNCTION km
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION fluxmod(z,e)
!     function for playing with flux inputs
      IMPLICIT NONE
         INTEGER(4) :: z,e
         REAL(8)::fluxmod,pi

         pi=3.1415926 

         IF (z.eq.1) then
            IF ((e.GT.100).OR.(e.LT.1)) then
               fluxmod = 0.00000001
            else
               fluxmod = 1.
            endif
         ELSE
            IF ((e.GT.100).OR.(e.LT.1)) then
               fluxmod = 0.00000001
            else
               fluxmod = 1.
            endif
         ENDIF
      END FUNCTION fluxmod
!                                                                      7
!----6-----------------------------------------------------------------2
       
