!      Code by Kyle Copeland from the ISO GCR source standard TS15390:2004
!      This is the same as the MSU and Nymmik GCR spectra, but uses
!      heliocentric potential for modulation
!
!      Data are from unit 8
!      Diagnostic print is to unit 1
!      Output is to unit 2
!      Input is from Keyboard for testing, CARI-7 for use
!
!-----------------------------------------------------------------------&
!     Rigidity from Kinetic Energy Per Nucleon
!     Checked 12/10/2014, for discussion see Smart and Shea, Chapt. 6, 
!                         Handbook of Geophysics, 1985
! 
      function Rigidity(T,Z)
      REAL(8)::R,T,Rigidity
      INTEGER(4)::I,J,N,Z
      REAL(8)::A, m 
      DIMENSION::A(28), m(28)

      Data m/0.938,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939&
     &,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939&
     &,0.939,0.939,0.939,0.939,0.939,0.939,0.939/
      Data A/1.,4.,6.9,9.,10.8,12.,14.,16.,19.,20.2,23.,24.3,27.,28.1,  &
     & 31.,32.1,35.4,39.9,39.1,40.1,44.9,47.9,50.9,52.,54.9,55.8,58.9,  &
     & 58.7/

	R=(A(Z)/ABS(Z))*SQRT(T*(T+2*m(Z)))
      Rigidity=R         
      end function 
!-----------------------------------------------------------------------&
! Below uses heliocentric potential to modulate ISO instead of Wolf numbers
       SUBROUTINE GET_HPD_ISO(t) !modified from ISO15390
!
       IMPLICIT NONE
!
       REAL(8)::R,t,phi,sigmaphi,F,sigmaF,varF,HP,E,K
       REAL(8)::Rigidity,FINDU,U
       REAL(8), DIMENSION(100)::KEPN
       INTEGER(4)::Z,J 

       CHARACTER(10)::INIVAR
       CHARACTER(12)::INIVAL
       CHARACTER(12)::VIEWER 
       CHARACTER(5)::OS
       CHARACTER(4)::OUTPUT
       CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
       INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      
!      
!      R = particle rigidity in GV
!      E = particle energy/nucleus in GeV
!      phi = nucleons/(GV-m2-sr-s)    
!      F = nucleons/(MeV-m2-sr-s)    
!      Z = nuclear charge number
!      sigma*=standard deviation of *
!      t = date and time in yyyy.xxxx format
!      KEPN = particle energy in GeV/nucleon
!      U and HP are dummy variables
!
!      KEPN values are energies matching the primary 
!      100 pt flux energy grid of the BO11 spectrum. 
       DATA (KEPN(J),J=1,100)/
     &1.000E-02,1.150E-02,1.322E-02,1.520E-02,1.748E-02,2.009E-02,      &
     &2.310E-02,2.656E-02,3.054E-02,3.511E-02,4.037E-02,4.642E-02,      &
     &5.337E-02,6.136E-02,7.055E-02,8.111E-02,9.326E-02,1.072E-01,      &
     &1.233E-01,1.417E-01,1.630E-01,1.874E-01,2.154E-01,2.477E-01,      &
     &2.848E-01,3.275E-01,3.765E-01,4.329E-01,4.977E-01,5.722E-01,      &
     &6.579E-01,7.565E-01,8.697E-01,1.,1.15,1.322,1.520,1.748,2.009,    &
     &2.310,2.656,3.054,3.511,4.037,4.642,5.337,6.136,7.055,8.111,      &
     &9.326,1.072E+01,1.233E+01,1.417E+01,1.630E+01,1.874E+01,2.154E+01,&
     &2.477E+01,2.848E+01,3.275E+01,3.765E+01,4.329E+01,4.977E+01,      &
     &5.722E+01,6.579E+01,7.565E+01,8.697E+01,1.000E+02,                &
     &1.150E+02,1.322E+02,1.520E+02,1.748E+02,2.009E+02,2.310E+02,      &
     &2.656E+02,3.054E+02,3.511E+02,403.7,464.2,533.7,613.6,705.5,811.1,&
     &932.6,1.072E+03,1.233E+03,1.417E+03,1.630E+03,1.874E+03,2.154E+03,&
     &2.477E+03,2.848E+03,3.275E+03,3.765E+03,4.329E+03,4.977E+03,      &
     &5.722E+03,6.579E+03,7.565E+03,8.697E+03,1.000E+04/      

      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
! WINDOWS/DOS PATHS     
       OPEN(2,file='GCR_MODELS\ISO_GCR2.OUT',status='unknown')
      ELSE !LINUX PATHS
       OPEN(2,file='GCR_MODELS/ISO_GCR2.OUT',status='unknown')
      ENDIF

       WRITE(2,80100) t 
       WRITE(2,*)'  Z       E         F' 

       U=(FINDU(t)+150)/1000. !+250 MV corrects to an almost perfect match with 
!       ICRU Rep 84 data at FL310-390, but makes the ER-2 comparison less satisfactory.       

       DO Z=1,28
          DO J=1,100
             K=KEPN(J)
             R=Rigidity(K,Z)
             Call ISO_GCR_FLUX2(Z,R,K,E,U,phi,F)
             varF=sigmaF**2
             WRITE(2,80102) Z,E,F 
          end do        
!             WRITE(*,*) Z,R,E,F
       end do
! end writing of output grid
80100  FORMAT(F11.6)               
80102  FORMAT(I4,2ES11.3)
       CLOSE(2)
       CLOSE(3)
       CLOSE(8)
       END SUBROUTINE GET_HPD_ISO
!-----------------------------------------------------------------------&
       SUBROUTINE ISO_GCR_FLUX2(Z,R,KEPN,KE,U,phi,F)
!
!      ISO LIS modulated by the heliocentric potential method used in LUIN
!      as explained in many papers going back to the 60's, e.g., some of
!      Keran O'Brien and Wally Friedberg's joint Environment International 
!      papers from the 90's   
!
!      Variables
!       U = heliocentric potential in GV
!       (Adding 250 MV makes H*(10) match the ICRU Rep 84 data at high
!        latitude almost perfectly, HOWEVER, equatorial results are 
!        about 20% high relative to ICRU for any reasonable value.) 
!       Z = nuclear charge #
!       R = Rigidity in GV at 1 AU, i.e. after modulation   
!       t = time
!       KEPN = kinetic energy per nucleon, GeV/n
!       KE = kinetic energy per nucleus, GeV
!       R0 = Rigidity before heliocentric modulation, dependent on U(t)
!       E0 = KEPN before heliocentric modulation, dependent on U(t)
!       m(Z) = rest mass per nucleon 
!       A(Z) = atomic mass #
!       Am(Z) = rest mass energy of particle
!       alpha(Z) = parameter, spectral index 
!       beta = velocity relative to lightspeed
!       gamma(Z) = parameter, spectral index
!       C = parameter dependent on Z, related to relative abundance
!
      IMPLICIT NONE
      REAL(8)::F,phi,R,t,KE,KEPN,beta
      REAL(8)::E0,R0,U
      REAL(8)::Rnaught,Rigidity
      INTEGER(4)::I,J,N,Z
      REAL(8)::A, m, Am, C, sigma, sigmaC, gamma, alpha(28)
      DIMENSION::A(28), m(28), Am(28), C(28), sigmaC(28), gamma(28)

      Data m/0.938,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939&
     &,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939,0.939&
     &,0.939,0.939,0.939,0.939,0.939,0.939,0.939/
      Data A/1.,4.,6.9,9.,10.8,12.,14.,16.,19.,20.2,23.,24.3,27.,28.1,  &
     & 31.,32.1,35.4,39.9,39.1,40.1,44.9,47.9,50.9,52.,54.9,55.8,58.9,  &
     & 58.7/
      Data C/18500.,3690.,19.5,17.7,49.2,103.0,36.7,87.4,3.19,16.4,4.43,&
     & 19.3,4.17,13.4,1.15,3.06,1.30,2.33,1.87,2.17,0.74,2.63,1.23,2.12,&
     & 1.14,9.32,0.10,0.49/
! Here in case I want to add sigmas for this spectrum later
      Data sigmaC/1300,220,1.5,1.3,1.6,3.0,1.2,2.1,0.28,0.70,0.28,0.70, &
     & 0.22,0.50,0.04,0.12,0.08,0.07,0.05,0.06,0.02,0.08,0.04,0.06,0.05,&
     & 0.24,0.08,0.02/	  
      Data gamma/2.74,2.77,2.82,3.05,2.96,2.76,2.89,2.70,2.82,2.76,2.84,&
     & 2.70,2.77,2.66,2.89,2.71,3.00,2.93,3.05,2.77,2.97,2.99,2.94,2.89,&
     & 2.84,2.63,2.63,2.63/
      Data alpha/2.85,3.12,3.41,4.30,3.93,3.18,3.77,3.11,4.05,3.11,3.14,&
     & 3.65,3.46,3.00,4.04,3.30,4.40,4.33,4.49,2.93,3.78,3.79,3.50,3.28,&
     & 3.29,3.01,4.25,3.52/

        Am(Z)=A(Z)*m(Z)
        E0=KEPN+((U)*Z)/INT(A(Z)+0.5)
        R0=Rigidity(E0,Z)
        beta=R0/SQRT(R0**2+(Am(Z)/Z)**2) 
        phi=((C(Z)*beta**alpha(Z))/R0**gamma(Z))*(R/R0)**2
        F=phi*(A(Z)/Z)*0.001/beta
!
!     Renormalize to expected units 
!        KE = T of nucleus and 
!        flux = nuclei/(GeV-m2-s-sr)
!      
!    KE per nucleus  
        KE=KEPN*INT(A(Z)+.5)
!    revised flux
        F=F*1000./INT(A(Z)+.5)

      RETURN
      END SUBROUTINE ISO_GCR_FLUX2