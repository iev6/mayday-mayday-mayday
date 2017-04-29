!
!   SKY_LIBS.FOR   !21000-21999 ???? 
!   Looks up vertical cutoffs from Shea and Smart grids, interpolates 
!   to find vertical cutoff for a particular date and time, then uses a 
!   library of routines to calculate quiet nonvertical magnetic rigidity 
!   cutoffs based on S. ROESLER'S approach. 
!
!   Uses the Geomagnetic storm adjustment method from Hani Anid's Dissertation. 
!
!   Coding by Kyle Copeland for CARI-7 and Dissertation, begun 29 July 2013
!   Mods since dissertation
! KC 2015 switch from 10X15 degree sky slices to floating size, this uses the same 
!         number of slices but creates a ~2x finer integration oover solid angle. Now 
!         all the slices are at least patially unblocked by the Earth. 
!         This increased dose rates relative to the thesis.  
! KC 20161001 added variable MinCut to sub FINDCUT. this is for experimenting with 
!         low-R GCR influence on High-alt dose rates and difference between my calcs 
!         and those of PHITS/EXPACS using ISO spectrum. At highest alt I am about 2x EXPACS
!         in effective dose rate despite almost identical FTDCCs.
! KC 20161030 EXPACS differences resolved using geometric factor
! KC 20170404 Fixed a bug in FT2GPCMS converting alts below sea level 
! 
!----------------------------------------------------------------------
!
      subroutine FINDCUT(NA,NZ,t,lat,lon,altkm,VC,SkyCut,SkyPass,SkyWt)
      IMPLICIT NONE
      INTEGER:: ue,le,passes,uyr,lyr,MD,i,ii,NA,NZ 
      REAL(8):: DATE1,dspan,MYDATE
      REAL(8)::t,lat,lon,altkm,NVC,VC,VC1,VC2,NVFACTOR
      REAL(8)::w1,w2,alat,along,total,zstep,MinCut
! KC 20161001 added mincut to experiment with low-R GCR influence on High-alt dose rates
      REAL(8)::READVC,VCALTADJ,ANID,linterp,Earth,SKY
      REAL(8), DIMENSION(NZ,NA)::SkyCut,SkyPass,SkyWt
      
      CHARACTER(10)::INIVAR
      CHARACTER(12)::INIVAL
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE,CARIDIAG
      INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,CARIDIAG,VIEWER,OUTPUT
      COMMON /INI/SP,GCR
      COMMON /SKYDIAG/DIAGNOSE

      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
      !WINDOS/DOS PATHS
        OPEN(201,file='DIAGNOSE\DIAGNOSE.SKY',STATUS='UNKNOWN')
        OPEN(11,file='CUTOFFS\WGRC1965.1X1',STATUS='OLD')
        OPEN(12,file='CUTOFFS\IGRF1980.1X1',STATUS='OLD')
        OPEN(13,file='CUTOFFS\DGRF1990.1X1',STATUS='OLD')
        OPEN(14,file='CUTOFFS\IGRF1995.1X1',STATUS='OLD')
        OPEN(15,file='CUTOFFS\IGRF2000.1X1',STATUS='OLD')
        OPEN(16,file='CUTOFFS\IGRF2010.1X1',STATUS='OLD')

      ELSE !UNIX PATHS

        OPEN(201,file='DIAGNOSE/DIAGNOSE.SKY',STATUS='UNKNOWN')
        OPEN(11,file='CUTOFFS/WGRC1965.1X1',STATUS='OLD')
        OPEN(12,file='CUTOFFS/IGRF1980.1X1',STATUS='OLD')
        OPEN(13,file='CUTOFFS/DGRF1990.1X1',STATUS='OLD')
        OPEN(14,file='CUTOFFS/IGRF1995.1X1',STATUS='OLD')
        OPEN(15,file='CUTOFFS/IGRF2000.1X1',STATUS='OLD')
        OPEN(16,file='CUTOFFS/IGRF2010.1X1',STATUS='OLD')

      ENDIF
      MinCut = 0.0 !Low rigidity cutoff for all incident particles KC
      ! Set to 0 for normal use, .4346 GV roughly mimics LUIN 100 MeV 

      DIAGNOSE='NO!'
             
      DATE1 = t
      dspan = 30./365.
      MYDATE = DATE1+dspan/2. !date of center of query period
      MD = INT(MYDATE+.5) 
      SELECT CASE(MD)
      CASE (2000:2009)
         ue = 16
         uyr= 2010
         le = 15
         lyr= 2000
         passes=2
      CASE (1995:1999)
         ue = 15
         uyr= 2000
         le = 14
         lyr= 1995
         passes=2
      CASE (1990:1994)
         ue = 14
         uyr= 1995
         le = 13
         lyr= 1990
         passes=2
      CASE (1980:1989)
         ue = 13
         uyr= 1990
         le = 12
         lyr= 1980
         passes=2
      CASE (1965:1979)
         ue = 12
         uyr= 1980
         le = 11
         lyr= 1965
         passes=2
      CASE (:1964)
         uyr= 1965
         ue = 11
         lyr= 1965
         le = 11
         passes=1
      CASE DEFAULT !2010 and beyond)  
         ue = 16
         uyr= 2010
         le = 16
         lyr= 2010
         passes=1
      END SELECT
      IF (DIAGNOSE.EQ.'YES') THEN 
         Write(201,*) 'Based on date,',t, 'aka', MD
         Write(201,*) 'There should be ',passes,' passes'
         Write(201,*) 'Upper year ',uyr,' lower year ', lyr
         Write(201,*) 'Longitude sought is', lon !Diagnostic
         Write(201,*) 'Latitude sought is', lat !Diagnostic
      ENDIF
! get vertical cutoffs
      VC1=READVC(ue,lat,lon)*VCALTADJ(ue,altkm)
      if (passes.eq.1) THEN
         VC=VC1
      else    
         VC2=READVC(le,lat,lon)*VCALTADJ(le,altkm)
         IF (DIAGNOSE.EQ.'YES') THEN
            WRITE(*,'(2I6,2ES11.3,F13.7)')uyr,lyr,VC1,VC2,t
         endif
         VC=linterp(REAL(uyr,8),REAL(lyr,8),VC1,VC2,t)
      endif
      IF (DIAGNOSE.EQ.'YES') THEN 
         WRITE(*,*)'Pre-ANID Vertical cutoff is ',VC !Diagnostic
         Write(201,*)'Pre-ANID Vertical cutoff is ',VC !Diagnostic
      ENDIF
! Correct for geomagnetic disturbances
      VC = ANID(VC,t)
      IF (DIAGNOSE.EQ.'YES') THEN  
         Write(*,*) 'ANID Vertical cutoff is ',VC !Diagnostic
         Write(201,*) 'ANID Vertical cutoff is ',VC !Diagnostic
      ENDIF
! generate nonvertical cutoffs
      total=0.0 !checksum
      DO i = 1,NZ
         alat=90.-(2*i-1)*zstep(NZ,altkm)
         DO ii=1,NA
            along=(360./NA)*ii
            SkyCut(i,ii)=NVC(VC,NZ,NA,alat,along)
            IF (SkyCut(i,ii).LT.MinCut) THEN
               SkyCut(i,ii)=MinCut
            ENDIF
            SkyPass(i,ii)=Earth(NZ,NA,alat,altkm) 
            SkyWt(i,ii)=SKY(NZ,NA,alat,altkm)*NVFACTOR(VC)
            total=total+SKY(NZ,NA,alat,altkm)
            IF (DIAGNOSE.EQ.'YES') THEN 
              Write(201,21122) alat, along, SkyCut(i,ii), SkyPass(i,ii),&!diagnostic
     &                          SkyWt(i,ii)
           ENDIF 
         enddo
      enddo
      IF (DIAGNOSE.EQ.'YES') Write(201,*)'total sky weight is=', total
21122 Format(5(2x,F10.5))             
      end subroutine FINDCUT
!----------------------------------------------------------------------
! This model is from Hani Al-Anid's RMC dissertation on SPEs
      function ANID(cut,dyear)
      IMPLICIT NONE
      INTEGER(4)::Kp, GETKP
      REAL(8)::VC,cut,ANID,dyear
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE

!     Kp=0 is standard condition, no ongoing geomagnetic storms
!	  Kp=0

        IF (DIAGNOSE.EQ.'YES') THEN 
           WRITE(201,*) 'Searching for Kp index for date', dyear
        ENDIF
        Kp = GETKP(dyear)       
        if (Kp.GE.5) then 
           VC = (cut+1.-0.54*EXP(-cut/2.9))/2.
        else
           VC = cut
        endif
        ANID=VC		
      end function ANID
!----------------------------------------------------------------------
! Altitude correction to vertical cutoff
      Function VCALTADJ(epoch,altkm)
      IMPLICIT NONE
      INTEGER(4)::epoch
      REAL(8)::RTAB,altkm,h,R0,x,VCALTADJ
      REAL(8), Dimension(6)::Tablekm
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE
!
      DATA Tablekm/20.,20.,20.,20.,20.,20./      
!     R0 = Radius of Earth, in km
!     RTAB = Radius of Earth + table alt, in km
!     alt = altitude of plane, in km
!     h = Radius of Earth + alt, in km
      if(DIAGNOSE.EQ.'YES')Write(201,*)'Calculating altitude adjusment.'
      R0=6371.
      RTAB=R0+tablekm(epoch-10)  
      h=R0+altkm  ! plane alt in km
      x=RTAB**2/h**2
      VCALTADJ=x
      if (DIAGNOSE.EQ.'YES') Write(201,*)'Altitude adjustment:',1/x !Diagnostic
      end function VCALTADJ
!----------------------------------------------------------------------
!     READVC reads Shea and Smart table of four surrounding values
!     and interpolates to find VC 
      function READVC(u,la,lon)
      IMPLICIT NONE
      REAL(8)::la,lon,dla,dlo
      REAL(8)::VMCR,ala,VC11,VC12,VC21,VC22,VC,READVC,fx,fy
      REAL(8), Dimension(189,360)::GRIDVC   
      INTEGER(4)::GEOLAT1,GEOLON1,GEOLAT2,GEOLON2
      INTEGER(4)::u,S
      CHARACTER(128)::headline
      CHARACTER(2)::BC
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE
! find surrounding coordinates
!  Range in tables is 89N to 89S, 0E to 359E 
      if (DIAGNOSE.EQ.'YES') then
         Write(201,*) 'Finding bracketing locations in table' !Diagnostic
         Write(201,*) 'Incoming Latitude ', la !Diagnostic
         Write(201,*) 'Incoming Longitude ', lon !Diagnostic
      endif
      dla=la
      dlo=lon
      if (la.GE.0.) then
         S=1
      else
         S=-1
      endif 
      ala=ABS(dla)      
      GEOLAT1 = INT(ala,kind=4)
      if (GEOLAT1.LT.89) then !bracketing values are in the table  
         GEOLAT2 = INT(ala+1,kind=4)
      else 
         GEOLAT1 = 89
         GEOLAT2 = GEOLAT1
      endif
         GEOLAT1 = S*GEOLAT1
         GEOLAT2 = S*GEOLAT2
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'Bracketing Latitudes are ',  &
     &    GEOLAT1, GEOLAT2!Diagnostic
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'Incoming Longitude ', lon !Diagnostic
      DO while (dlo.lt.0.)  
         dlo=dlo+360.
      enddo
      DO while (dlo.ge.360.)  
         dlo=dlo-360.
      enddo
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'In-range Longitude ', dlo !Diagnostic
      GEOLON1=INT(dlo,kind=4)       
      if (GEOLON1.eq.359) then
         GEOLON2=0
      else
         GEOLON2=GEOLON1+1
      endif  
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'Bracketing Longitudes are ', &
     & GEOLON1,GEOLON2!Diagnostic
! get VMCR values for lat,lon pairs in table
      if (DIAGNOSE.EQ.'YES') Write(201,*)'Finding nearest table cutoffs' !Diagnostic
      VC11=VMCR(u,GEOLAT1,GEOLON1)    
      if (DIAGNOSE.EQ.'YES') Write(201,*)'In-range Longitude ', dlo !Diagnostic
      if (DIAGNOSE.EQ.'YES') Write(201,*)'Bracketing Longitudes are ',  &
     &   GEOLON1,GEOLON2!Diagnostic
      VC12=VMCR(u,GEOLAT1,GEOLON2)    
      if (DIAGNOSE.EQ.'YES') Write(201,*)'In-range Longitude ', dlo !Diagnostic
      GEOLON1=INT(dlo,kind=4) !refresh longitudes       
      if (GEOLON1.eq.359) then
         GEOLON2=0
      else
         GEOLON2=GEOLON1+1
      endif
      if (DIAGNOSE.EQ.'YES') Write(201,*)'Bracketing Longitudes are ',  &
     & GEOLON1,GEOLON2!Diagnostic
      VC21=VMCR(u,GEOLAT2,GEOLON1)    
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'In-range Longitude ', dlo !Diagnostic
      GEOLON1=INT(dlo,kind=4) !refresh longitudes      
      if (GEOLON1.eq.359) then
         GEOLON2=0
      else
         GEOLON2=GEOLON1+1
      endif
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'Bracketing Longitudes are ', &
     & GEOLON1,GEOLON2!Diagnostic
      VC22=VMCR(u,GEOLAT2,GEOLON2)
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'In-range Longitude ', dlo !Diagnostic
      GEOLON1=INT(dlo,kind=4) !refresh longitudes      
      if (GEOLON1.eq.359) then
         GEOLON2=0
      else
         GEOLON2=GEOLON1+1
      endif 
      if (DIAGNOSE.EQ.'YES') then 
         Write(201,*) 'Bracketing Longitudes are thus', GEOLON1,GEOLON2!Diagnostic
         Write(201,21206)'Box cutoffs are:', VC11,VC12,VC21,VC22 !Diagnostic
      endif
      fx = dlo-GEOLON1
      fy = ala-ABS(GEOLAT1)
      CALL BOXINTERP(VC11,VC12,VC21,VC22,fx,fy,VC)
      if (DIAGNOSE.EQ.'YES') Write(201,*)'VCut after boxinterp', VC
      READVC=VC
21206 Format(a16,4f7.3)           
      end function READVC
!----------------------------------------------------------------------
! Linear Interpolation inside a 1x1 box
      SUBROUTINE BOXINTERP(A,B,C,D,fx,fy,VALatP)
      REAL(8)::A,B,C,D,VALatP,fx,fy
      REAL(8)::U,L,R
      
      U=(1.-fx)*C+fx*D
      L=(1.-fx)*A+fx*B
      R=(1.-fy)*L+fy*U
      VALatP=R
                
      end SUBROUTINE BOXINTERP
!----------------------------------------------------------------------
! Get data from Shea and Smart Table (Vertical Magnetic Cutoff Rigidity) 
      FUNCTION VMCR(fnum,lat,lon)

      INTEGER(4)::fnum,lat,lon,i,skip
      INTEGER(4)::VMCRPAGE,VMCRCOL,VMCRLINE
      REAL(8)::cut
      CHARACTER(128)::filler
      CHARACTER(2)::bs
      INTEGER(4)::RLAT
      REAL(8),DIMENSION(20)::tval
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE

      REWIND(fnum)
      skip=((VMCRPAGE(lon)-1)*216+VMCRLINE(lat)-1)
      DO i=1,skip
        Read(fnum,'(A128)') filler
      ENDDO
      Read(fnum,21001) RLAT,bs,tval(1),tval(2),tval(3),tval(4),tval(5), &
     & tval(6),tval(7),tval(8),tval(9),tval(10),tval(11),tval(12),      &
     & tval(13),tval(14),tval(15),tval(16),tval(17),tval(18),tval(19),  &
     & tval(20)
      if (DIAGNOSE.EQ.'YES') then 
         Write(201,*) 'Called VMCR', lat, lon, RLAT
         DO i =1,20
            Write(201,*) tval(i)	  
         ENDDO
      endif
21001 format(I5,A2,20F6.2)
      cut = tval(VMCRCOL(lon))
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'VMCR =', cut
      VMCR=cut
      END FUNCTION VMCR
!----------------------------------------------------------------------
! Finds page containing the Vertical Magnetic Cutoff Rigidity from longitude
      function VMCRPAGE(L)
      IMPLICIT NONE
      INTEGER(4)::L,page,VMCRPAGE
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE

      page = 1 + L / 20       
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'page=',page
      VMCRPAGE=page

      end function VMCRPAGE
!----------------------------------------------------------------------
! Finds line on page containing the VMCR datum from latitude
      function VMCRLINE(L)
      IMPLICIT NONE
      INTEGER(4)::L,VMCRLINE,row
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE
      IF (L.GE.80) then
      row = 5+90-L
      elseif((L.GE.70).AND.(L.LT.80)) then
      row = 16+(80-L)
      elseif((L.GE.60).AND.(L.LT.70)) then
      row = 27+(70-L)
      elseif((L.GE.50).AND.(L.LT.60)) then
      row = 38+(60-L)
      elseif((L.GE.40).AND.(L.LT.50)) then
      row = 49+(50-L)
      elseif((L.GE.30).AND.(L.LT.40)) then
      row = 65+(40-L)
      elseif((L.GE.20).AND.(L.LT.30)) then
      row = 76+(30-L)
      elseif((L.GE.10).AND.(L.LT.20)) then
      row = 87+(20-L)
      elseif((L.GE.0).AND.(L.LT.10)) then
      row = 98+(10-L)
      elseif((L.LT.0).AND.(L.GT.-10)) then
      row = 114+(-L)
      elseif((L.LE.-10).AND.(L.GT.-20)) then
      row = 115+(-L)
      elseif((L.LE.-20).AND.(L.GT.-30)) then
      row = 116+(-L)
      elseif((L.LE.-30).AND.(L.GT.-40)) then
      row = 117+(-L)
      elseif((L.LE.-40).AND.(L.GT.-50)) then
      row = 123+(-L)
      elseif((L.LE.-50).AND.(L.GT.-60)) then
      row = 124+(-L)
      elseif((L.LE.-60).AND.(L.GT.-70)) then
      row = 125+(-L)
      elseif((L.LE.-70).AND.(L.GT.-80)) then
      row = 126+(-L)
      else! if((L.LE.-80).AND.(L.GT.-90)) then
      row = 127+(-L)
      endif
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'row=',row
      VMCRLINE=row
      end function VMCRLINE
!----------------------------------------------------------------------
! Find column of date for specific longitude point
      function VMCRCOL(L)
      IMPLICIT NONE
      INTEGER(4)::L,col,VMCRCOL
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE
      Do while (L.GE.0)
         L=L-20
      End Do
      if (L.LT.0) L=L+20
      col=L+1
      if (DIAGNOSE.EQ.'YES') Write(201,*) 'column=',col
      VMCRCOL=col 
      end function VMCRCOL
!----------------------------------------------------------------------
!     Nonvertical cutoffs
      function NVC_SSS(VC,NZ,NA,doalat,doalong)
      IMPLICIT NONE 
      INTEGER::NZ,NA
! nonvertical cutoffs based on Stormer theory using Smart and Shea formula
! ADVANCES IN SPACE RESEARCH, Vol 32, No. 1, pp 95-102, 2003.  
! Fixed type assignment error that was causing integer-like values 5Jan2015
      REAL(8)::NVC_SSS
      Real(8)::VC,NVC,doalat,doalong,cut, MAXVC
      REAL(8)::PI,rlat,rlon,rmlat,x
      REAL(8)::epsilon,phi,GML
      LOGICAL :: TF
      CHARACTER(1)::NVRC
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE
      COMMON /NVC_TF/TF,NVRC 

      IF ((NVRC.EQ.'1').OR.(NVRC.EQ.'3')) THEN
 
!     VC = vertical cutoff
!     MAXVC = highest VC from tables
!     NVC = nonvertical cutoff
!     doalat = direction of approach, latitude (epsilon=zenith angle) 
!              of approaching particles
!     doalong = direction of approach, longitude (phi=angle from north)
!              of approaching particles
!     GML = Geomagnetic latitude of location in radians
!      phi=azimuth measured clockwise from magnetic north

!     Fortran default functions are in radians so must convert angles
         PI=3.14159265358979                                     
         phi=(doalong)*PI/180.
!         phi=(doalong+90.)*PI/180.
         epsilon=(90.-doalat)*PI/180.
         GML = ACOS((VC/MAXVC(VC))**(.25))
         x=1.-SIN(epsilon)*SIN(phi)*COS(GML)**3
         cut = 4.*VC/(1.+SQRT(x))**2
         NVC_SSS=cut
!        THESIS EQS 16-18 derived from a solution by Shea and Smart to
!        Stormer theory normalizing local nonvertical cutoffs to the local
!        vertical cutoff  
      ELSE 
! If user requested non-vertical transport option '0' or '2' in .INI file.
         NVC_SSS=VC
!        THIS OPTION USES THE VERTICAL CUTOFF FOR THE WHOLE SKY
      ENDIF              
      end function NVC_SSS
!     Nonvertical cutoffs
!---------------------------------------------------------------------72
      function NVC(VC,NZ,NA,doalat,doalong) 
      IMPLICIT NONE
      INTEGER::NZ,NA
      Real(8)::VC,NVC,NVC_SSS,doalat,doalong 

         NVC=NVC_SSS(VC,NZ,NA,doalat,doalong) !Shea and Smart from Stormer

!       OR ALWAYS USES THE VERTICAL CUTOFF FOR THE WHOLE SKY
!       NVC=VC
                    
      end function NVC
!----------------------------------------------------------------------
      function MAXVC(VC) 
      IMPLICIT NONE
      Real(8)::VC,MAXVC 

         MAXVC=15.0 !DEFAULT VALUE from Stormer theory/geomagnetic latitude 
         IF (VC.GT.MAXVC) THEN
            MAXVC=VC !DRIFTS WITH MAGNETIC FIELD EPOCH
         ENDIF
                    
      end function MAXVC
!----------------------------------------------------------------------
      function NVFACTOR(VC) 
! correction factor for dose rate based on cutoff
      IMPLICIT NONE
      Real(8)::VC,NVFACTOR,PI,GML,MAXVC
      LOGICAL :: TF
      CHARACTER(3)::DIAGNOSE
      CHARACTER(1)::NVRC
      COMMON /SKYDIAG/DIAGNOSE
      COMMON /NVC_TF/TF,NVRC 

! nonvertical cutoff reduction of dose rate based on Figure 1 from Felsberger et al, 
! Radiation Protection Dosimetry (2009), Vol. 136, No. 4, pp. 267–273 
! doi:10.1093/rpd/ncp128
! compared with Shea and Smart vertical cutoffs.
!      IF (TF) THEN 
!         NVFACTOR=0.01*(100.-VC*1.5) !original from Felsberger et al, Figure 1
         NVFACTOR=0.01*(100.-VC) ! better fit with ISO spectrum at high VC
!      else
!         NVFACTOR=1. !normal in CARI-7A, which uses a slower, more formal
!                     approach using NVC and Slant subs.
!      endif
            
      end function NVFACTOR
!----------------------------------------------------------------------
      Function SKY(NZ,NA,doalat,altkm)
!     This function calculates the fraction of the celectial sphere 
!     surrounding the aircraft in this segment of sky influenced
!     by the cutoff block in question based on the particle's direction
!     of approach.
!     The sky is broken into NZxNA blocks above the horizon.

      IMPLICIT NONE
      INTEGER::NZ,NA
      REAL(8)::doalat,doalong,f,SKY,zstep,altkm
      REAL(8)::Su,Sl,urlat,lrlat,PI
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE

!      doalat = direction of approach in latitude
!      doalong = direction of approach in longitude
!      f,SKY = fraction of celestial sphere in this segment

      PI=3.14159265358979

      urlat = (PI/180.)*(doalat+zstep(NZ,altkm)) !convert to radians
      lrlat = (PI/180.)*(doalat-zstep(NZ,altkm))
      Su = 0.5*(1-SIN(urlat))  
      Sl = 0.5*(1-SIN(lrlat))
      f=(1./NA)*(Sl-Su)
      SKY=f  
      end function SKY
!----------------------------------------------------------------------
      function Earth(NZ,NA,doalat,altkm)
! Function to calculate the percentage of occultation by the Earth for 
! each block of sky.
! KC 201603 Since I now calculate the angle to horizon and integrate up from it,
! this function ot longer serves it original purpose (the results is always =1),
! but in the future it may be useful for defining anisotropic shields.
   
      IMPLICIT NONE 
      INTEGER::NZ,NA
      Real(8)::doalat,zeta,f,PI,Earth,zstep
      REAL(8)::Sm,Su,Sl,urlat,lrlat,altkm,Re
      CHARACTER(3)::DIAGNOSE
      COMMON /SKYDIAG/DIAGNOSE

!      PI=3.14159265358979
!      
!      Re=6374.1
! If angle of approach is below horizon there is no transmission of 
! particles       
!      urlat = (PI/180.)*(doalat+zstep(NZ,altkm))  
!      lrlat = (PI/180.)*(doalat-zstep(NZ,altkm))
!      if (urlat.le.-zeta(altkm)) then 
!         f=0. ! no transmission through solid Earth
!      elseif (lrlat.gt.-zeta(altkm)) then
!         f=1. ! Earth is not in the way
!      else !partial blockage of access to from   
!         Sm=0.5*(1-Sin(urlat))
!         Su=0.5*(1-Sin(-zeta(altkm)))
!         Sl=0.5*(1-Sin(lrlat))
!         f=(Su-Sl)/(Sm-Sl)
!      endif
!      Earth = f

!     CARI-7 only has 1 section of sky
      Earth=1.  
      end function Earth
!----------------------------------------------------------------------
      FUNCTION zeta(altkm)
      implicit none
      REAL(8)::zeta,A,PI,Re,altkm
!     
! Calculate angle below horizontal to true horizon from vehicle altitude
!
      PI=3.14159265358979
      Re=6374.1
      A=Re+altkm
      zeta=(PI/2.)-ASIN(Re/A) !result in radians

      END FUNCTION 
!----------------------------------------------------------------------
      FUNCTION zstep(NZ,altkm)
      implicit none
      INTEGER::NZ
      REAL(8)::zeta,A,PI,Re,altkm,zd,zstep
!     
! Calculate 1/2 bin width for zenith steps in degrees
      PI=3.14159265358979
      zd=(180./PI)*(zeta(altkm))
      zstep = (90.+zd)/(NZ*2) !should range with altkm from 2.5 to ~2.8 degrees
      END FUNCTION 
!----------------------------------------------------------------------
      function GETKP(dyear)
      IMPLICIT NONE
!
! Function to get Kp index from file based on date and time from dyear 
!  
      REAL(8)::dyear
      INTEGER(4)::yyyy,mm,dd,h,SEARCH,octal,k,avap
      INTEGER(4), DIMENSION(8)::kp, ap
      INTEGER(4)::yr,mo,da,hr,GETKP
      CHARACTER*56::HEADLINE
      CHARACTER*1::porm
      CHARACTER(3)::DIAGNOSE
      LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /SKYDIAG/DIAGNOSE

      IF (MONTHLYAVE.OR.YEARLYAVE.OR.DAILYAVE) THEN !i.e. we have no specific date and time
         GETKP=0 ! Assumes no GEOMAGNETIC STORMS effects in averages
         RETURN
      ENDIF
      
      CALL Dates(dyear,yr,mo,da,hr)
      octal = 1+(h/3) !kp is 3hr index
      OPEN(UNIT=21,FILE='KP_INDEX/KP_INDEX.TXT',STATUS='OLD')
      SEARCH=1
      READ(21,21201) HEADLINE
      DO WHILE ((SEARCH.EQ.1).AND.(.NOT.EOF(21)))
        READ(21,21202) yyyy,mm,dd,porm,kp(1),porm,kp(2),porm,kp(3),porm,&
     &       kp(4),porm,kp(5),porm,kp(6),porm,kp(7),porm,kp(8),         &
     &       ap(1),ap(2),ap(3),ap(4),ap(5),ap(6),ap(7),ap(8),avap
!         Write(201,21202) yyyy,mm,dd,porm,kp(1),porm,kp(2),porm,kp(3),porm,&
!     &       kp(4),porm,kp(5),porm,kp(6),porm,kp(7),porm,kp(8),         &
!     &       ap(1),ap(2),ap(3),ap(4),ap(5),ap(6),ap(7),ap(8),avap
!         Write(201,*) 'testing', yyyy, mm, dd, ' vs ', yr, mo, da
        if ((yyyy.EQ.yr).AND.(mm.EQ.mo).and.(dd.EQ.da)) then
          SEARCH=0
          k=kp(octal)
         IF (DIAGNOSE.EQ.'YES') THEN
          Write(201,21203)'Using Kp=',k,' for ',yyyy,mm,dd,' at ',h,' h'
         endif
        ENDIF
      END DO
      IF (SEARCH.EQ.1) THEN !date not found in file, User for Kp=0
         WRITE(201,*) 'Digiyear', dyear
         WRITE(201,*) 'Year:',yr
         WRITE(201,*) 'Month:',mo
         WRITE(201,*) 'Day:',da
         WRITE(201,*) 'Hour:',h
         Write(201,*) 'Kp index not found in database for this date'
         Write(201,*) 'Using standard, quiet conditions (Kp=0-4)   '
         k=0
      ENDIF
      CLOSE(21)
      GETKP=k
21201  FORMAT(A56)
21202  FORMAT(I4,I2,I2,8(A1,I1),9A3)
21203  FORMAT(A9,I2,A5,I4,2I3,A4,I2,A2)
      end function GETKP
!----------------------------------------------------------------------
      subroutine Dates(dyear,y,m,d,h)
!converts fractional year into y,m,d,h coords
      real(8)::A,dyear,r,diy
      integer(4)::i,y,m,d,h,B
      integer(4)::dimo,tdaeom

      y=INT(dyear)
      r=dyear-y
      m=0
      i=0
      A=diy(y)*r
      B=INT(A)
      h=INT(((A-INT(A))*24.)+0.49)
      DO i=1,12
         if (tdaeom(y,i).GT.A) then 
            m=m+0
         else
            m=i !last month that has fully passed
         endif
      ENDDO
      m=m+1
      d=B-(tdaeom(y,m)-dimo(y,m))
      end subroutine
!----------------------------------------------------------------------
! THIS CONVERTS CALENDAR DATE AND TIME TO YEAR.FRACTION FORM
      subroutine dt2dyear(y,m,d,uth,utm,dyear)
      REAL(8)::dyear,diy
      INTEGER(4)::y,m,d,uth,utm,dimo,tdaeom

      dyear=y+((tdaeom(y,m)-dimo(y,m)+d-1)+(uth/24)+utm/(24*60))/diy(y)
!      PRINT*, dyear, y,m,
      end subroutine
!----------------------------------------------------------------------
      function diy(y) 
! days in year y
      integer(4)::leap,y
      real(8)::diy,dummy
      IF (leap(y).EQ.1) THEN
         dummy=366.
      ELSE
         dummy=365.
      ENDIF
!      Print*, dummy
      diy=dummy
      end function
!----------------------------------------------------------------------
      function dimo(y,m)
! days in month m of year y
      integer(4)::dummy,dimo,m,y

      select case (m) 
      case(1,3,5,7,8,10,12)
        dummy=31
      case(2) 
        if (leap(y).eq.1) then 
          dummy=29
        else
          dummy=28
        endif
      case(4,6,9,11)
        dummy=30
      case DEFAULT
        dummy=0 
      end select
!      Print*, y,m,dummy
      dimo=dummy  
      end function 
!----------------------------------------------------------------------
      function tdaeom(y,m)
      !total days passed in year at the end of month m for year y
      integer(4)::tdaeom,y,m,dummy,leap,dimo,i
        dummy=0
        do i=1,m
           dummy=dummy+dimo(y,i)
        end do
!        Print*, dummy, y, m
        tdaeom=dummy
      end function 
!----------------------------------------------------------------------
      function leap(y)
      integer(4)::y,leap
        IF(MOD(y,400).EQ.0) then
           leap = 1
           RETURN
        ELSEIF (MOD(y,100).EQ.0) THEN
           leap = 0
           RETURN
        ELSEIF(MOD(y,4).EQ.0) THEN
           leap = 1
           RETURN
        ELSE
           leap = 0
        ENDIF
      end function 
!----------------------------------------------------------------------
! This sub loads Forbush decrease data into memory
      SUBROUTINE LOADFORBUSH

      REAL(8),Dimension(1957:2056,12,31,0:23)::ff
      REAL(8)::MONTHLY,HOURLY,RATIO
      
      INTEGER(4)::Y,M,D,H,j 

      CHARACTER(10)::INIVAR
      CHARACTER(12)::INIVAL
      CHARACTER(12)::VIEWER 
      CHARACTER(5)::OS
      CHARACTER(4)::OUTPUT
      CHARACTER(3)::MENUS,DISPLAY,DIAGNOSE
      INTEGER(4)::SP,GCR
 
      COMMON /INIT/MENUS,OS,DISPLAY,DIAGNOSE,VIEWER,OUTPUT
      COMMON /INI/SP,GCR

      COMMON /FDDATA/ff

      ff=1. ! Initialize ff to assumed value of 1.0  
      WRITE(*,*) 'Loading Forbush data'
! Then assign values at dates and times where data are available

      IF ((OS(1:3).EQ.'WIN').OR.(OS(1:3).EQ.'DOS')) THEN
         OPEN(131,file='FORBUSH\FORBUSH.DAT',STATUS='OLD')
      ELSE
         OPEN(131,file='FORBUSH/FORBUSH.DAT',STATUS='OLD')
      ENDIF

      Do while (.NOT.EOF(131))
         READ(131,13101) Y,M,D,H,j,MONTHLY,HOURLY,RATIO
         ff(Y,M,D,H)=RATIO 
!         WRITE(*,*) Y,M,D,H,RATIO
      ENDDO
      CLOSE(131)
13101 format(5I4,F12.4,I6,F12.4)
      END SUBROUTINE 
!
!----------------------------------------------------------------------
      SUBROUTINE FORBUSH(t,fcf,fa)
      IMPLICIT NONE
      REAL(8)::t,MONTHLY,HOURLY,RATIO,fa,TEST,S
      REAL(8), DIMENSION(28,100)::fcf
      REAL(8), DIMENSION(1957:2056,12,31,0:23)::ff
      INTEGER(4)::z,e,i,j,Y,M,D,H,N
      CHARACTER*3::DIAGNOSE
      LOGICAL::YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /USEAVES/YEARLYAVE,MONTHLYAVE,DAILYAVE
      COMMON /SKYDIAG/DIAGNOSE
      COMMON /FDDATA/ff 


!      WRITE(*,*) 'Time into Forbush is:',t
      IF ((t.LT.1958).OR.(t.GT.2056.999)) THEN
         fa=1.0
         WRITE(*,*) 'Date is out of range for Forbush data'
         WRITE(*,*) 'Warning: FORBUSH FACTOR IS ASSUMED AS ',fa
         WRITE(40,*) 'Date is out of range for Forbush data'
         WRITE(40,*) 'Warning: FORBUSH FACTOR IS ASSUMED AS ',fa
         Return
      ENDIF       
!      WRITE(*,*) 'Y:',YEARLYAVE,' M:', MONTHLYAVE, ' D:',DAILYAVE
      IF (MONTHLYAVE.OR.YEARLYAVE) THEN !i.e. we have no specific date and time
         IF (DIAGNOSE.eq.'YES') WRITE(40,*)' Date is an average, no FDs'
         fa=1.0
         fcf=fcf*fa ! Assumes no Forbush effects in averages
         RETURN
      ENDIF
      IF (DAILYAVE) THEN !get data for day
         CALL Dates(t,Y,M,D,H)
         S=0
         fa=1.0
         DO N=0,23
             RATIO=ff(Y,M,D,N)
             IF (DIAGNOSE.eq.'YES') WRITE(40,13102)Y,M,D,N,RATIO
             IF (RATIO.GT.1.3)THEN
                S=S+1.0 !DATA IS BAD FOR THIS TIME
             ELSE
                S=S+RATIO ! directly scale effect with NM shifts
             ENDIF
         ENDDO
         fa=S/24
         IF (DIAGNOSE.eq.'YES') WRITE(40,*)'FORBUSH correction is:',fa
      ELSE ! LOOK UP FORBUSH DECREASE ACTIVITY FOR THIS DATE AND TIME
         CALL Dates(t,Y,M,D,H) 
         RATIO=ff(Y,M,D,H) ! 
         IF (DIAGNOSE.eq.'YES') WRITE(40,13102)Y,M,D,H,RATIO
         IF (RATIO.GT.1.3)THEN
             fa=1.0 !DATA IS BAD FOR THIS TIME, PROBABLE GLE
         ELSE
             fa=RATIO ! directly scale effect with NM shifts
         ENDIF
         IF (DIAGNOSE.eq.'YES') WRITE(40,*)'FORBUSH correction is:',fa
      ENDIF

      fcf=fcf*fa 
13102 format(4I4,F12.4)

      END SUBROUTINE FORBUSH
!                                                                      7
!----6-----------------------------------------------------------------2
!     CONVERSION from geopotential altitude in feet to atmospheric depth 
!     in g/sq.cm. 
!
!     Unit 29 Data from 
!     NOAA,NASA,and USAF. U.S. Standard Atmosphere, 1976. NOAA S/T 76-152, 
!     US GOV PRINTING OFFICE, WASHINTON d.c. 
      SUBROUTINE FT2GPCMS(F,G)
      IMPLICIT NONE        
         INTEGER(4) :: I,J
         REAL(8), INTENT(IN) :: F
         REAL(8), INTENT(OUT) :: G
         REAL(8), DIMENSION(553) :: KFT, GRAMS
         CHARACTER(10)::HEADER1,HEADER2
         CHARACTER(3)::DIAGNOSE='no ' !'YES'
         OPEN(UNIT=29,FILE='ATMOSPHERE/FT-GM.DAT',STATUS='OLD')          
         READ(29,9000) HEADER1, HEADER2
         DO I = 1,553
            READ(29,*) KFT(I), GRAMS(I)
            KFT(I)=KFT(I)*1000 !TABLE IS IN 1000'S OF FEET
            IF (DIAGNOSE.EQ.'YES') WRITE(40,*) F, KFT(I), GRAMS(I)
         ENDDO
         IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'FEET IN:', F
         IF ((F .LE. 400000.).AND.(F.GE.-600)) THEN
            CALL TDSPLINE(553,KFT,GRAMS,1.0E+30,1.0E+30,F,G)
            IF (DIAGNOSE.EQ.'YES') WRITE(40,*)'GRAMS OUT:', G
               CLOSE(29) 
               RETURN
         ELSE !ALTITUDE IS OUT OF RANGE
            CALL EPITATH(' ALTITUDE CONVERSION OUT OF RANGE ',35)
         ENDIF
9000  FORMAT(A10,A10)
      CLOSE(29)
      END SUBROUTINE
!                                                                      7
!----6-----------------------------------------------------------------2
! Assigns corrections of depth to shower data. Two corrections, the 
! first (NMF) for transport without a magnetic field, the second approximates  
! a Chapman function by reassigning showers depths for off-vertically
! approaching showers. For explanations see notes below. 
! 
      FUNCTION NVALT(Zetastep,ALT,I,Z)

      IMPLICIT NONE

      INTEGER(4)::I, Zetastep, Z 
! I is zenith index, Z is nuclear charge, Zetastep is 1/2 zenith stepsize in degrees
      REAL(8)::ALT, NVALT, PI, BETA, ALPHA, PHI, NMF, C, SLANT, SDA
      REAL(8)::ZN
      LOGICAL :: TF
      CHARACTER(3)::DIAGNOSE
      CHARACTER(1)::NVRC
      COMMON /SKYDIAG/DIAGNOSE
      COMMON /NVC_TF/TF,NVRC 

      IF (NVRC.EQ.'0') THEN !no depth corrections, no angular cutoffs
         NVALT=ALT 
         RETURN
      ENDIF
      IF (NVRC.EQ.'1') THEN 
         NVALT=ALT
      ENDIF
      IF (NVRC.EQ.'2') THEN
         NVALT=ALT*SLANT(Zetastep,ALT,I,Z)
      ENDIF
      IF (NVRC.EQ.'3') THEN
         NVALT=ALT*SLANT(Zetastep,ALT,I,Z)  
      ENDIF   

! Add slant depth adjustment, SDA
      PI=3.141592685359
      ALPHA = (2*I-1)*Zetastep 
      !Zetastep = zstep (NZ,km(ALT))
      ! Alpha is effective zenith angle from vertical in degrees by index
      beta = ALPHA*PI/180. ! to radians
! Optional corrections for slant depth effect expect approximately 
! 1/cos(beta) for low beta, less for high beta based on Chapman function.
! I.e., SDA=ABS(COS(beta))**(-PHI) where PHI=1 is simple slant function, PHI = 0.0 
! is no correction for nonvertical
!
      PHI=1.0
      !Next bit to avoids inverse cos(beta) exploding or being negative 
      !close to the horizon 
!      IF (beta.GT.(85*PI/180)) THEN
!         beta=85*PI/180 
!      ENDIF
      SDA = ABS(COS(beta))**(-PHI)
      IF ((NVRC.EQ.'2').OR.(NVRC.EQ.'3')) THEN
         NVALT=SDA*NVALT
      ENDIF

! SDA (slant depth adjustment) The ISO shower represents the effect of beams 
! of numerous depths, with vertical being the shortest component. SDA is an 
! attempt to assign a more correct depth to the ISO shower when using it
! at different slant depths.
!
!      IF (NVALT.GT.1035.0) THEN ! stricly avoids extrapolation in DAP
!         NVALT = 1035.0
!      ENDIF
! KC 20160915 switched back to allowing large depths 
      IF (NVALT.GT.5000.0) THEN ! DAP will extrapolate the high depth 
      ! using eponential attenuation scaling based on ICRU data
         NVALT = 5000.0
      ENDIF

      END FUNCTION 
!                                                                      7
!----6-----------------------------------------------------------------2
      FUNCTION SLANT(Zetastep,ALT,I,Z)

! A slant angle adjustment

      IMPLICIT NONE

      INTEGER(4)::I, Zetastep, Z 
! I is zenith index, NZ is zenith steps, Z is nuclear charge
      REAL(8)::PI, SR, ALPHA, C, SLANT, ALT
      LOGICAL :: TF
      CHARACTER(3)::DIAGNOSE
      CHARACTER(1)::NVRC
      COMMON /SKYDIAG/DIAGNOSE
      COMMON /NVC_TF/TF,NVRC 

      SLANT=1.0 !1=NO MODULATION

!      PI=3.141592685359
!      ALPHA = (2*I-1)*Zetastep! effective angle from vertical in degrees by index
!      SR = ALPHA*PI/180. ! to radians
!      SLANT=0.5+0.5*ABS(COS(SR)) !flux renormalization based on MCNPX F2 tally normalization

!      
! The ISO shower represents the effect of beams 
! of numerous depths, with vertical being the shortest component. This is an 
! attempt to assign a more correct intensity to the ISO shower when using it
! at different slant angles.

      END FUNCTION 
!                                                                      7
!----6-----------------------------------------------------------------2

