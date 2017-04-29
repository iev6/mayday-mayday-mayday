! Adapted for CARI by KC from codes freely distributed by the U.S. NOAA- 
! NGS (National Geodetic Survey)
! https://www.ngs.noaa.gov/
!
! b::forward
!
      subroutine forward(slat, slon, foraz, edist, elat, elon)
!
!********1*********2*********3*********4*********5*********6*********7**
!
! name:      forward
! version:   200208.19 
! author:    stephen j. frakes
! purpose:   to compute a geodetic forward (direct problem)
!            and then display output information
!
! input parameters:
! -----------------
!
!c output parameters:
!c ------------------
!c
!c local variables and constants:
!c ------------------------------
!c answer           user prompt response
!c arc              meridional arc distance latitude p1 to p2 (meters)      
!c b                semiminor axis polar (in meters)
!c baz              azimuth back (in radians)
!c blimit           geodetic distance allowed on ellipsoid (in meters)
!c buff             input char buffer array
!c dd,dm,ds         temporary values for degrees, minutes, seconds
!c dmt,d_mt         char constants for units (in meters)        
!c dd_max           maximum ellipsoid distance -1 (in meters)        
!c edist            ellipsoid distance (in meters)
!c elips            ellipsoid choice
!c esq              eccentricity squared for reference ellipsoid
!c faz              azimuth forward (in radians)
!c filout           output file name
!c finv             reciprocal flattening
!c hem              hemisphere flag for lat & lon entry  
!c ierror           error condition flag with d,m,s conversion
!c lgh              length of buff() array
!c option           user prompt response             
!c r1,r2            temporary variables    
!c ss               temporary value for ellipsoid distance
!c tol              tolerance for conversion of seconds    
!c
!c name1            name of station one
!c ld1,lm1,sl1      latitude  sta one - degrees,minutes,seconds
!c ald1,alm1,sl1    latitude  sta one - degrees,minutes,seconds
!c lat1sn           latitude  sta one - sign (+/- 1)
!c d_ns1            latitude  sta one - char ('N','S')
!c lod1,lom1,slo1   longitude sta one - degrees,minutes,seconds
!c alod1,alom1,slo1 longitude sta one - degrees,minutes,seconds
!c lon1sn           longitude sta one - sign (+/- 1)
!c d_ew1            longitude sta one - char ('E','W')
!c iaz1,maz1,saz1   forward azimuth   - degrees,minutes,seconds
!c isign1           forward azimuth   - flag  (+/- 1)
!c azd1,azm1,saz1   forward azimuth   - degrees,minutes,seconds
!c iazsn            forward azimuth   - flag  (+/- 1)
!c glat1,glon1      station one       - (lat & lon in radians )
!c
!c name2            name of station two
!c ld2,lm2,sl2      latitude  sta two - degrees,minutes,seconds
!c lat2sn           latitude  sta two - sign (+/- 1)
!! d_ns2            latitude  sta two - char ('N','S')
! lod2,lom2,slo2   longitude sta two - degrees,minutes,seconds
! lon2sn           longitude sta two - sign (+/- 1)
! d_ew2            longitude sta two - char ('E','W')
! iaz2,maz2,saz2   back azimuth      - degrees,minutes,seconds
! isign2           back azimuth      - flag  (+/- 1)
! glat2,glon2      station two       - (lat & lon in radians )
!
! global variables and constants:
! -------------------------------
! a                semimajor axis equatorial (in meters)
! f                flattening
! pi               constant 3.14159....
! rad              constant 180.0/pi  
!
!    this module called by:  n/a
!
!    this module calls:      elipss, getrad, dirct1, todmsp
!    gethem, trimm,   bufdms, gvalr8, gvali4, fixdms, gpnarc
!    datan,  write,  read,   dabs,   open,   stop
!
!    include files used:     n/a
!
!    common blocks used:     const, elipsoid
!
!    references:             see comments within subroutines
!
!    comments:
!
!********1*********2*********3*********4*********5*********6*********7**
!::modification history
!::1990mm.dd, sjf, creation of program           
!::199412.15, bmt, creation of program on viper
!::200203.08, crs, modified by c.schwarz to correct spelling of Clarke
!::                at request of Dave Doyle
!::200207.18, rws, modified i/o & standardized program documentation
!::                added subs trimm, bufdms, gethem, gvali4, gvalr8      
!::200207.23, rws, added sub gpnarc
!::200208.15, rws, fixed an error in bufdms
!::              - renamed ellips to elipss "common error" with dirct2
!::              - added FAZ & BAZ to printed output 
!::200208.19, rws, added more error flags for web interface 
!::              - added logical nowebb                      
!::200208.xx, sjf, program version number 2.0                   
!********1*********2*********3*********4*********5*********6*********7**
!e::forward
!
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
!
      logical  nowebb
!
      character*1  answer,option,dmt,buff(50),hem
      character*6  d_ns1, d_ew1, d_ns2, d_ew2, d_mt
      character*30 filout,name1,name2,elips
!
      integer*4    ierror
      integer*4    lgh
!
      common/const/pi,rad
      common/elipsoid/a,f
!
!     ms_unix      0 = web version
!                  1 = ms_dos or unix
!
      parameter   ( ms_unix = 0 )
!
      pi=4.d0*datan(1.d0)
      rad=180.d0/pi
      dmt='m'
      d_mt='Meters'
!
      if (ms_unix.eq.1) then
        nowebb = .true.
      else
        nowebb = .false.
      endif
!
        option='1'
        a=6378137.d0
        f=1.d0/298.25722210088d0
        elips='GRS80 / WGS84  (NAD83)'
!
      esq = f*(2.0d0-f)
!
      r1  = 0.0d0
      r2  = pi/2.0d0
      call gpnarc ( a, f, esq, pi, r1, r2, arc )
!
!     compute the geodetic limit distance (blimit), it is equal
!     to twice the distance to the pole minus one meter
!
      blimit = 2.0d0*arc-1.0d0
!
!     maximum distance allowed on ellipsoid
!
      dd_max = blimit                  

      glat1=slat/rad
      glon1=slon/rad
      faz=foraz !/rad
      ss = edist
      ss = dabs(ss)
!
      if (ss.lt.dd_max) then
        edist  = ss
        irdst1 = 0
      else
        irdst1 = 1
! This should never be called during a flight! 
!        write(*,*) ' Invalid Distance ... Please re-enter it '
!        write(*,*) '  '
        edist  = 0.001d0
      endif
!
      call dirct1 (glat1,glon1,glat2,glon2,faz,baz,edist)
      elat=glat2*rad
      elon=glon2*rad
!
      end

      subroutine bufdms (buff,lgh,hem,dd,dm,ds,ierror)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
!
      logical     done,flag
!
      character*1 buff(*),abuf(21)
      character*1 ch
      character*1 hem
      integer*4   ll,lgh
      integer*4   i4,id,im,is,icond,ierror
      real*8      x(5)
!
!     set the "error flag" 
!
      ierror = 0
      icond  = 0
!
!     set defaults for dd,dm,ds
!
      dd = 0.0d0
      dm = 0.0d0
      ds = 0.0d0
!
!     set default limits for "hem" flag
!
      if     (hem.eq.'N' .or. hem.eq.'S' ) then
        ddmax = 90.0d0
      elseif (hem.eq.'E' .or. hem.eq.'W' ) then
        ddmax = 360.0d0
      elseif (hem.eq.'A' ) then
        ddmax = 360.0d0
      elseif (hem.eq.'Z' ) then
        ddmax = 180.0d0
      elseif( hem.eq.'*' ) then
        ddmax  = 0.0d0
        ierror = 1
      else
        ddmax = 360.0d0
      endif
!
      do 1 i=1,5
        x(i) = 0.0d0
    1 continue
!
      icolon = 0
      ipoint = 0
      icount = 0
      flag   = .true.
      jlgh   = lgh
!
      do 2 i=1,jlgh
        if( buff(i).eq.':' ) then
          icolon = icolon+1
        endif
        if( buff(i).eq.'.' ) then
          ipoint = ipoint+1
          flag   = .false.
        endif
        if( flag ) then
          icount = icount+1
        endif
    2 continue
!
      if( ipoint.eq.1 .and. icolon.eq.0 ) then
!
!       load temp buffer
!
        do 3 i=1,jlgh
          abuf(i) = buff(i)
    3   continue
        abuf(jlgh+1) = '$'
        ll = jlgh
!
        call gvalr8 (abuf,ll,r8,icond)
!
        if( icount.ge.5 ) then
!
!         value is a packed decimal of ==>  DDMMSS.sssss       
!
          ss = r8/10000.0d0
          id = idint( ss )
!
          r8 = r8-10000.0d0*dble(float(id))
          ss = r8/100.0d0
          im = idint( ss )
!
          r8 = r8-100.0d0*dble(float(im))
        else
!
!         value is a decimal of ==>  .xx   X.xxx   X.  
!
          id = idint( r8 )
          r8 = (r8-id)*60.0d0
          im = idint( r8 )
          r8 = (r8-im)*60.0d0
        endif
!
!       account for rounding error
!
        is = idnint( r8*1.0d5 )
        if( is.ge.6000000 ) then
           r8 = 0.0d0
           im = im+1
        endif
!
        if( im.ge.60 ) then
          im = 0
          id = id+1
        endif
!
        dd = dble( float( id ) )
        dm = dble( float( im ) )
        ds = r8
      else
!
!       buff() value is a d,m,s of ==>  NN:NN:XX.xxx    
!
        k    = 0
        next = 1
        done = .false.
        ie   = jlgh
!
        do 100 j=1,5
          ib = next
          do 90 i=ib,ie
            ch   = buff(i)
            last = i
            if( i.eq.jlgh .or. ch.eq.':' ) then
              if( i.eq.jlgh ) then
                done = .true.
              endif
              if( ch.eq.':' ) then
                last = i-1
              endif
              goto 91
            endif
   90     continue
          goto 98
!
   91     ipoint = 0
          ik     = 0
          do 92 i=next,last
            ik = ik+1
            ch = buff(i)
            if( ch.eq.'.' ) then
              ipoint = ipoint+1
            endif
            abuf(ik) = buff(i) 
   92     continue
          abuf(ik+1) = '$' 
!
          ll = ik
          if( ipoint.eq.0 ) then
            call gvali4 (abuf,ll,i4,icond)
            r8 = dble(float( i4 )) 
          else
            call gvalr8 (abuf,ll,r8,icond)
          endif
!
          k    = k+1
          x(k) = r8
!
   98     if( done ) then
            goto 101
          endif
!
          next = last
   99     next = next+1     
          if( buff(next).eq.':' ) then
            goto 99
          endif
  100   continue
!
!       load dd,dm,ds
!
  101   if( k.ge.1 ) then
          dd = x(1)
        endif
!
        if( k.ge.2 ) then
          dm = x(2)
        endif
!
        if( k.ge.3 ) then
          ds = x(3)
        endif
      endif
!
      if( dd.gt.ddmax  .or. dm.ge.60.0d0 .or. ds.ge.60.0d0 ) then
        ierror = 1
        dd = 0.0d0
        dm = 0.0d0
        ds = 0.0d0
      endif
!
      if( icond.ne.0 ) then
        ierror = 1
      endif
!
      return
      end


      SUBROUTINE DIRCT1(GLAT1,GLON1,GLAT2,GLON2,FAZ,BAZ,S)
!
! *** SOLUTION OF THE GEODETIC DIRECT PROBLEM AFTER T.VINCENTY
! *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS
! *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL
!
! *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID
! *** F IS THE FLATTENING OF THE REFERENCE ELLIPSOID
! *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST
! *** AZIMUTHS IN RADIANS CLOCKWISE FROM NORTH
! *** GEODESIC DISTANCE S ASSUMED IN UNITS OF SEMI-MAJOR AXIS A
!
! *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 20FEB75
! *** MODIFIED FOR SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 750608
!
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CONST/PI,RAD
      COMMON/ELIPSOID/A,F
      DATA EPS/0.5D-13/
      R=1.-F
      TU=R*DSIN(GLAT1)/DCOS(GLAT1)
      SF=DSIN(FAZ)
      CF=DCOS(FAZ)
      BAZ=0.
      IF(CF.NE.0.) BAZ=DATAN2(TU,CF)*2.
      CU=1./DSQRT(TU*TU+1.)
      SU=TU*CU
      SA=CU*SF
      C2A=-SA*SA+1.
      X=DSQRT((1./R/R-1.)*C2A+1.)+1.
      X=(X-2.)/X
      C=1.-X
      C=(X*X/4.+1)/C
      D=(0.375D0*X*X-1.)*X
      TU=S/R/A/C
      Y=TU
  100 SY=DSIN(Y)
      CY=DCOS(Y)
      CZ=DCOS(BAZ+Y)
      E=CZ*CZ*2.-1.
      C=Y
      X=E*CY
      Y=E+E-1.
      Y=(((SY*SY*4.-3.)*Y*CZ*D/6.+X)*D/4.-CZ)*SY*D+TU
      IF(DABS(Y-C).GT.EPS)GO TO 100
      BAZ=CU*CY*CF-SU*SY
      C=R*DSQRT(SA*SA+BAZ*BAZ)
      D=SU*CY+CU*SY*CF
      GLAT2=DATAN2(D,C)
      C=CU*CY-SU*SY*CF
      X=DATAN2(SY*SF,C)
      C=((-3.*C2A+4.)*F+4.)*C2A*F/16.
      D=((E*CY*C+CZ)*SY*C+Y)*SA
      GLON2=GLON1+X-(1.-C)*D*F
      BAZ=DATAN2(SA,BAZ)+PI
      RETURN
      END

!B::GPNARC
!
      SUBROUTINE GPNARC (AMAX,FLAT,ESQ,PI,P1,P2,ARC)
!
!********1*********2*********3*********4*********5*********6*********7*
!
! NAME:        GPNARC
! VERSION:     200005.26
! WRITTEN BY:  ROBERT (Sid) SAFFORD
! PURPOSE:     SUBROUTINE TO COMPUTE THE LENGTH OF A MERIDIONAL ARC 
!              BETWEEN TWO LATITUDES
!
! INPUT PARAMETERS:
! -----------------
! AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
! FLAT         FLATTENING (0.0033528 ... )
! ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
! PI           3.14159...
! P1           LAT STATION 1
! P2           LAT STATION 2
!
! OUTPUT PARAMETERS:
! ------------------
! ARC          GEODETIC DISTANCE 
!
! LOCAL VARIABLES AND CONSTANTS:
! ------------------------------
! GLOBAL VARIABLES AND CONSTANTS:
! -------------------------------
!
!    MODULE CALLED BY:    GENERAL 
!
!    THIS MODULE CALLS:   
!       LLIBFORE/ OPEN,   CLOSE,  READ,   WRITE,  INQUIRE
!                 DABS,   DBLE,   FLOAT,  IABS,   CHAR,   ICHAR
!
!    INCLUDE FILES USED:
!    COMMON BLOCKS USED:  
!
!    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
!                MS-DOS Operating System
!    COMMENTS:
!********1*********2*********3*********4*********5*********6*********7*
!!:MODIFICATION HISTORY
!::197507.05, RWS, VER 00 TENCOL RELEASED FOR FIELD USE
!::198311.20, RWS, VER 01 MTEN   RELEASED TO FIELD
!::198411.26, RWS, VER 07 MTEN2  RELEASED TO FIELD
!::1985xx.xx, RWS, CODE   CREATED               
!::198506.10, RWS, WRK    ENHANCEMENTS RELEASED TO FIELD
!::198509.01, RWS, VER 11 MTEN3  RELEASED TO FIELD
!::198512.18, RWS, CODE   MODIFIED FOR MTEN3
!::198708.10, RWS, CODE   MODIFIED TO USE NEW MTEN4 GPN RECORD FORMAT
!::199112.31, RWS, VER 20 MTEN4 RELEASED TO FIELD
!::200001.13, RWS, VER 21 MTEN4 RELEASED TO FIELD
!::200005.26, RWS, CODE   RESTRUCTURED & DOCUMENTATION ADDED             
!::200012.31, RWS, VER 23 MTEN5 RELEASED                                 
!********1*********2*********3*********4*********5*********6*********7*
!E::GPNARC
! ---------------------------
!     M T E N  (VERSION 3)
!     M T E N  (VERSION 5.23)
! ---------------------------
! 
      IMPLICIT REAL*8 (A-H,O-Z)
!
      LOGICAL  FLAG
!
      DATA TT/5.0D-15/
!
!     CHECK FOR A 90 DEGREE LOOKUP
!
      FLAG = .FALSE.
!
      S1 = DABS(P1)
      S2 = DABS(P2)
!
      IF( (PI/2.0D0-TT).LT.S2 .AND. S2.LT.(PI/2.0D0+TT) ) then
        FLAG = .TRUE.
      ENDIF
!
      IF( S1.GT.TT ) then
        FLAG = .FALSE.
      ENDIF
!
      DA = (P2-P1)
      S1 = 0.0D0
      S2 = 0.0D0
!
!     COMPUTE THE LENGTH OF A MERIDIONAL ARC BETWEEN TWO LATITUDES
!
      E2 = ESQ
      E4 = E2*E2
      E6 = E4*E2
      E8 = E6*E2
      EX = E8*E2
!
      T1 = E2*(003.0D0/4.0D0)
      T2 = E4*(015.0D0/64.0D0)
      T3 = E6*(035.0D0/512.0D0)
      T4 = E8*(315.0D0/16384.0D0)
      T5 = EX*(693.0D0/131072.0D0)
!
      A  = 1.0D0+T1+3.0D0*T2+10.0D0*T3+35.0D0*T4+126.0D0*T5
!
      IF( FLAG ) then
        GOTO 1
      ENDIF
!
      B  = T1+4.0D0*T2+15.0D0*T3+56.0D0*T4+210.0D0*T5
      C  = T2+06.0D0*T3+28.0D0*T4+120.0D0*T5
      D  = T3+08.0D0*T4+045.0D0*T5
      E  = T4+010.0D0*T5
      F  = T5
!
      DB = DSIN(P2*2.0D0)-DSIN(P1*2.0D0)
      DC = DSIN(P2*4.0D0)-DSIN(P1*4.0D0)
      DD = DSIN(P2*6.0D0)-DSIN(P1*6.0D0)
      DE = DSIN(P2*8.0D0)-DSIN(P1*8.0D0)
      DF = DSIN(P2*10.0D0)-DSIN(P1*10.0D0)
!
!     COMPUTE THE S2 PART OF THE SERIES EXPANSION
!
      S2 = -DB*B/2.0D0+DC*C/4.0D0-DD*D/6.0D0+DE*E/8.0D0-DF*F/10.0D0
!
!     COMPUTE THE S1 PART OF THE SERIES EXPANSION
!
    1 S1 = DA*A
!
!     COMPUTE THE ARC LENGTH
!
      ARC = AMAX*(1.0D0-ESQ)*(S1+S2)
!
      RETURN
      END

      subroutine gvali4 (buff,ll,vali4,icond)
      implicit     integer (i-n)
!
      logical      plus,sign,done,error
      character*1  buff(*)
      character*1  ch
!
!     integer*2    i
!     integer*2    l1
!
      integer*4    ich,icond
      integer*4    ll    
      integer*4    vali4
!
      l1    = ll
      vali4 = 0
      icond = 0
      plus  = .true.
      sign  = .false.
      done  = .false.
      error = .false.
!
      i = 0
   10 i = i+1
      if( i.gt.l1 .or. done ) then
        go to 1000
      else
        ch  = buff(i)
        ich = ichar( buff(i) )
      endif
!
      if(     ch.eq.'+' ) then
!
!       enter on plus sign
!
        if( sign ) then
          goto 150
        else 
          sign = .true.
          goto 10
        endif
      elseif( ch.eq.'-' ) then
!
!       enter on minus sign
!
        if( sign ) then
          goto 150
        else
          sign = .true.
          plus = .false.
          goto 10
        endif
      elseif( ch.ge.'0' .and. ch.le.'9' ) then
        goto 100
      elseif( ch.eq.' ' ) then
!
!       enter on space -- ignore leading spaces
!
        if( .not.sign ) then
          goto 10
        else
          buff(i) = '0'
          ich = 48
          goto 100
        endif
      elseif( ch.eq.':' ) then
!
!       enter on colon -- ignore 
!
        if( .not.sign ) then
          goto 10
        else
          goto 1000
        endif
      elseif( ch.eq.'$' ) then
!
!       enter on dollar "$"      
!
        done = .true.
        goto 10
      else
!
!       something wrong
!
        goto 150
      endif
!
!     enter on numeric
!
  100 vali4 = 10*vali4+(ich-48)
      sign  = .true.
      goto 10
!
!     treat illegal character
!
  150 buff(i) = '0'
      vali4 = 0
      icond = 1
!
 1000 if( .not.plus ) then
        vali4 = -vali4
      endif
!
      return
      end
      subroutine gvalr8 (buff,ll,valr8,icond)
      implicit     integer (i-n)
!
      logical      plus,sign,dpoint,done
!
      character*1  buff(*)
      character*1  ch
!
!     integer*2    i, ip
!     integer*2    l1
!     integer*2    nn, num, n48
!
      integer*4    ich,icond
      integer*4    ll
!
      real*8       ten
      real*8       valr8
      real*8       zero
!
      data zero,ten/0.0d0,10.0d0/
!
      n48     =  48
      l1      =  ll
      icond   =   0
      valr8   =  zero  
      plus    = .true.
      sign    = .false.
      dpoint  = .false.
      done    = .false.
!
!     start loop thru buffer
!
      i = 0
   10 i = i+1
      if( i.gt.l1 .or. done ) then
        go to 1000
      else 
        ch  = buff(i)
        nn  = ichar( ch )
        ich = nn
      endif 
!
      if(     ch.eq.'+' ) then
!
!       enter on plus sign
!
        if( sign ) then
          goto 150
        else
          sign = .true.
          goto 10
        endif
      elseif( ch.eq.'-' ) then
!
!       enter on minus sign
!
        if( sign ) then
          goto 150
        else
          sign = .true.
          plus = .false.
          goto 10
        endif
      elseif( ch.eq.'.' ) then
!
!       enter on decimal point
!
        ip     = 0
        sign   = .true.
        dpoint = .true.
        goto 10
      elseif( ch.ge.'0' .and. ch.le.'9' ) then
        goto 100
      elseif( ch.eq.' ' ) then
!
!       enter on space
!
        if( .not.sign ) then
          goto 10
        else
          buff(i) = '0'
          ich = 48
          goto 100
        endif
      elseif( ch.eq.':' .or. ch.eq.'$' ) then
!
!       enter on colon or "$" sign
!
        done = .true.
        goto 10
      else
!
!       something wrong
!
        goto 150
      endif
!
!     enter on numeric
!
  100 sign = .true.
      if( dpoint ) then
        ip = ip+1
      endif
!
      num   = ich
      valr8 = ten*valr8+dble(float( num-n48 ))
      goto 10
!
!     treat illegal character
!
  150 buff(i) = '0'
      valr8   =  0.0d0
      icond   =  1
!
 1000 if( dpoint ) then
        valr8 =  valr8/(ten**ip)
      endif
!
      if( .not.plus ) then
        valr8 = -valr8
      endif
!
      return
      end

      subroutine todmsp(val,id,im,s,isign)
 
!** convert position radians to deg,min,sec
!** range is [-pi to +pi]
 
      implicit double precision(a-h,o-z)
      common/const/pi,rad
 
    1 if(val.gt.pi) then
        val=val-pi-pi
        go to 1
      endif
 
    2 if(val.lt.-pi) then
        val=val+pi+pi
        go to 2
      endif
 
      if(val.lt.0.d0) then
        isign=-1
      else
        isign=+1
      endif
 
      s=dabs(val*rad)
      id=idint(s)
      s=(s-id)*60.d0
      im=idint(s)
      s=(s-im)*60.d0
 
!** account for rounding error
 
      is=idnint(s*1.d5)
      if(is.ge.6000000) then
        s=0.d0
        im=im+1
      endif
      if(im.ge.60) then
        im=0
        id=id+1
      endif
 
      return
      end


!b::inverse
!
!      PROGRAM HOW_FAR
!     front end for inverse added and internal menus stripped by 
!     Kyle Copeland for use by CARI
!
!     Program now returns distance and azimuths given origin 
!     and destination latitudes and longitudes
! 
!      double precision::LA1,LO1,LA2,LO2,FAZ,BAZ,NMILES,meters
!
!      NMILES        nautical miles (6080 ft, 1853 meters)
!      
!      LA1=30.0d0
!      LA2=26.0d0
!      LO1=-105.0d0
!      LO2=155.0d0
!      CALL INVERSE(LA1,LO1,LA2,LO2,FAZ,BAZ,meters)
!      NMILES=meters/1853.0d0
!      print*, LA1,' N ', LO1, ' E to '
!      print*, LA2, ' N ',LO2, ' E is'
!      print*, NMILES, ' nautical miles '
!      print*, FAZ, ' forward azimuth '
!      
!      END PROGRAM HOW_FAR
      SUBROUTINE inverse (gla1,glo1,gla2,glo2,faz,baz,edist)
!
!********1*********2*********3*********4*********5*********6*********7**
!
! name:      inverse
! version:   200208.19
! author:    stephen j. frakes
! purpose:   to compute a geodetic inverse  
!            and then display output information
!
! input parameters:
! -----------------
!
! output parameters:
! ------------------
!
! local variables and constants:
! ------------------------------
! answer           user prompt response
! b                semiminor axis polar (in meters)
! baz              azimuth back (in radians)
! buff             input char buffer array
! dd,dm,ds         temporary values for degrees, minutes, seconds
! dlon             temporary value for difference in longitude (radians)   
! dmt,d_mt         char constants for meter units         
! edist            ellipsoid distance (in meters)
! elips            ellipsoid choice
! esq              eccentricity squared for reference ellipsoid
! faz              azimuth forward (in radians)
! filout           output file name
! finv             reciprocal flattening
! hem              hemisphere flag for lat & lon entry  
! ierror           error condition flag with d,m,s conversion
! lgh              length of buff() array
! option           user prompt response             
! r1,r2            temporary variables    
! ss               temporary variable     
! tol              tolerance for conversion of seconds
!
! name1            name of station one
! ld1,lm1,sl1      latitude  sta one - degrees,minutes,seconds
! ald1,alm1,sl1    latitude  sta one - degrees,minutes,seconds
! lat1sn           latitude  sta one - sign (+/- 1)
! d_ns1            latitude  sta one - char ('N','S')
! lod1,lom1,slo1   longitude sta one - degrees,minutes,seconds
! alod1,alom1,slo1 longitude sta one - degrees,minutes,seconds
! lon1sn           longitude sta one - sign (+/- 1)
! d_ew1            longitude sta one - char ('E','W')
! iaz1,maz1,saz1   forward azimuth   - degrees,minutes,seconds
! isign1           forward azimuth   - flag  (+/- 1)
! glat1,glon1      station one       - (lat & lon in radians )
! p1,e1            standpoint one    - (lat & lon in radians )
!
! name2            name of station two
! ld2,lm2,sl2      latitude  sta two - degrees,minutes,seconds
! ald2,alm2,sl2    latitude  sta two - degrees,minutes,seconds
! lat2sn           latitude  sta two - sign (+/- 1)
! d_ns2            latitude  sta two - char ('N','S')
! lod2,lom2,slo2   longitude sta two - degrees,minutes,seconds
! alod2,alom2,slo2 longitude sta one - degrees,minutes,seconds
! lon2sn           longitude sta two - sign (+/- 1)
! d_ew2            longitude sta two - char ('E','W')
! iaz2,maz2,saz2   back azimuth      - degrees,minutes,seconds
! isign2           back azimuth      - flag  (+/- 1)
! glat2,glon2      station two       - (lat & lon in radians )
! p2,e2            forepoint two     - (lat & lon in radians )
!
! global variables and constants:
! -------------------------------
! a                semimajor axis equatorial (in meters)
! f                flattening
! pi               constant 3.14159....
! rad              constant 180.0/pi  
!
!    this module called by:  n/a
!
!    this module calls:      elipss, getrad, inver1, todmsp
!    gethem, trimm,   bufdms, gvalr8, gvali4, fixdms, gpnhri
!    datan,  write,  read,   dabs,   open,   stop
!
!    include files used:     n/a
!
!    common blocks used:     const, elipsoid
!
!    references:             see comments within subroutines
!
!    comments:
!
!********1*********2*********3*********4*********5*********6*********7**
!::modification history
!::1990mm.dd, sjf, creation of program           
!::199412.15, bmt, creation of program on viper
!::200203.08, crs, modified by c.schwarz to correct spelling of Clarke
!::200207.15, rws, modified i/o & standardized program documentation
!::                added subs trimm, bufdms, gethem, gvali4, gvalr8      
!::200207.23, rws, replaced sub inver1 with gpnarc, gpnloa, gpnhri
!::200208.15, rws, fixed an error in bufdms
!::              - renamed ellips to elipss "common error" with dirct2
!::              - added FAZ & BAZ to printed output      
!::200208.19, rws, added more error flags for web interface code
!::              - added logical nowebb                             
!::200208.xx, sjf, program version number 2.0 
!::201206.20, kac, removed duplicate subs from inverse+forward 
!                  collection now called GEODESIC.FOR
!                 -renamed 'trim' to 'trimm' to avoid duplicating name
!                  of the intrinsic fortran function TRIM
!********1*********2*********3*********4*********5*********6*********7**
!e::inverse
!
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
!
      logical  nowebb
!
      character*1  answer,option,dmt,buff(50),hem
      character*6  d_ns1, d_ew1, d_ns2, d_ew2, d_mt
      character*30 filout,name1,name2,elips
!
      integer*4    ierror
      integer*4    lgh
!
      common/const/pi,rad
      common/elipsoid/a,f
!
!     ms_unix      0 = web version
!                  1 = ms_dos or unix version
!
      parameter   ( ms_unix = 0 )
!
      pi   = 4.d0*datan(1.d0)
      rad  = 180.d0/pi
      dmt  = 'm'
      d_mt = 'Meters'
!
      if( ms_unix.eq.1 ) then
        nowebb = .true. 
      else
        nowebb = .false.
      endif
!
!      write(*,*) '  Ellipsoid options : '
!      write(*,*) '  '
!      write(*,*) '  1) GRS80 / WGS84  (NAD83) '
!      write(*,*) '  2) Clarke 1866    (NAD27) '
!      write(*,*) '  3) Any other ellipsoid '
!      write(*,*) '  '
!      write(*,*) '  Enter choice : '
      option = '1'
!
        a=6378137.d0
        f=1.d0/298.25722210088d0
        elips='GRS80 / WGS84  (NAD83)'
!
      esq = f*(2.0d0-f)
!
      glat1 = gla1/rad
      glon1 = glo1/rad
      glat2 = gla2/rad
      glon2 = glo2/rad
!
      p1 = glat1
      e1 = glon1
      p2 = glat2
      e2 = glon2
!
      if( e1.lt.0.0d0 ) then
        e1 = e1+2.0d0*pi
      endif
      if( e2.lt.0.0d0 ) then
        e2 = e2+2.0d0*pi
      endif
!
!     compute the geodetic inverse
!
! ************************************************************
! *   replaced subroutine inver1 with gpnhri
! *  
! *   call inver1 (glat1,glon1,glat2,glon2,faz,baz,edist)
! *
! ************************************************************
!
      call gpnhri (a,f,esq,pi,p1,e1,p2,e2,faz,baz,edist)
!      print*, p1,' N ', e1, ' E to '
!      print*, p2, ' N ',e2, ' E is'
!      print*, edist, ' meters '
!
!     check for a non distance ... p1,e1 & p2,e2 equal zero ?
!
      if( edist.lt.0.00005d0 ) then
        faz = 0.0d0
        baz = 0.0d0
      endif
! 
!     stop
      end

!cb::gpnhri
!
      subroutine gpnhri (a,f,esq,pi,p1,e1,p2,e2,az1,az2,s)      
!!
!********1*********2*********3*********4*********5*********6*********7*
!
! name:        gpnhri
! version:     200208.09
! written by:  robert (sid) safford
! purpose:     subroutine to compute helmert rainsford inverse problem 
! 
!     solution of the geodetic inverse problem after t. vincenty
!     modified rainsford's method with helmert's elliptical terms
!     effective in any azimuth and at any distance short of antipocal
!     from/to stations must not be the geographic pole.
!     parameter a is the semi-major axis of the reference ellipsoid
!     finv=1/f is the inverse flattening of the reference ellipsoid
!     latitudes and longitudes in radians positive north and west
!     forward and back azimuths returned in radians clockwise from south
!     geodesic distance s returned in units of semi-major axis a
!     programmed for ibm 360-195   09/23/75
!
!     note - note - note -
!     1. do not use for meridional arcs and be careful on the equator.
!     2. azimuths are from north(+) clockwise and 
!     3. longitudes are positive east(+) 
!
! input parameters:
! -----------------
! a            semi-major axis of reference ellipsoid      meters
! f            flattening (0.0033528...)
! esq          eccentricity squared 
! pi           3.14159...
! p1           lat station 1                               radians
! e1           lon station 1                               radians
! p2           lat station 2                               radians
! e2           lon station 2                               radians
!
! output parameters:
! ------------------
! az1          azi at sta 1 -> sta 2                       radians
! az2          azi at sta 2 -> sta 1                       radians
! s            geodetic dist between sta(s) 1 & 2          meters
!
! local variables and constants:
! ------------------------------
! aa               constant from subroutine gpnloa                    
! alimit           equatorial arc distance along the equator   (radians)
! arc              meridional arc distance latitude p1 to p2 (in meters)      
! az1              azimuth forward                          (in radians)
! az2              azimuth back                             (in radians)
! bb               constant from subroutine gpnloa                    
! dlon             temporary value for difference in longitude (radians)   
! equ              equatorial distance                       (in meters)
! r1,r2            temporary variables    
! s                ellipsoid distance                        (in meters)
! sms              equatorial - geodesic distance (S - s) "Sms"       
! ss               temporary variable     
! tol0             tolerance for checking computation value         
! tol1             tolerance for checking a real zero value         
! tol2             tolerance for close to zero value  
! twopi            two times constant pi               
!
! global variables and constants:
! -------------------------------
!
!    module called by:    general 
!
!    this module calls:   gpnarc, gpnloa
!       llibfore/ dsin,   dcos,   dsqrt,  dabs,  datan2, write
!
!    include files used:
!    common blocks used:  
!
!    references: microsoft fortran 4.10 optimizing compiler, 1988
!                ms-dos operating system
!    comments:
!********1*********2*********3*********4*********5*********6*********7*
!::modification history
!::197507.05, rws, ver 00 tencol released for field use
!::198311.20, rws, ver 01 mten   released to field
!::198411.26, rws, ver 07 mten2  released to field
!::198506.10, rws, wrk    enhancements released to field
!::198507.22, rws, code   modified for mten3
!::198509.01, rws, ver 11 mten3  released to field
!::198708.10, rws, code   modified to use new mten4 gpn record format
!::199112.31, rws, ver 20 mten4 released to field
!::200001.13, rws, ver 21 mten4 released to field
!::200005.26, rws, code   restructured & documentation added             
!::200012.31, rws, ver 23 mten5 released                                 
!::200104.09, rws, code   added to calblin program                       
!::200208.09, rws, code   added subroutines gpnarc & gpnloa              
!********1*********2*********3*********4*********5*********6*********7*
!e::gpnhri
!  -------------------------------
!     m t e n  (version 3)
!              (version 4.22)
!              (version 5.23)
!  -------------------------------
!
      implicit real*8 (a-h,o-z)
!
      data tol0 /5.0d-15/
      data tol1 /5.0d-14/
      data tol2 /7.0d-03/
!
      twopi = 2.0d0*pi
!
!     test the longitude difference with tol1
!     tol1 is approximately 0.000000001 arc seconds
!
      ss = e2-e1
      if( dabs(ss).lt.tol1 ) then
        e2 = e2+tol1
!        write(*,*) ' longitudal difference is near zero '
!                 
        r2 = p2
        r1 = p1
        call gpnarc ( a, f, esq, pi, r1, r2, arc )
        s  = dabs( arc )
!
        if( p2.gt.p1 ) then
          az1 = 0.0d0
          az2 = pi
        else
          az1 = pi   
          az2 = 0.0d0
        endif
        return 
      endif
!
!     test for longitude over 180 degrees
!
      dlon = e2-e1
!
      if( dlon.ge.0.0d0 ) then
        if( pi.le.dlon .and. dlon.lt.twopi ) then
          dlon = dlon-twopi
        endif
      else
        ss = dabs(dlon)
        if( pi.le.ss .and. ss.lt.twopi ) then
          dlon = dlon+twopi
        endif
      endif
!
      ss = dabs( dlon )
      if( ss.gt.pi ) then
!::     write(*,*) '  '
!::     write(*,*) ' Longitude difference over 180 degrees  '  
!::     write(*,*) ' Turn it around '
        ss = twopi-ss
      endif
!
!     compute the limit in longitude (alimit), it is equal 
!     to twice the distance from the equator to the pole,
!     as measured along the equator (east/ewst)
!
      alimit = pi*(1.0d0-f)
!
!     test for anti-nodal difference      
!
      if( ss.ge.alimit ) then
        r1 = dabs(p1)
        r2 = dabs(p2)
!
!       latitudes r1 & r2 are not near the equator
!
        if( r1.gt.tol2 .and. r2.gt.tol2 ) then
          goto 60
        endif
!
!       longitude difference is greater than lift-off point
!       now check to see if  "both"  r1 & r2 are on equator
!
        if( r1.lt.tol1 .and. r2.gt.tol2 ) then
          goto 60
        endif
        if( r2.lt.tol1 .and. r1.gt.tol2 ) then
          goto 60
        endif
!
!       check for either r1 or r2 just off the equator but < tol2
!
        if (r1.gt.tol1 .or. r2.gt.tol1) then
          az1 = 0.0d0
          az2 = 0.0d0
          s   = 0.0d0
          return 
        endif
!
!       compute the azimuth to anti-nodal point
!
!::     write(*,*) '  '
!::     write(*,*) ' Longitude difference beyond lift-off point '  
!::     write(*,*) '  '
!
        call gpnloa (a,f,esq,pi,dlon,az1,az2,aa,bb,sms)
!
!       compute the equatorial distance & geodetic
!
        equ = a*dabs(dlon)
        s   = equ-sms
        return 
      endif
!
   60 continue
!
      f0   = (1.0d0-f)
      b    = a*f0
      epsq = esq/(1.0d0-esq)
      f2   = f*f     
      f3   = f*f2    
      f4   = f*f3    
!
!     the longitude difference 
!
      dlon  = e2-e1   
      ab    = dlon      
      kount = 0    
!
!     the reduced latitudes    
!
      u1    = f0*dsin(p1)/dcos(p1)     
      u2    = f0*dsin(p2)/dcos(p2)
!
      u1    = datan(u1)
      u2    = datan(u2)
!
      su1   = dsin(u1)    
      cu1   = dcos(u1)    
!
      su2   = dsin(u2)
      cu2   = dcos(u2)
!
!     counter for the iteration operation
!
    1 kount = kount+1     
!
      clon  = dcos(ab)   
      slon  = dsin(ab)   
!
      csig  = su1*su2+cu1*cu2*clon  
      ssig  = dsqrt((slon*cu2)**2+(su2*cu1-su1*cu2*clon)**2)  
!
      sig   = datan2(ssig,csig)
      sinalf=cu1*cu2*slon/ssig
!
      w   = (1.0d0-sinalf*sinalf)
      t4  = w*w   
      t6  = w*t4   
!
!     the coefficients of type a      
!
      ao  = f-f2*(1.0d0+f+f2)*w/4.0d0+3.0d0*f3*(1.0d0+                  &
     &        9.0d0*f/4.0d0)*t4/16.0d0-25.0d0*f4*t6/128.0d0
      a2  = f2*(1.0d0+f+f2)*w/4.0d0-f3*(1.0d0+9.0d0*f/4.0d0)*t4/4.0d0+  &
     &        75.0d0*f4*t6/256.0d0
      a4  = f3*(1.0d0+9.0d0*f/4.0d0)*t4/32.0d0-15.0d0*f4*t6/256.0d0
      a6  = 5.0d0*f4*t6/768.0d0
!
!     the multiple angle functions    
!
      qo  = 0.0d0
      if( w.gt.tol0 ) then
        qo = -2.0d0*su1*su2/w
      endif     
!
      q2  = csig+qo
      q4  = 2.0d0*q2*q2-1.0d0    
      q6  = q2*(4.0d0*q2*q2-3.0d0)      
      r2  = 2.0d0*ssig*csig      
      r3  = ssig*(3.0d0-4.0d0*ssig*ssig) 
!
!     the longitude difference 
!
      s   = sinalf*(ao*sig+a2*ssig*q2+a4*r2*q4+a6*r3*q6)    
      xz  = dlon+s   
!
      xy  = dabs(xz-ab)    
      ab  = dlon+s   
!
      if( xy.lt.0.5d-13 ) then
        goto 4
      endif
!
      if( kount.le.7 ) then
        goto 1
      endif
!
!     the coefficients of type b      
!
    4 z   = epsq*w
!
      bo  = 1.0d0+z*(1.0d0/4.0d0+z*(-3.0d0/64.0d0+z*(5.0d0/256.0d0-     &
     &         z*175.0d0/16384.0d0)))      
      b2  = z*(-1.0d0/4.0d0+z*(1.0d0/16.0d0+z*(-15.0d0/512.0d0+         &
     &         z*35.0d0/2048.0d0)))  
      b4  = z*z*(-1.0d0/128.0d0+z*(3.0d0/512.0d0-z*35.0d0/8192.0d0))
      b6  = z*z*z*(-1.0d0/1536.0d0+z*5.0d0/6144.0d0)    
!
!     the distance in meters   
!
      s   = b*(bo*sig+b2*ssig*q2+b4*r2*q4+b6*r3*q6) 
!
!     first compute the az1 & az2 for along the equator
!
      if( dlon.gt.pi ) then
        dlon = (dlon-2.0d0*pi)
      endif
!
      if( dabs(dlon).gt.pi ) then
        dlon = (dlon+2.0d0*pi)
      endif
!
      az1 = pi/2.0d0
      if( dlon.lt.0.0d0 ) then
        az1 = 3.0d0*az1
      endif
!
      az2 = az1+pi
      if( az2.gt.2.0d0*pi ) then
        az2 = az2-2.0d0*pi
      endif
!
!     now compute the az1 & az2 for latitudes not on the equator
!
      if( .not.(dabs(su1).lt.tol0 .and. dabs(su2).lt.tol0) ) then
        tana1 =  slon*cu2/(su2*cu1-clon*su1*cu2)  
        tana2 =  slon*cu1/(su1*cu2-clon*su2*cu1)  
        sina1 =  sinalf/cu1
        sina2 = -sinalf/cu2      
!
!       azimuths from north,longitudes positive east  
!
        az1   = datan2(sina1,sina1/tana1)   
        az2   = pi-datan2(sina2,sina2/tana2)
      endif
!
      if( az1.lt.0.0d0 ) then
        az1 = az1+2.0d0*pi   
      endif
!
      if( az2.lt.0.0d0 ) then
        az2 = az2+2.0d0*pi
      endif
!
      return     
      end 

!B::GPNLOA
!
      SUBROUTINE GPNLOA (AMAX,FLAT,ESQ,PI,DL,AZ1,AZ2,AO,BO,SMS)
!
!********1*********2*********3*********4*********5*********6*********7*
!
! NAME:        GPNLOA
! VERSION:     200005.26
! WRITTEN BY:  ROBERT (Sid) SAFFORD
! PURPOSE:     SUBROUTINE TO COMPUTE THE LIFF-OFF-AZIMUTH CONSTANTS
!
! INPUT PARAMETERS:
! -----------------
! AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
! FLAT         FLATTENING (0.0033528 ... )
! ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
! PI           3.14159...
! DL           LON DIFFERENCE
! AZ1          AZI AT STA 1 -> STA 2
!
! OUTPUT PARAMETERS:
! ------------------
! AZ2          AZ2 AT STA 2 -> STA 1
! AO           CONST
! BO           CONST
! SMS          DISTANCE ... EQUATORIAL - GEODESIC  (S - s)   "SMS"
!
! LOCAL VARIABLES AND CONSTANTS:
! ------------------------------
! GLOBAL VARIABLES AND CONSTANTS:
! -------------------------------
!
!    MODULE CALLED BY:    GENERAL 
!
!    THIS MODULE CALLS:   
!       LLIBFORE/ DSIN,   DCOS,   DABS,   DASIN 
!
!    INCLUDE FILES USED:
!    COMMON BLOCKS USED:  
!
!    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
!                MS-DOS Operating System
!    COMMENTS:
!********1*********2*********3*********4*********5*********6*********7*
!::MODIFICATION HISTORY
!::1985xx.xx, RWS, CODE   CREATED               
!::198506.10, RWS, WRK    ENHANCEMENTS RELEASED TO FIELD
!::198509.01, RWS, VER 11 MTEN3  RELEASED TO FIELD
!::198512.18, RWS, CODE   MODIFIED FOR MTEN3
!::198708.10, RWS, CODE   MODIFIED TO USE NEW MTEN4 GPN RECORD FORMAT
!::199112.31, RWS, VER 20 MTEN4 RELEASED TO FIELD
!::200001.13, RWS, VER 21 MTEN4 RELEASED TO FIELD
!::200005.26, RWS, CODE   RESTRUCTURED & DOCUMENTATION ADDED             
!::200012.31, RWS, VER 23 MTEN5 RELEASED                                 
!********1*********2*********3*********4*********5*********6*********7*
!E::GPNLOA
! ---------------------------
!     M T E N  (VERSION 3)
!              (VERSION 4.22)
!              (VERSION 5.23)
! ---------------------------
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DATA TT/5.0D-13/
!
      DLON = DABS(DL)
      CONS = (PI-DLON)/(PI*FLAT)
      F    = FLAT
!
!     COMPUTE AN APPROXIMATE AZ
!
      AZ   = DASIN(CONS)
!
      T1   =    1.0D0
      T2   =  (-1.0D0/4.0D0)*F*(1.0D0+F+F*F)
      T4   =    3.0D0/16.0D0*F*F*(1.0D0+(9.0D0/4.0D0)*F)
      T6   = (-25.0D0/128.0D0)*F*F*F
!
      ITER = 0
    1 ITER = ITER+1
      S    = DCOS(AZ)
      C2   = S*S
!
!     COMPUTE NEW AO
!
      AO   = T1 + T2*C2 + T4*C2*C2 + T6*C2*C2*C2
      CS   = CONS/AO
      S    = DASIN(CS)
      IF( DABS(S-AZ).LT.TT ) then
        GOTO 2
      ENDIF
!
      AZ   = S
      IF( ITER.LE.6 ) then
        GOTO 1
      ENDIF
!
    2 AZ1  = S
      IF( DL.LT.0.0D0 ) then
        AZ1 = 2.0D0*PI-AZ1
      ENDIF
!
      AZ2  = 2.0D0*PI-AZ1
!
!     EQUATORIAL - GEODESIC  (S - s)   "SMS"
!
      ESQP = ESQ/(1.0D0-ESQ)
      S    = DCOS(AZ1)
!
      U2   = ESQP*S*S
      U4   = U2*U2
      U6   = U4*U2
      U8   = U6*U2
!
      T1   =     1.0D0
      T2   =    (1.0D0/4.0D0)*U2
      T4   =   (-3.0D0/64.0D0)*U4
      T6   =    (5.0D0/256.0D0)*U6
      T8   = (-175.0D0/16384.0D0)*U8
!
      BO   = T1 + T2 + T4 + T6 + T8
      S    = DSIN(AZ1)
      SMS  = AMAX*PI*(1.0D0 - FLAT*DABS(S)*AO - BO*(1.0D0-FLAT))
!
      RETURN
      END
