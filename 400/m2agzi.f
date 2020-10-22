      SUBROUTINE m2agzi(MODE)
!     ***********************
!======================================================================!
!======================================================================!
!    Z boson exchange contribution to low angle Bhabha scattering,     !
!    both exponentiated and non-exponentiated up to O(alpha^1).        !
!======================================================================!
!======================================================================!
! IMPORTANT: The non-exponentiated version should be run with a photon !
!            removal option switched ON, i.e. KEYREM=0 !!!             !
!======================================================================!
! The following weights are implemented (Z-contributions only!):       !
! NOTE: For practical calculations of the Z exchange contribution      !
!       in the luminosity measurement we RECOMMEND to use the weight:  !
!       wtset(12)                                                      !
!       ^^^^^^^^^                                                      !
!======================================================================!
!    This routine fills wtset(10-29) and wtset(80-99)                  !
!----------------------------------------------------------------------!
!    1) STANDARD version with EXPONENTIATION                           !
!      wtset(10)  -  O(alpha^0)exp                                     !
!      wtset(11)  -  O(alpha^1)exp, no vacuum polarization             !
! ==>  wtset(12)  -  O(alpha^1)exp, with vacuum polarization    <===   !
!----------------------------------------------------------------------!
!    2) SPECIAL version with EXPONENTIATION                            !
!    with explicit Z self energy correction                            !
!    (for some tests, mainly comparisons with ALIBABA):                !
!      wtset(20)  -  O(alpha^0)exp                                     !
!      wtset(21)  -  O(alpha^1)exp, no vac. pol., no Z self en.        !
!      wtset(22)  -  O(alpha^1)exp, with vac. pol. and Z self en.      !
! Note: Version 2 gives generally the same results as version 1,       !
!       but it is kept for the sake of some tests.                     !
!======================================================================!
!    Non-exponentiated case:                                           !
!    3) STANDARD version                                               !
!      wtset(80)  -  O(alpha^0)                                        !
!      wtset(81)  -  O(alpha^1), no vacuum polarization                !
! -->  wtset(82)  -  O(alpha^1), with vacuum polarization       <---   !
!----------------------------------------------------------------------!
!    4) SPECIAL version                                                !
!    with O(alpha^1) ansatz by W.Beenakker     !
!    & B.Pietrzyk (for some tests, mainly comparisons with BABAMC):    !
!      wtset(90)  -  O(alpha^0)                                        !
!      wtset(91)  -  O(alpha^1), no vac. pol., no Z self en.           !
!      wtset(92)  -  O(alpha^1), with vac. pol. and Z self en.         !
!======================================================================!
!    In order to make the program more efficient we introduced         !
!    a number of keys which allow the user to choose a desired         !
!    version of the Z- contribution calculation from the list          !
!    above (any combination of them can be chosen as well).            !
!    Those keys are provided trough the integer array: KeyZex(4)       !
!    and have the meaning as follows:                                  !
!      KeyZex(1) =0/1   version no 1 (see above) switched OFF/ON       !
!      KeyZex(2) =0/1   version no 2 (see above) switched OFF/ON       !
!      KeyZex(3) =0/1   version no 3 (see above) switched OFF/ON       !
!      KeyZex(4) =0/1   version no 4 (see above) switched OFF/ON       !
!    The default set up is: DATA KeyZex /1,0,0,0/                      !
!======================================================================!
!    Written by:  Wieslaw Placzek          Knoxville, March 1995       !
!    Last update: 12.10.1995               by: W. Placzek              !
!**********************************************************************!
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0, ALFINV=137.03604D0)
      PARAMETER( ALFA= 1D0/ALFINV, ALF1 = 1D0/PI/ALFINV)
      COMMON / BHPAR1 / DEL,EPSCM,THMIN,XMIVIS
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR3 / KEYRAD,KEYOPT
      COMMON / BHPARZ / AMAZ,GAMMZ,SINW2,GV,GA
      COMMON / TRANSR / TRAN,TRMIN,TRMAX
      COMMON / MOMS1  / TRANP,P1(4),P2(4),PHOT1(100,4),PHSU1(4),NPHOT1
      COMMON / MOMS2  / TRANQ,Q1(4),Q2(4),PHOT2(100,4),PHSU2(4),NPHOT2
      COMMON / MOMZ1  / AL1(50),BE1(50),WTM1(50),MK1(50)
      COMMON / MOMZ2  / AL2(50),BE2(50),WTM2(50),MK2(50)
      COMMON / PSIPHI / TH1,EXT1,EXB1,PSI1,EXW1,
     $                  TH2,EXT2,EXB2,PSI2,EXW2,EXE2,FIF,PHI
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      SAVE   / BHPAR2 /, / BHPAR3 /, / TRANSR /
      SAVE   / MOMS1  /, / MOMS2  /, / MOMZ1  /, / MOMZ2  /
      SAVE   / WGTALL /, / BHPARZ /, / PSIPHI /
      DIMENSION RP(4),RQ(4)
      DIMENSION KeyZex(4)
! Array of keys for Z-contributions switch on/off and a default set-up
      DATA KeyZex /1,0,1,0/
      SAVE KeyZex
! --------------------------------------------------------------------
      DATA icont /0/
      SAVE icont
! ------------------ Inline functions ------------------
! Dot products for some Lorentz invariants
      sdot(x1,y1,z1,e1,x2,y2,z2,e2)
     $            = (e1+e2)**2-(x1+x2)**2-(y1+y2)**2-(z1+z2)**2
      tdot(x1,y1,z1,e1,x2,y2,z2,e2)
     $            = (e1-e2)**2-(x1-x2)**2-(y1-y2)**2-(z1-z2)**2
! --------------------------------------------------------------------
      icont=icont+1
! .........................................
      DO k=1,4
        RP(k)= P2(K)+PHSU1(K)
        RQ(k)= Q2(K)+PHSU2(K)
      ENDDO
! Some invariants
      tr = -TRAN
      t  = -TRANP
      t1 = -TRANQ
      s  = CMSEne**2
      s1 = sdot(P2(1),P2(2),P2(3),P2(4) ,Q2(1),Q2(2),Q2(3),Q2(4))
      u  = tdot(P1(1),P1(2),P1(3),P1(4) ,Q2(1),Q2(2),Q2(3),Q2(4))
      u1 = tdot(P2(1),P2(2),P2(3),P2(4) ,Q1(1),Q1(2),Q1(3),Q1(4))
      xpp  = P2(4)*P1(4)-P2(3)*P1(3)-P2(2)*P1(2)-P2(1)*P1(1)
      xpr  = P1(4)*RP(4)-P1(3)*RP(3)-P1(2)*RP(2)-P1(1)*RP(1)
      xqq  = Q2(4)*Q1(4)-Q2(3)*Q1(3)-Q2(2)*Q1(2)-Q2(1)*Q1(1)
      xqr  = Q1(4)*RQ(4)-Q1(3)*RQ(3)-Q1(2)*RQ(2)-Q1(1)*RQ(1)
!.....Variables for O(alpha1) virtual+soft corrections
      wdz = 0.5*dlog(TRAN/AMEL**2)
      vdz = 0.5*dlog(s/AMEL**2)
      vdz1= 0.5*dlog(s1/AMEL**2)
      udz = 0.5*dlog(abs(4*AMEL**2-s -tr)/AMEL**2)
      udz1= 0.5*dlog(abs(4*AMEL**2-s1-tr)/AMEL**2)
!*****************************************************
!  Crude MC distribution (S-tilde factors omitted)
      crude  =  s**2/tran**2
      deltp  =  AMEL**2/TRANP
      deltq  =  AMEL**2/TRANQ
      gamp   =  2*alf1 *( dlog(1/deltp)-1 )
      gamq   =  2*alf1 *( dlog(1/deltq)-1 )
!==================================================================!
      b10uza =0
      b10uzv =0
      b10lza =0
      b10lzv =0
      s10uza =0
      s10uzv =0
      s10lza =0
      s10lzv =0
! Vacuum polarization factors
      vpft  = 1/(1+RePi(-TRAN))       ! <===
      vpfac = 1/(1+RePi(  s )) *vpft
      vpfac1= 1/(1+RePi(  s1)) *vpft

!        *************************************************
!        *************************************************
!        **              EXPONENTIATION                 **
!        *************************************************
!        *************************************************
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!>>> If exp. version is chosen
      IF (KeyZex(1) .EQ. 1 .OR. KeyZex(2) .EQ. 1) THEN
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!... Z mass
      amZ  = AMAZ
      amZ2 = amZ**2
!... s-shift of the Z-width
      GamZ = GAMMZ *s /amZ2
      GamZ1= GAMMZ *s1/amZ2
!... Z coupling const.
      Zcc = gv**2 + ga**2
!... kinematical factors for Z contribution
      fasc  = Zcc *(t+t1)*(u**2+u1**2)/(4*t*t1)
      facZZ = Zcc**2 *(u**2+u1**2+t**2+t1**2)/4
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!>>> If chosen version no 1
      IF (KeyZex(1) .EQ. 1) THEN
*--------------------------------------------------------------------
*--------------------------------------------------------------------
      prki = fasc/( (s -amZ2)**2 + (amZ*gamZ )**2 )
      prki1= fasc/( (s1-amZ2)**2 + (amZ*gamZ1)**2 )
!... Real and imaginary contributions at the lowest order
      delZR = (s -amZ2) *prki
      delZI =-gamZ *amZ *prki
!... with s-shift (s --> s1)
      delZR1= (s1-amz2) *prki1
      delZI1=-gamZ1*amZ *prki1
!... with vacuum polar. corr.
      devZR = delZR *vpfac   ! <===
      devZI = delZI *vpfac   ! <===
      devZR1= delZR1*vpfac1  ! <===
      devZI1= delZI1*vpfac1  ! <===
!... Z-gamma interference, soft + virt. corrections to beta_0
      psis = dacos( (amZ2-s )/sqrt((s -amZ2)**2+(amZ*GamZ )**2))
      psis1= dacos( (amZ2-s1)/sqrt((s1-amZ2)**2+(amZ*GamZ1)**2))
      delRe =2*(wdz-udz)*dlog(s/(4*(sqrt(s)-amZ)**2+GamZ**2))
     $      +2*(vdz**2-udz**2 ) +4*wdz*(udz -vdz) +vdz +wdz
     $      +2*udz  +dilogy(1-TRAN/s) -dilogy(TRAN/s) -5*pi**2/6 -2
      delIm = pi*( 2*(wdz-vdz) + 4*(udz -wdz)*psis/pi + 1.5 )
      delRe1=2*(wdz-udz1)*dlog(s1/(4*(sqrt(s1)-amZ)**2+gamZ1**2))
     $      +2*(vdz1**2-udz1**2) +4*wdz*(udz1-vdz1)+ vdz1 +wdz
     $      +2*udz1 +dilogy(1-TRAN/s1) -dilogy(TRAN/s1) -5*pi**2/6 -2
      delIm1=pi*( 2*(wdz-vdz1) + 4*(udz1-wdz)*psis1/pi + 1.5 )
      gZRe = 1 +alf1*delRe
      gZIm = alf1*delIm
      gZRe1= 1 +alf1*delRe1
      gZIm1= alf1*delIm1
!=================================================================
!=================================================================
!                           beta0
!=================================================================
!=================================================================
!... no vacuum pol. cor.
      delZa = delZR *gZRe  + delZI *gZIm
      delZa1= delZR1*gZRe1 + delZI1*gZIm1
!... with vacuum pol. cor.
      delZv = devZR *gZRe  + devZI *gZIm   ! <===
      delZv1= devZR1*gZRe1 + devZI1*gZIm1  ! <===
!... Born Z(s)-Z(s) contribution, no vac. pol. corr.
      deZZa = facZZ/( (s -amZ2)**2 + (amZ*gamZ )**2 )
      deZZa1= facZZ/( (s1-amZ2)**2 + (amZ*gamZ1)**2 )
!... with vac. pol. corr.
      deZZv = deZZa *(vpfac /vpft)**2  ! <===
      deZZv1= deZZa1*(vpfac1/vpft)**2  ! <===
!-----
      bt01Za = (delZa+delZa1)/2 +(deZZa+deZZa1)/2
      bt01Zv = (delZv+delZv1)/2 +(deZZv+deZZv1)/2  ! <===
!-----
!... some factors for calculation of hard photon contrib.
      Zpro = (s -amZ2)**2 +(amZ*GamZ )**2
      Zpro1= (s1-amZ2)**2 +(amZ*GamZ1)**2
      faRe = Zcc*(s -amZ2)*(u**2+u1**2)
      faRe1= zcc*(s1-amz2)*(u**2+u1**2)
      faIm = 2*gv*ga *amZ*gamZ *(u**2-u1**2)*s
      faIm1= 2*gv*ga *amZ*gamZ1*(u**2-u1**2)*s
      ENDIF
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!>>> If chosen version no 2
      IF (KeyZex(2) .EQ. 1) THEN
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!... the same as above but with the explicit Z self energy correction
      seki = fasc/( (s -amZ2)**2 + (amZ*GAMMZ)**2 )
      seki1= fasc/( (s1-amZ2)**2 + (amZ*GAMMZ)**2 )
!... Real and imaginary contributions at the lowest order
      selZR = (s -amZ2) *seki
      selZI =-GAMMZ*amZ *seki
      selZR1= (s1-amZ2) *seki1
      selZI1=-GAMMZ*amZ *seki1
!... Z self energy corr.; from W. Beenakker & B. Pietrzyk, CERN-TH.6649/92
      Zg = GAMMZ/amZ
      Zy = (s -amZ2)/(amZ*GAMMZ )
      Zy1= (s1-amZ2)/(amZ*GAMMZ)
      PiZ= RePi(amZ2)
      Zsfac = ( 1 - 2*Zy *Zg/(1+Zy**2 ) )/(1+PiZ)
      Zsfac1= ( 1 - 2*Zy1*Zg/(1+Zy1**2) )/(1+PiZ)
!... with vac. pol. and Z self en. corr.
      sevZR = selZR *vpft*Zsfac
      sevZI = selZI *vpft*Zsfac
      sevZR1= selZR1*vpft*Zsfac1
      sevZI1= selZI1*vpft*Zsfac1
!... Z-gamma interference, Born + virt. corrections to beta_0
      psise = dacos( (amZ2-s )/sqrt((s -amZ2)**2+(amZ*GAMMZ)**2))
      psise1= dacos( (amZ2-s1)/sqrt((s1-amZ2)**2+(amZ*GAMMZ)**2))
      selRe =2*(wdz-udz)*dlog(s/(4*(sqrt(s)-amZ)**2+GAMMZ**2))
     $      +2*(vdz**2-udz**2 ) +4*wdz*(udz -vdz) +vdz +wdz
     $      +2*udz  +dilogy(1-TRAN/s) -dilogy(TRAN/s) -5*pi**2/6 -2
      selIm = pi*( 2*(wdz-vdz) + 4*(udz -wdz)*psise/pi + 1.5 )
      selRe1=2*(wdz-udz1)*dlog(s1/(4*(sqrt(s1)-amZ)**2+GAMMZ**2))
     $      +2*(vdz1**2-udz1**2) +4*wdz*(udz1-vdz1)+ vdz1 +wdz
     $      +2*udz1 +dilogy(1-TRAN/s1) -dilogy(TRAN/s1) -5*pi**2/6 -2
      selIm1=pi*( 2*(wdz-vdz1) + 4*(udz1-wdz)*psise1/pi + 1.5 )
      sZRe = 1 +alf1*selRe
      sZIm = alf1*selIm
      sZRe1= 1 +alf1*selRe1
      sZIm1= alf1*selIm1
!=================================================================
!=================================================================
!                           beta0
!=================================================================
!=================================================================
!... no vacuum pol., no Z self en.
      selZa = selZR *sZRe  + selZI *sZIm
      selZa1= selZR1*sZRe1 + selZI1*sZIm1
C... with vacuum pol. and Z self en.
      selZv = sevZR *sZRe  + sevZI *sZIm
      selZv1= sevZR1*sZRe1 + sevZI1*sZIm1
C... Born Z(s)-Z(s) contribution, no Z-self energy corr.
      seZZa = facZZ/( (s -amZ2)**2 + (amZ*GAMMZ)**2 )
      seZZa1= facZZ/( (s1-amZ2)**2 + (amZ*GAMMZ)**2 )
C... with Z-self energy corr.
      seZZv = seZZa *Zsfac**2
      seZZv1= seZZa1*Zsfac1**2
!-----
      st01Za = (selZa+selZa1)/2 +(seZZa+seZZa1)/2
      st01Zv = (selZv+selZv1)/2 +(seZZv+seZZv1)/2
!-----
!... some factors for calculation of hard photon contrib.
      spro = (s -amZ2)**2 +(amZ*GAMMZ)**2
      spro1= (s1-amZ2)**2 +(amZ*GAMMZ)**2
      saRe = Zcc*(s -amZ2)*(u**2+u1**2)
      saRe1= Zcc*(s1-amZ2)*(u**2+u1**2)
      saIm = 2*gv*ga *amZ*GAMMZ*(u**2-u1**2)*s
      ENDIF

!=================================================================
!=================================================================
!              Contributions from beta1 UPPER line
!=================================================================
!=================================================================
      DO i=1,nphot1
! numerically safe variables
        a = al1(i)
        b = be1(i)
        y = a +b*deltp
        z = b +a*deltp
! one photon bremss. distribution
        p1k = xpr*z
        p2k = xpr*y
        q1k=Q1(4)*PHOT1(i,4)
        q2k=Q2(4)*PHOT1(i,4)
        DO ii=1,3
          q1k=q1k -Q1(ii)*PHOT1(i,ii)
          q2k=q2k -Q2(ii)*PHOT1(i,ii)
        ENDDO
! soft factor
        sfcc = 2*xpp/(p1k*p2k)
        epsp2k=P2(1)*PHOT1(i,2)-P2(2)*PHOT1(i,1)
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!>>> If chosen version no 1
      IF (KeyZex(1) .EQ. 1) THEN
*--------------------------------------------------------------------
*--------------------------------------------------------------------
! no s-shift
        au0=( faRe *(sfcc +s1/(p2k*q2k) +u/(p1k*q2k) +4/p2k)
     $       +faIm *epsp2k/(p1k*p2k*q2k) )/Zpro /t1/4
! with s-shift
        au1=( faRe1*(sfcc +s/(p1k*q1k) +u1/(p2k*q1k) -4/p1k)
     $       -faIm1*epsp2k/(p1k*p2k*q1k) )/Zpro1/t1/4
! no vac. pol.
        au1a=au0 +au1
! with vac. pol.
        au1v= au0*vpfac + au1*vpfac1          ! <===
! beta1 O(alf1)
        bu10Za = au1a/sfcc -(delZR+delZR1)/2
        bu10Zv = au1v/sfcc -(devZR+devZR1)/2  ! <===
        b10uZa = b10uZa +bu10Za
        b10UZv = b10uZv +bu10Zv  ! <===
      ENDIF
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!>>> If chosen version no 2
      IF (KeyZex(2) .EQ. 1) THEN
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!... as above but for Z self energy case
! no s-shift
        su0=( saRe *(sfcc +s1/(p2k*q2k) +u/(p1k*q2k) +4/p2k)
     $       +saIm *epsp2k/(p1k*p2k*q2k) )/spro /t1/4
! with s-shift
        su1=( saRe1*(sfcc +s/(p1k*q1k) +u1/(p2k*q1k) -4/p1k)
     $       -saIm *epsp2k/(p1k*p2k*q1k) )/spro1/t1/4
! no vac. pol., no Z self en.
        su1a=su0 +su1
! with vac. pol. and Z self en.
        su1v=(su0*Zsfac + su1*Zsfac1) *vpft
! beta1 O(alf1)
        su10Za = su1a/sfcc -(selZR+selZR1)/2
        su10Zv = su1v/sfcc -(sevZR+sevZR1)/2
        s10uZa = s10uZa  +su10Za
        s10uZv = s10uZv  +su10Zv
      ENDIF
*--------------------------------------------------------------------
*--------------------------------------------------------------------
      ENDDO
!=================================================================
!=================================================================
!             Contributions from beta1 LOWER line
!=================================================================
!=================================================================
      DO i=1,nphot2
! numerically safe variables
        a = al2(i)
        b = be2(i)
        y = a+b*deltq
        z = b+a*deltq
! one photon bremss. distribution
        q1k = xqr*z
        q2k = xqr*y
        p1k=P1(4)*PHOT2(i,4)
        p2k=P2(4)*PHOT2(i,4)
        DO ii=1,3
          p1k=p1k-P1(ii)*PHOT2(i,ii)
          p2k=p2k-P2(ii)*PHOT2(i,ii)
        ENDDO
! soft factor
        sfcc = 2*xqq/(q1k*q2k)
        epsp2k=P2(1)*PHOT2(i,2)-P2(2)*PHOT2(i,1)
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!>>> If chosen version no 1
      IF (KeyZex(1) .EQ. 1) THEN
*--------------------------------------------------------------------
*--------------------------------------------------------------------
! no s-shift
        ad0=( faRe *(sfcc +s1/(p2k*q2k) +u1/(p2k*q1k) +4/q2k)
     $       +faIm *epsp2k/(p2k*q1k*q2k) )/Zpro /t/4
! with s-shift
        ad1=( faRe1*(sfcc +s/(p1k*q1k)  +u /(p1k*q2k) -4/q1k)
     $       -faIm1*epsp2k/(p1k*q1k*q2k) )/Zpro1/t/4
! no vac. pol.
        ad1a=ad0 +ad1
! with vac. pol.
        ad1v=ad0*vpfac + ad1*vpfac1           ! <===
! beta1 O(alf1)
        bl10Za = ad1a/sfcc -(delZR+delZR1)/2
        bl10Zv = ad1v/sfcc -(devZR+devZR1)/2  ! <===
        b10lZa = b10lZa  +bl10Za
        b10lZv = b10lZv  +bl10Zv  ! <===
      ENDIF
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!>>> If chosen version no 2
      IF (KeyZex(2) .EQ. 1) THEN
*--------------------------------------------------------------------
*--------------------------------------------------------------------
!... as above but for Z self energy case
! no s-shift
        sd0=( saRe *(sfcc +s1/(p2k*q2k) +u1/(p2k*q1k) +4/q2k)
     $       +saIm *epsp2k/(p2k*q1k*q2k) )/spro /t/4
! with s-shift
        sd1=( saRe1*(sfcc +s/(p1k*q1k)  +u /(p1k*q2k) -4/q1k)
     $       -saIm *epsp2k/(p1k*q1k*q2k) )/spro1/t/4
! no vac. pol., no Z self en.
        sd1a=sd0 +sd1
! with vac. pol.
        sd1v=(sd0*Zsfac + sd1*Zsfac1) *vpft
! beta1 O(alf1)
        sl10Za = sd1a/sfcc -(selZR+selZR1)/2
        sl10Zv = sd1v/sfcc -(sevZR+sevZR1)/2
        s10lZa = s10lZa  +sl10Za
        s10lZv = s10lZv  +sl10Zv
      ENDIF
*--------------------------------------------------------------------
*--------------------------------------------------------------------
      ENDDO

!        *************************************************
!        *************************************************
!        **              EXPONENTIATION                 **
!        **          Definitions of MC weights          **
!        *************************************************
!        *************************************************
!              ---------------------------------
!              !  UPPER + LOWER line,    (A)   !
!              ---------------------------------
! zero order
      wtset( 10) = ( (delZR+delZR1)/2 +(deZZa+deZZa1)/2 )/crude
! first order - no vac. pol.
      wtset( 11) = (bt01Za + b10uZa + b10lZa) /crude
! first order - with vac. pol.                        ! <=== RECOMMENDED !!!
      wtset( 12) = (bt01Zv + b10uZv + b10lZv) /crude  ! <=== RECOMMENDED !!!
!... as above but with explicit Z self energy corr.
!    (for some tests, e.g. comparisons with ALIBABA)
! zero order
      wtset( 20) = ( (selZR+selZR1)/2 +(seZZa+seZZa1)/2 )/crude
! first order - no vac. pol., no Z self en.
      wtset( 21) = (st01Za + s10uZa + s10lZa) /crude
! first order - with vac. pol. and Z self energy
      wtset( 22) = (st01Zv + s10uZv + s10lZv) /crude
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      ENDIF
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



!==================================================================
!==================================================================
!                  NO-EXPONENTIATION version                      !
!==================================================================
!==================================================================

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!>>> If non-exp. version is chosen
      IF (KeyZex(3) .EQ. 1 .OR. KeyZex(4) .EQ. 1) THEN
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      delZB  =0
      delZ1a =0
      delZ1v =0
      selZ0  =0
      selZ1a =0
      selZ1v =0
      sevac  =0
      seZse  =0
!... Z mass
      amZ  = AMAZ
      amZ2 = amZ**2
!... s-shift of the Z-width
      GamZ = GAMMZ *s /amZ2
      GamZ1= GAMMZ *s1/amZ2
!... Z-width without shift
      GZex = GAMMZ
!... Z-width as in original BABAMC
      GZba= 2.3098
!... sin^2(theta_W) the same as in BABAMC
      sw2 = 0.2273
      ga1 = -1/(4*sqrt(sw2*(1-sw2)))
      gv1 = ga1*(1-4*sw2)
      Zcc1 = gv1**2 + ga1**2
!
      IF(nphot1 .EQ. 0.and.nphot2 .EQ. 0) THEN
!===================================================================
!===================================================================
!                        [00] No photons
!===================================================================
!===================================================================
!------
! O(alf0,1) distribution
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 3
        IF (KeyZex(3) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
        fki = Zcc *u**2/t
        fkZ = Zcc**2 *(u**2+t**2)/2
        Zpro =(s-amZ2)**2 + (GamZ*amZ)**2
! Born gamma(t)-Z(s) contribution
        dgZts = fki*(s-amZ2)/Zpro
! Born Z(s)-Z(s) contribution
        dZZss =fkZ/Zpro
! Lowest order contrib. with Z-width shift
        delZ0 = dgZts + dZZss
! Total Born Z contribution
        delZB = delZ0*Zpro/( (s-amZ2)**2 + (GAMMZ*amZ)**2 )
        ENDIF
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 4
        IF (KeyZex(4) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!... for comparisons with BABAMC
        ski = Zcc1 *u**2/t
        skZ = Zcc1**2 *(u**2+t**2)/2
        sproe =(s-amZ2)**2 + (GZex*amZ)**2
        sprob =(s-amZ2)**2 + (GZba*amZ)**2
! Born gamma(t)-Z(s) contribution
        sgZtse = ski*(s-amZ2)/sproe
        sgZtsb = ski*(s-amZ2)/sprob
! Born Z(s)-Z(s) contribution
        sZZsse=skZ/sproe
        sZZssb=skZ/sprob
! Total Born Z contribution
        selZ0 = sgZtse + sZZsse
! for the original BABAMC Z-width
        selZ0b= sgZtsb+ sZZssb
        ENDIF
*-------------------------------------------------
	Eb  = CMSEne/2
	xkm = epsCM*Eb             ! soft foton cut-off in CMS
!####################################################################
! case for CMSEne=amZ
	IF (abs(CMSEne-amZ).lt.1d-10) THEN
!####################################################################
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 3
          IF (KeyZex(3) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
          gZRe = 2*datan(2*xkm/GamZ)*(1+2*udz-2*vdz-2*wdz)
          gZIm = -pi*(2*(udz-vdz)+1.5)
! virtual+soft corretion, no vac. pol. corr.
          delZa = delZ0 + alf1*fki/amZ/GamZ*(gZRe+gZIm)
! with vac. pol. corr.
          delZv =( delZa + dZZss*(vpfac/vpft**2 -1) )*vpfac ! <---
          delZ1a =delZa
          delZ1v =delZv ! <---
          ENDIF
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 4
          IF (KeyZex(4) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
C... for comparisons with BABAMC
          sZRe = 2*datan(2*xkm/GZex)*(1+2*udz-2*vdz-2*wdz)
          sZIm = -pi*(2*(udz-vdz)+1.5)
! virtual+soft corretion, no Z self en. corr.
          selZa = selZ0 + alf1*ski/amZ/GZex*(sZRe+sZIm)
! Z self en. corr
          Zsefa = RePi(amZ2)
! virtual+soft corretion, with Z self en. corr.
          selZv = selZa + 2*sZZssb*Zsefa + selZ0b-selZ0
          selZ1a =selZa
          selZ1v =selZv
          ENDIF
!####################################################################
! otherwise (i.e. away from the Z peak)
        ELSE
!####################################################################
          hit = -t/s
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 3
          IF (KeyZex(3) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
          R0  = (CMSEne    -amZ)**2 +GamZ**2/4
          Rkm = (CMSEne-xkm-amZ)**2 +GamZ**2/4
          fiZ = datan( (-2*(CMSEne-amZ)+2*xkm)/GamZ )
     $         -datan(  -2*(CMSEne-amZ)       /GamZ )
          psis= dacos( (amZ2-s)/sqrt((s-amZ2)**2+(amZ*GamZ)**2) )
          gZRi = -2*( 1+2*(udz-vdz-wdz) )*( dlog(xkm/Eb)
     $                                     -gamZ/2/(CMSEne-amZ)*fiZ )
     $           -(1-2*vdz)*dlog((xkm/Eb)**2*R0/Rkm)
     $           -2*(udz-wdz)*dlog(xkm**2/Rkm)
     $           -2*udz**2 +2*vdz**2 +4*wdz*(udz-vdz) +3*(vdz+wdz)
     $           +dilogy(1-hit) -dilogy(hit) +pi**2/6 -4
     $           -GamZ*amZ/(s-amz2)*pi*(2*(wdz-vdz)+1.5
     $                                +4*(udz-wdz)*psis/pi)
! virtual+soft corretion, bo vac. pol. corr.
          delZa = delZ0 + alf1*dgZts*gZRi
!  with vac. pol. corr.
          delZv =( delZa +dZZss*(vpfac/vpft**2 -1) )*vpfac  ! <---
          delZ1a =delZa
          delZ1v =delZv    ! <---
          ENDIF
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 4
          IF (KeyZex(4) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!... for comparisons with BABAMC
          sR0  = (CMSene    -amZ)**2 +gZex**2/4
          sRkm = (CMSEne-xkm-amZ)**2 +gZex**2/4
          sfis = datan( (-2*(CMSEne-amZ)+2*xkm)/gZex )
     $          -datan(  -2*(CMSEne-amZ)       /gZex )
          spsis= dacos( (amZ2-s)/sqrt((s-amZ2)**2+(amZ*gZex)**2) )
          sZRi = -2*( 1+2*(udz-vdz-wdz) )*( dlog(xkm/Eb)
     $                                     -gZex/2/(CMSEne-amZ)*sfis )
     $           -(1-2*vdz)*dlog((xkm/Eb)**2*sR0/sRkm)
     $           -2*(udz-wdz)*dlog(xkm**2/sRkm)
     $           -2*udz**2 +2*vdz**2 +4*wdz*(udz-vdz) +3*(vdz+wdz)
     $           +dilogy(1-hit) -dilogy(hit) +pi**2/6 -4
     $           -gZex*amZ/(s-amz2)*pi*(2*(wdz-vdz)+1.5
     $                                +4*(udz-wdz)*spsis/pi)
! virtual+soft corretion
          selZa = selZ0 + alf1*sgZtse*sZRi
! vacuum polarization correction
          sevac =-sgZtse *RePi(-TRAN)
! Z self energy corr.; taken from W. Beenakker & B. Pietrzyk, CERN-TH.6649/92
          Zy = (s-amZ2)/(amz*gZba)
          Zg = gZba/amZ
          Zsefa = ( RePi(amZ2)*(1-Zy**2)-2*Zy*Zg )/(1+Zy**2)
          seZse = sgZtsb*Zsefa + 2*sZZssb*Zsefa + selZ0b-selZ0
          selZv = selZa + sevac + seZse
          selZ1a =selZa
          selZ1v =selZv
!--------------------------------------------------------------------
!--------------------------------------------------------------------
          ENDIF
!####################################################################
        ENDIF
!===================================================================
!===================================================================
!              [10] One upper line photon
!===================================================================
!===================================================================
      ELSEIF(nphot1 .EQ. 1.and.nphot2 .EQ. 0) THEN
! numerically safe variables
        a = al1(1)
        b = be1(1)
        y = a +b*deltp
        z = b +a*deltp
! one photon bremss. distribution
        p1k = xpr*z
        p2k = xpr*y
        q1k=Q1(4)*PHOT1(1,4)
        q2k=Q2(4)*PHOT1(1,4)
        DO ii=1,3
          q1k=Q1K-Q1(ii)*PHOT1(1,ii)
          q2k=Q2K-Q2(ii)*PHOT1(1,ii)
        ENDDO
! soft factor
        sfcc = 2*xpp/(p1k*p2k)
        sfc  = sfcc *a*b/(y*z)
        epsp2k=P2(1)*PHOT1(1,2)-P2(2)*PHOT1(1,1)
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 3
        IF (KeyZex(3) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
        Zpro = (s -amZ2)**2 +(amZ*GamZ )**2
        Zpro1= (s1-amZ2)**2 +(amZ*GamZ1)**2
        faRe = Zcc*(s -amZ2)*(u**2+u1**2)
        faRe1= Zcc*(s1-amZ2)*(u**2+u1**2)
        faIm = 2*gv1*ga1*amZ*GamZ *(u**2-u1**2)*s
        faIm1= 2*gv1*ga1*amZ*GamZ1*(u**2-u1**2)*s
        fau  = 4*u**2 /(u**2+u1**2)
        fau1 = 4*u1**2/(u**2+u1**2)
! no s-shift
        au0=( faRe *(sfcc +s1/(p2k*q2k) +u/(p1k*q2k) +4/p2k
     $       -fau *(amel/p2k)**2)
     $       +faIm *epsp2k/(p1k*p2k*q2k) )/Zpro /t1/4
! with s-shift
        au1=( faRe1*(sfcc +s/(p1k*q1k) +u1/(p2k*q1k) -4/p1k
     $       -fau1*(amel/p1k)**2)
     $       -faIm1*epsp2k/(p1k*p2k*q1k) )/Zpro1/t1/4
        au1a = au0 +au1
! with vacuum polarization
        au1v=au0*vpfac +au1*vpfac1  ! <---
        delZa = au1a/sfc
        delZv = au1v/sfc
        delZ1a= delZa
        delZ1v= delZv               ! <---
        ENDIF
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 4
        IF (KeyZex(4) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
C... for comparisons with BABAMC
        spro = (s -amZ2)**2 +(amZ*GZex)**2
        spro1= (s1-amZ2)**2 +(amZ*GZex)**2
        saRe = Zcc1*(s -amZ2)*(u**2+u1**2)
        saRe1= Zcc1*(s1-amZ2)*(u**2+u1**2)
        saIm = 2*gv1*ga1*amZ*GZex *(u**2-u1**2)*s
        sau  = 4*u**2 /(u**2+u1**2)
        sau1 = 4*u1**2/(u**2+u1**2)
! no s-shift
        su0=( saRe *(sfcc +s1/(p2k*q2k) +u/(p1k*q2k) +4/p2k
     $       -sau *(amel/p2k)**2)
     $       +saIm*epsp2k/(p1k*p2k*q2k) )/spro /t1/4
! with s-shift
        su1=( saRe1*(sfcc +s/(p1k*q1k) +u1/(p2k*q1k) -4/p1k
     $       -sau1*(amel/p1k)**2)
     $       -saIm *epsp2k/(p1k*p2k*q1k) )/spro1/t1/4
        su1a = su0 +su1
        su1v = su1a
        selZa = su1a/sfc
        selZv = su1v/sfc
        selZ1a= selZa
        selZ1v= selZv
        ENDIF
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      ELSEIF(nphot1 .EQ. 0.and.nphot2 .EQ. 1) THEN
!===================================================================
!===================================================================
!                [01] One lower line photon
!===================================================================
!===================================================================
! numerically safe variables
        a = al2(1)
        b = be2(1)
        y = a+b*deltq
        z = b+a*deltq
! one photon bremss. distribution
        q1k = xqr*z
        q2k = xqr*y
        p1k=P1(4)*PHOT2(1,4)
        p2k=P2(4)*PHOT2(1,4)
        DO ii=1,3
          p1k=p1k-P1(ii)*PHOT2(1,ii)
          p2k=p2k-P2(ii)*PHOT2(1,ii)
        ENDDO
! soft factor
        sfcc = 2*xqq/(q1k*q2k)
        sfc  = sfcc *a*b/(y*z)
        epsp2k=P2(1)*PHOT2(1,2)-P2(2)*PHOT2(1,1)
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 3
        IF (KeyZex(3) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
        Zpro = (s -amZ2)**2 +(amZ*GamZ )**2
        Zpro1= (s1-amZ2)**2 +(amZ*GamZ1)**2
        faRe = Zcc*(s -amZ2)*(u**2+u1**2)
        faRe1= Zcc*(s1-amZ2)*(u**2+u1**2)
        faIm = 2*gv1*ga1*amZ*GamZ *(u**2-u1**2)*s
        faIm1= 2*gv1*ga1*amZ*GamZ1*(u**2-u1**2)*s
        fau  = 4*u**2 /(u**2+u1**2)
        fau1 = 4*u1**2/(u**2+u1**2)
! no s-shift
        ad0=( faRe *(sfcc +s1/(p2k*q2k) +u1/(p2k*q1k) +4/q2k
     $       -fau1*(amel/q2k)**2)
     $       +faIm *epsp2k/(p2k*q1k*q2k) )/Zpro /t/4
! with s-shift
        ad1=( faRe1*(sfcc +s/(p1k*q1k) +u /(p1k*q2k) -4/q1k
     $       -fau *(amel/q1k)**2)
     $       -faIm1*epsp2k/(p1k*q1k*q2k) )/Zpro1/t/4
        ad1a = ad0 +ad1
! with vacuum polarization
        ad1v=ad0*vpfac +ad1*vpfac1   ! <---
        delZa = ad1a/sfc
        delZv = ad1v/sfc
        delZ1a= delZa
        delZ1v= delZv                ! <---
        ENDIF
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>>> If chosen version no 4
        IF (KeyZex(4) .EQ. 1) THEN
!--------------------------------------------------------------------
!--------------------------------------------------------------------
C... for comparisons with BABAMC
        spro = (s -amZ2)**2 +(amZ*GZex )**2
        spro1= (s1-amZ2)**2 +(amZ*GZex)**2
        saRe = Zcc1*(s -amZ2)*(u**2+u1**2)
        saRe1= Zcc1*(s1-amZ2)*(u**2+u1**2)
        saIm = 2*gv1*ga1*amZ*GZex*(u**2-u1**2)*s
        sau  = 4*u**2 /(u**2+u1**2)
        sau1 = 4*u1**2/(u**2+u1**2)
! no s-shift
        sd0=( saRe *(sfcc +s1/(p2k*q2k) +u1/(p2k*q1k) +4/q2k
     $       -sau1*(amel/q2k)**2)
     $       +saIm *epsp2k/(p2k*q1k*q2k) )/spro /t/4
! with s-shift
        sd1=( saRe1*(sfcc +s/(p1k*q1k) +u /(p1k*q2k) -4/q1k
     $       -sau *(amel/q1k)**2)
     $       -saIm *epsp2k/(p1k*q1k*q2k) )/spro1/t/4
        sd1a = sd0 +sd1
        sd1v = sd1a
        selZa = sd1a/sfc
        selZv = sd1v/sfc
        selZ1a= selZa
        selZ1v= selZv
        ENDIF
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      ENDIF
!===================================================================
!===================================================================
!
!        ***************************************************
!        ****             NO EXPONENTIATON              ****
!        ****         Definitions of MC weights         ****
!        ***************************************************
!          ---------------------------------
!          /////  UPPER + LOWER line   /////
!          ---------------------------------
! YFS form-factor (to be divided off)
      Emin = EPSCM*CMSENE/2
      epsp = SQRT( Emin**2/P1(4)/P2(4) )
      epsq = SQRT( Emin**2/Q1(4)/Q2(4) )
      fYFSui = EXP( -gamp*log(1/epsp)  +gamp/4 -alf1/2 )
      fYFSli = EXP( -gamq*log(1/epsq)  +gamq/4 -alf1/2 )
      fYFSi  = fYFSui*fYFSli
! O(alf0)
      wtset(80) = delZB  /fYFSi/crude
! O(alf1) - no vacuum pol.
      wtset(81) = delZ1a /fYFSi/crude
! O(alf1) - with vac. pol.
      wtset(82) = delZ1v /fYFSi/crude    ! <---
*-----------------------------------------------------
C... weights below are mainly for comparisons with BABAMC
! O(alf0)
      wtset(90) = selZ0  /fYFSi/crude
! O(alf1) - no vacuum pol.
      wtset(91) = selZ1a /fYFSi/crude
! O(alf1) - with vac. pol.
      wtset(92) = selZ1v /fYFSi/crude
! Vacuum polarization correction
      wtset(97) = sevac  /fYFSi/crude
! Z self energy correction
      wtset(98) = seZse  /fYFSi/crude
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      ENDIF
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      END


