      SUBROUTINE MODL2A(MODE)
!     ***********************
!======================================================================!
!                                                                      !
!    OLD type Matrix element type (A)                                  !
!    similar as in BHLUMI version. 2.01  (Identical in first order)    !
!    but now also in the second order     (January. 94).               !
!    Exponentiated and not exponentiated  (January. 95)                !
!                                                                      !
!======================================================================!
!     ************************
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
      DIMENSION ph1(4),ph2(4)
      DIMENSION sfu(100),sfmu(100),bu10(100),bu11(100),bu11g(100)
      DIMENSION sfl(100),sfml(100),bl10(100),bl11(100),bl11g(100)
! ------------------ Inline functions ------------------
! Elements of single bremss. distribution
      CHI(x)= (1+(1-x)**2)/2
      XDELS(a,b,delt) =
     $ delt*(a**2+b**2)*(a**2+b**2)/((1-a)**2+(1-b)**2)/(a*b)
!==================================================================!
      DATA icont /0/
      icont=icont+1

      s = cmsene**2
      zeta = abs(tran)/s
!*****************************************************
!  Crude MC distribution (S-tilde factors omitted)
      CRUDE  =  s**2/tran**2
      facpq  =  s**2/(tranp*tranq)
      delta  =  amel**2/tran
      deltp  =  amel**2/tranp
      deltq  =  amel**2/tranq
      gam    =  2*alf1 *( dlog(1/delta)-1 )
      gamp   =  2*alf1 *( dlog(1/deltp)-1 )
      gamq   =  2*alf1 *( dlog(1/deltq)-1 )
      dis0   =  facpq *dix1(p1,p2,q1,q2)
!==================================================================!
      bt10u =0
      bt10l =0
      bt11u =0
      bt11l =0
      bt20u =0
      bt20l =0
      bt2ul =0
! ------------------
! beta0
! ------------------
      bt00   = dis0
      bt01   = bt00*(1 + gam )
      bt02   = bt00*(1 + gam +0.5*gam**2)
! ------------------------------------
! Contributions from beta1 UPPER line
! ------------------------------------
      DO i=1,nphot1
! numerically safe variables
         a   = al1(i)
         b   = be1(i)
         v   = a+b-a*b
         y   = a +b*deltp
         z   = b +a*deltp
! soft factor
         sfc  = a*b/(y*z)**2
         sfu(i) = sfc
! one photon bremss. distribution
         wm = a*b/(y*z) +deltp*(y**2+z**2)/((1-y)**2+(1-z)**2)*(y/z+z/y)
         sfmu(i) = wm/(y*z)
         dis10   = sfmu(i)*dis0
         dis11g  = dis10 *(1+gam + gam/2*log((1-b)/(1-a)))
         dis11   = dis10 *(1+gam + gam/4*log((1-b)/(1-a)))
! beta1 O(alf1),O(alf2)
         bu10(i) =  dis10   -bt00*sfc
         bu11(i) =  dis11   -bt01*sfc
         bu11g(i)=  dis11g  -bt01*sfc ! with LL bug
         bt10u   =  bt10u   +bu10(i)/sfc
         bt11u   =  bt11u   +bu11(i)/sfc
      ENDDO
! -------------------------------------
!  Contributions from beta1 LOWER line
! -------------------------------------
      DO i=1,nphot2
         a   = al2(i)
         b   = be2(i)
         v   = a+b-a*b
         y   = a+b*deltq
         z   = b+a*deltq
! soft factor
         sfc  = a*b/(y*z)**2
         sfl(i) = sfc
! one photon bremss. distribution
         wm = a*b/(y*z) +deltq*(y**2+z**2)/((1-y)**2+(1-z)**2)*(y/z+z/y)
         sfml(i) = wm/(y*z)
         dis10   = sfml(i)*dis0
         dis11g  = dis10 *(1+gam + gam/2*log((1-b)/(1-a)))
         dis11   = dis10 *(1+gam + gam/4*log((1-b)/(1-a)))
! beta1 O(alf1),O(alf2)
         bl10(i) =  dis10   -bt00*sfc
         bl11(i) =  dis11   -bt01*sfc
         bl11g(i)=  dis11g  -bt01*sfc ! with LL bug
         bt10l   =  bt10l   +bl10(i)/sfc
         bt11l   =  bt11l   +bl11(i)/sfc
      ENDDO

!------------------
! beta2 upper line
!------------------
      DO i=1,nphot1
      DO j=i+1,nphot1
! Define photon four-vectors for further local use
      DO  k=1,4
         ph1(k) = phot1(i,k)
         ph2(k) = phot1(j,k)
      ENDDO
! Basic variables first photon
      a1  = al1(i)
      b1  = be1(i)
      y1  = a1 +b1*deltp
      z1  = b1 +a1*deltp
      v1  = a1+b1-a1*b1
! Basic variables second photon
      a2  = al1(j)
      b2  = be1(j)
      y2  = a2 +b2*deltp
      z2  = b2 +a2*deltp
      v2  = a2+b2-a2*b2
! Soft factors, starred variables
      sf1  =  sfu(i)
      sf2  =  sfu(j)
      a1st = a1/(1-v2)
      b1st = b1/(1-v2)
      v1st = v1/(1-v2)
      a2st = a2/(1-v1)
      b2st = b2/(1-v1)
      v2st = v2/(1-v1)
      a1st = min(a1st,1-zeta)
      b1st = min(b1st,1-zeta)
      v1st = min(v1st,1-zeta)
      a2st = min(a2st,1-zeta)
      b2st = min(b2st,1-zeta)
      v2st = min(v2st,1-zeta)
! H-functions exact and Leading Log
      gs1    = sf1*chi(v1)
      gs2    = sf2*chi(v2)
      gs1st  = sf1*chi(v1st)
      gs2st  = sf2*chi(v2st)
! The actual two photon distributions
      IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
        IF ((a1-b1).gt.0d0) THEN
!          FIRST photon in INITIAL state, second in the final
           hs1f   = DIX2(2,ph2,p1,p2,q1,q2) *sfmu(i)
           hs2i   = DIX2(1,ph1,p1,p2,q1,q2) *sfmu(j)
           IF(v1.gt.v2) THEN
               d2uu  = facpq  *hs1f  *gs2
           ELSE
               d2uu  = facpq  *gs1st *hs2i
           ENDIF
        ELSE
!          SECOND photon in INITIAL state, first in the final
           hs1i   = DIX2(1,ph2,p1,p2,q1,q2) *sfmu(i)
           hs2f   = DIX2(2,ph1,p1,p2,q1,q2) *sfmu(j)
           IF(v1.gt.v2) THEN
               d2uu  = facpq  *hs1i *gs2st
           ELSE
               d2uu  = facpq  *gs1  *hs2f
           ENDIF
        ENDIF
      ELSE
!       SAME directions two photons
!         Both in INITIAL  state or Both in FINAL  state
          IF(v1.gt.v2) THEN
            hs1i  = DIX2(1,ph2,p1,p2,q1,q2) *sfmu(i)
            hs1f  = DIX2(2,ph2,p1,p2,q1,q2) *sfmu(i)
            d2uu  = facpq/2 *( hs1f*gs2st + hs1i*gs2)
          ELSE
            hs2i  = DIX2(1,ph1,p1,p2,q1,q2) *sfmu(j)
            hs2f  = DIX2(2,ph1,p1,p2,q1,q2) *sfmu(j)
            d2uu  = facpq/2 *( hs2f*gs1st + hs2i*gs1)
          ENDIF
      ENDIF
      bu20  = d2uu-bu10(i)*sfu(j)-bu10(j)*sfu(i)-bt00*sfu(j)*sfu(i)
      bt20u = bt20u + bu20/sfu(i)/sfu(j)
      ENDDO
      ENDDO
! ------------------
! beta2 lower line
! ------------------
      DO i=1,nphot2
      DO j=i+1,nphot2
! Define photon four-vectors for further local use
      DO  k=1,4
         ph1(k) = phot2(i,k)
         ph2(k) = phot2(j,k)
      ENDDO
! Basic variables first photon
      a1  = al2(i)
      b1  = be2(i)
      y1  = a1 +b1*deltq
      z1  = b1 +a1*deltq
      v1  = a1+b1-a1*b1
! Basic variables second photon
      a2  = al2(j)
      b2  = be2(j)
      y2  = a2 +b2*deltq
      z2  = b2 +a2*deltq
      v2  = a2+b2-a2*b2
! soft factors, starred variables
      sf1  =  sfl(i)
      sf2  =  sfl(j)
      a1st = a1/(1-v2)
      b1st = b1/(1-v2)
      v1st = v1/(1-v2)
      a2st = a2/(1-v1)
      b2st = b2/(1-v1)
      v2st = v2/(1-v1)
      a1st = min(a1st,1-zeta)
      b1st = min(b1st,1-zeta)
      v1st = min(v1st,1-zeta)
      a2st = min(a2st,1-zeta)
      b2st = min(b2st,1-zeta)
      v2st = min(v2st,1-zeta)
! H-functions exact and Leading Log
      gs1    = sf1*chi(v1)
      gs2    = sf2*chi(v2)
      gs1st  = sf1*chi(v1st)
      gs2st  = sf2*chi(v2st)
! The actual two photon distributions
      IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
        IF ((a1-b1).gt.0d0) THEN
!          FIRST photon in INITIAL state, second in the final
           hs1f   = DIX2(4,ph2,p1,p2,q1,q2) *sfml(i)
           hs2i   = DIX2(3,ph1,p1,p2,q1,q2) *sfml(j)
           IF(v1.gt.v2) THEN
               d2ll  = facpq  *hs1f  *gs2
           ELSE
               d2ll  = facpq  *gs1st*hs2i
           ENDIF
        ELSE
!          SECOND photon in INITIAL state, first in the final
           hs1i   = DIX2(3,ph2,p1,p2,q1,q2) *sfml(i)
           hs2f   = DIX2(4,ph1,p1,p2,q1,q2) *sfml(j)
           IF(v1.gt.v2) THEN
               d2ll  = facpq  *hs1i *gs2st
           ELSE
               d2ll  = facpq  *gs1 *hs2f
           ENDIF
        ENDIF
      ELSE
!       SAME directions two photons
!         Both in INITIAL  state or Both in FINAL  state
          IF(v1.gt.v2) THEN
            hs1i  = DIX2(3,ph2,p1,p2,q1,q2) *sfml(i)
            hs1f  = DIX2(4,ph2,p1,p2,q1,q2) *sfml(i)
            d2ll  = facpq/2 *( hs1f*gs2st + hs1i*gs2)
          ELSE
            hs2i  = DIX2(3,ph1,p1,p2,q1,q2) *sfml(j)
            hs2f  = DIX2(4,ph1,p1,p2,q1,q2) *sfml(j)
            d2ll  = facpq/2 *( hs2f*gs1st + hs2i*gs1)
          ENDIF
      ENDIF
      bl20  = d2ll-bl10(i)*sfl(j)-bl10(j)*sfl(i)-bt00*sfl(j)*sfl(i)
      bt20l = bt20l + bl20/sfl(i)/sfl(j)
      ENDDO
      ENDDO
!=================================================
! [11] One upper and one lower line photons
! Note that bt2ul=0 in the LL.
      DO i=1,nphot1
      DO j=1,nphot2
! Basic variables first photon
        a1  = al1(i)
        b1  = be1(i)
        v1  = a1 +b1-a1*b1
! Basic variables second photon
        a2  = al2(j)
        b2  = be2(j)
        v2  = a2+b2-a2*b2
! soft factors
        sf1  =  sfu(i)
        sf2  =  sfl(j)
        IF(v1.gt.v2) THEN
           d2ul= dis0 *sfmu(i) *sfl(j)
         ELSE
           d2ul= dis0 *sfu(i)  *sfml(j)
        ENDIF
        b2ul  = d2ul-bu10(i)*sfl(j)-sfu(i)*bl10(j)-bt00*sfu(i)*sfl(j)
        bt2ul = bt2ul + b2ul/sfu(i)/sfl(j)
      ENDDO
      ENDDO
     

!     **************************************
!     **  MC weights, with exponentiation **
!     **************************************
!          ---------------------------------
!          !    UPPER line ONLY,    (A)    !
!          ---------------------------------
! Case of upper line only has separate entries in WTSET because it is 
! analysed in the separate run with separate plotting programs 
! with different entries in vvrho, differerent histogram names etc.
!
! All beta's:  TOTAL O(alf0),O(alf1),O(alf2)
      WTSET( 30) =    bt00                /crude
      WTSET( 31) =   (bt01+bt10u)         /crude
      WTSET( 32) =   (bt02+bt11u +bt20u)  /crude
! Special bet0-type distr. (without chi0) for tests only 
      WTSET( 34) =   ( s**2/(tranp*tranq) )/crude
! Individual beta's in various orders.
! O(alf1)
      WTSET( 35) =   bt01  /crude                  ! beta0
      WTSET( 36) =   bt10u /crude                  ! beta1
! O(alf2)
      WTSET( 37) =   bt02  /crude                  ! beta0
      WTSET( 38) =   bt11u /crude                  ! beta1
      WTSET( 39) =   bt20u /crude                  ! beta2
!          ---------------------------------
!          !  UPPER + LOWER line,    (A)   !
!          ---------------------------------
! All beta's:  O(alf0),O(alf1),O(alf2)
      WTSET( 40) =    bt00  /crude
      WTSET( 41) =   (bt01+bt10u+bt10l) /crude
      WTSET( 42) =   (bt02+bt11u+bt11l+bt20u+bt20l+bt2ul) /crude
! Individual beta's in various orders.
! O(alf1)
      WTSET( 43) =    bt01   /crude                ! beta0
      WTSET( 44) =   (bt10u+bt10l) /crude          ! beta1
      WTSET( 45) =    bt10u  /crude                ! beta1 components
      WTSET( 46) =    bt10l  /crude
! O(alf2)
      WTSET( 47) =    bt02   /crude                ! beta0
      WTSET( 48) =   (bt11u+bt11l)/crude           ! beta1
      WTSET( 49) =   (bt20u+bt20l+bt2ul)/crude     ! beta2
      WTSET( 50) =    bt11u  /crude                ! beta1 components
      WTSET( 51) =    bt11l  /crude
      WTSET( 52) =    bt2ul  /crude
      WTSET( 53) =    bt20u  /crude                ! beta2 components
      WTSET( 54) =    bt20l  /crude

!==================================================================
!==================================================================
!                  Non-exponentiated version                      !
!==================================================================
!==================================================================
      pdel  = del*bcud(p1,p2,phsu1)
      qdel  = del*bcud(q1,q2,phsu2)
      fyfsu = exp(-gamp*dlog(1/pdel)  +gamp/4 -alf1/2)
      fyfsl = exp(-gamq*dlog(1/qdel)  +gamq/4 -alf1/2)
!
      bb0u    =0
      bb1u    =0
      bb2u    =0
      dd0u    =0
      dd1u    =0
      dd2u    =0
      IF(nphot1.eq.0) THEN
        bb0u   = dis0
        bb1u   = bb0u*(1 -gam*dlog(1/pdel))
        bb2u   = bb0u*(1 -gam*dlog(1/pdel) +1./2*(gam*dlog(1/pdel))**2)
        virt1  = 3./4*gam -1./2*alf1
        bb0u   = dis0
        dd1u   = dis0 *(1 -gam*dlog(1/pdel) +virt1)
      ELSEIF(nphot1.eq.1) THEN
        a   = al1(1)
        b   = be1(1)
        y   = a +b*deltp
        z   = b +a*deltp
        sfc = a*b/(y*z)**2
        bb1u = dis0
        wm  = a*b/(y*z) +deltp*(y**2+z**2)/((1-y)**2+(1-z)**2)*(y/z+z/y)
        sfm = wm/(y*z)
        dist1  = dis0 *sfm
        dd1u    = dist1/sfc
      ENDIF
!
      bb1  =0
      bb2  =0
!        **************************************
!        ****   Definitions of MC weights  ****
!        **************************************
! Case of upper-line-only has separate entries in WTSET because it is 
! analysed in the separate run with separate plotting programs 
! with different entries in vvrho, differerent histogram names etc.
!
! YFS formfactor fyfsu has to be removed in the unexponentiated case!
! Baseline
      WTSET(56) =    bb1u               /fyfsu/crude ! 1-photon
      WTSET(57) =    bb2u               /fyfsu/crude ! 2-photons
      WTSET(58) =    bb1                /fyfsu/crude ! 1-photon
      WTSET(59) =    bb2                /fyfsu/crude ! 2-photons
!          ---------------------------------
!          /////    UPPER line ONLY    /////
!          ---------------------------------
! TOTAL O(alf0),O(alf1),O(alf2)
      WTSET(60) =    dd0u                /fyfsu/crude
      WTSET(61) =    dd1u                /fyfsu/crude
      WTSET(62) =    dd2u                /fyfsu/crude
! Baseline
      WTSET(63) =    0
      WTSET(64) =    0
!************************************************************************

!    ////////////////////////////////////////////
!    /////  UNEXP   UPPER + LOWER line      /////
!    ////////////////////////////////////////////

! Entire 0,1,2-photon distributions
      dis0   =0
      dis1   =0
      dis2   =0
!
      virt1y = 6./4*gam -alf1
      virt2y = 0.5*(6./4*gam)**2 -3/2d0*alf1*gam
      h0     =  dix1(p1,p2,q1,q2)
!=================================================
      IF(nphot1.eq.0.and.nphot2.eq.0) THEN
! [00] No photons
!------
! O(alf0,1,2) Entire distribution 
        dis0= facpq*h0
        dis1= facpq*h0*(1+2*gam*dlog(del))
     $       +facpq*h0*virt1y
        dis2= facpq*h0*(1+2*gam*dlog(del)+0.5*(2*gam*dlog(del))**2)
     $       +facpq*h0*(1+2*gam*dlog(del))*virt1y
     $       +facpq*h0* virt2y
!=================================================
      ELSEIF(nphot1.eq.1.and.nphot2.eq.0) THEN
! [10] One upper line photon 
!  O(alf1,2) Entire distribution
!  Note that 1/2 in the log(1-v) term = {1/2 from two tree diagrams} 
        dis1   = facpq*h0*(1+xdels(al1(1),be1(1),deltp))
        dis2g  = dis1*(1 +(gamp+gamq)*dlog(del) +virt1y
     $        -1/2d0*gam*log(1-al1(1)) -1/2d0*gam*log(1-be1(1)) )
        dis2   = dis1*(1 +(gamp+gamq)*dlog(del) +virt1y
     $        -1/4d0*gam*log(1-al1(1)) -3/4d0*gam*log(1-be1(1))
     $       +1/4d0*(gamp-gam) -(gamp-gam)*log(1-be1(1))  ! emulate 4.02a.exp
     $    )
!=================================================
      ELSEIF(nphot1.eq.0.and.nphot2.eq.1) THEN
! [01] One lower line photon 
!------
! O(alf1,2) Entire distribution
! Note that 1/2 in the log(1-v) term = {1/2 from two tree diagrams} 
        dis1= facpq*h0*(1+xdels(al2(1),be2(1),deltq))
        dis2g= dis1*(1 +(gamp+gamq)*dlog(del) +virt1y
     $        -1/2d0*gam*log(1-al2(1)) -1/2d0*gam*log(1-be2(1)) )
        dis2 = dis1*(1 +(gamp+gamq)*dlog(del) +virt1y
     $        -1/4d0*gam*log(1-al2(1)) -3/4d0*gam*log(1-be2(1))
     $       +1/4d0*(gamq-gam) -(gamq-gam)*log(1-be2(1))   ! emulate 4.02a.exp
     $       )
!=================================================
      ELSEIF(nphot1.eq.2.and.nphot2.eq.0) THEN
! [20] Two upper line photons 
! Basic distribution
        a1  = al1(1)
        b1  = be1(1)
        a2  = al1(2)
        b2  = be1(2)
        v1  = a1+b1-a1*b1
        v2  = a2+b2-a2*b2
! Define photon four-vectors for further local use
        DO k=1,4
          ph1(k) = phot1(1,k)
          ph2(k) = phot1(2,k)
        ENDDO
        a1st = min(a1/(1-v2),1-zeta)
        b1st = min(b1/(1-v2),1-zeta)
        a2st = min(a2/(1-v1),1-zeta)
        b2st = min(b2/(1-v1),1-zeta)
        v1st = min(v1/(1-v2),1-zeta)
        v2st = min(v2/(1-v1),1-zeta)
!       H-functions exact
        hs1i = DIX2(1,ph2,p1,p2,q1,q2)*(1+xdels(a1,b1,deltp))
        hs1f = DIX2(2,ph2,p1,p2,q1,q2)*(1+xdels(a1,b1,deltp))
        hs2i = DIX2(1,ph1,p1,p2,q1,q2)*(1+xdels(a2,b2,deltp))
        hs2f = DIX2(2,ph1,p1,p2,q1,q2)*(1+xdels(a2,b2,deltp))
!       The actual two photon distributions
        IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
          IF (a1.gt.b1) THEN
!          FIRST photon in INITIAL state, second in the final
             IF(v1.gt.v2) dis2  = facpq* hs1f      *chi(v2)
             IF(v2.ge.v1) dis2  = facpq* chi(v1st) *hs2i
          ELSE
!          SECOND photon in INITIAL state, first in the final
             IF(v1.gt.v2) dis2  = facpq* hs1i    *chi(v2st)
             IF(v2.ge.v1) dis2  = facpq* chi(v1) *hs2f
          ENDIF
        ELSE
!       SAME directions two photons
          IF(v1.gt.v2) dis2 = facpq/2 *(hs1f*chi(v2st) + hs1i*chi(v2))
          IF(v2.ge.v1) dis2 = facpq/2 *(hs2f*chi(v1st) + hs2i*chi(v1))
        ENDIF
!=================================================
      ELSEIF(nphot1.eq.0.and.nphot2.eq.2) THEN
! [02] Two lower line photons 
! Basic distribution
        a1  = al2(1)
        b1  = be2(1)
        a2  = al2(2)
        b2  = be2(2)
        v1  = a1+b1-a1*b1
        v2  = a2+b2-a2*b2
        a1st = min(a1/(1-v2),1-zeta)
        b1st = min(b1/(1-v2),1-zeta)
        a2st = min(a2/(1-v1),1-zeta)
        b2st = min(b2/(1-v1),1-zeta)
        v1st = min(v1/(1-v2),1-zeta)
        v2st = min(v2/(1-v1),1-zeta)
! Define photon four-vectors for further local use
        DO k=1,4
          ph1(k) = phot2(1,k)
          ph2(k) = phot2(2,k)
        ENDDO
!       H-functions exact
        hs1i = DIX2(3,ph2,p1,p2,q1,q2)*(1+xdels(a1,b1,deltq))
        hs1f = DIX2(4,ph2,p1,p2,q1,q2)*(1+xdels(a1,b1,deltq))
        hs2i = DIX2(3,ph1,p1,p2,q1,q2)*(1+xdels(a2,b2,deltq))
        hs2f = DIX2(4,ph1,p1,p2,q1,q2)*(1+xdels(a2,b2,deltq))
!       The actual two photon distributions
        IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
          IF (a1.gt.b1) THEN
!          FIRST photon in INITIAL state, second in the final
             IF(v1.gt.v2) dis2  = facpq* hs1f      *chi(v2)
             IF(v2.ge.v1) dis2  = facpq* chi(v1st) *hs2i
          ELSE
!          SECOND photon in INITIAL state, first in the final
             IF(v1.gt.v2) dis2  = facpq* hs1i    *chi(v2st)
             IF(v2.ge.v1) dis2  = facpq* chi(v1) *hs2f
          ENDIF
        ELSE
!       SAME directions two photons
          IF(v1.gt.v2) dis2 = facpq/2 *(hs1f*chi(v2st) + hs1i*chi(v2))
          IF(v2.ge.v1) dis2 = facpq/2 *(hs2f*chi(v1st) + hs2i*chi(v1))
        ENDIF
!=================================================
      ELSEIF(nphot1.eq.1.and.nphot2.eq.1) THEN
! [11] One upper and lower line photons 
        a1  = al1(1)
        b1  = be1(1)
        a2  = al2(1)
        b2  = be2(1)
        v1  = a1+b1-a1*b1
        v2  = a2+b2-a2*b2
        IF(v1.gt.v2) THEN
          dis2 = facpq *dix1(p1,p2,q1,q2)*(1+xdels(a1,b1,deltp))
        ELSE
          dis2 = facpq *dix1(p1,p2,q1,q2)*(1+xdels(a2,b2,deltq))
        ENDIF
      ENDIF
!     =====

!        **************************************
!        ****   Definitions of MC weights  ****
!        ****    Non-exponentiated case    ****
!        **************************************
!          ---------------------------------
!          /////  UPPER + LOWER line   /////
!          ---------------------------------
 
      fyfs = fyfsu*fyfsl
! TOTAL O(alf0),O(alf1),O(alf2)
      WTSET(70) =    dis0                /fyfs/crude
      WTSET(71) =    dis1                /fyfs/crude
      WTSET(72) =    dis2                /fyfs/crude
!
      END

      FUNCTION DIX1(p1,p2,q1,q2)
!     ***************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION p1(4),p2(4),q1(4),q2(4)
! ------------------ Inline functions ------------------
      pdot(x1,y1,z1,e1,x2,y2,z2,e2) = 2*(e1*e2-x1*x2-y1*y2-z1*z2)
!
      s  = pdot(p1(1),p1(2),p1(3),p1(4) ,q1(1),q1(2),q1(3),q1(4))
      s1 = pdot(p2(1),p2(2),p2(3),p2(4) ,q2(1),q2(2),q2(3),q2(4))
      u  = pdot(p1(1),p1(2),p1(3),p1(4) ,q2(1),q2(2),q2(3),q2(4))
      u1 = pdot(p2(1),p2(2),p2(3),p2(4) ,q1(1),q1(2),q1(3),q1(4))
      dix1   = (s**2+u**2+s1**2+u1**2)/(4*s**2)
      END

      FUNCTION DIX2(noga,ph,p1,p2,q1,q2)
!     **********************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ph(4),p1(4),p2(4),q1(4),q2(4)
      DIMENSION w1(4),w2(4),v1(4),v2(4)
! ------------------ Inline functions ------------------
      pdot(x1,y1,z1,e1,x2,y2,z2,e2) = 2*(e1*e2-x1*x2-y1*y2-z1*z2)
!
      IF(noga.eq.1) THEN
        DO 1 k=1,4
 1      w1(k) = p1(k) - ph(k)
        s  = pdot(w1(1),w1(2),w1(3),w1(4) ,q1(1),q1(2),q1(3),q1(4))
        s1 = pdot(p2(1),p2(2),p2(3),p2(4) ,q2(1),q2(2),q2(3),q2(4))
        u  = pdot(w1(1),w1(2),w1(3),w1(4) ,q2(1),q2(2),q2(3),q2(4))
        u1 = pdot(p2(1),p2(2),p2(3),p2(4) ,q1(1),q1(2),q1(3),q1(4))
      ELSEIF(noga.eq.2) THEN
        DO 2 k=1,4
 2      w2(k) = p2(k) + ph(k)
        s  = pdot(p1(1),p1(2),p1(3),p1(4) ,q1(1),q1(2),q1(3),q1(4))
        s1 = pdot(w2(1),w2(2),w2(3),w2(4) ,q2(1),q2(2),q2(3),q2(4))
        u  = pdot(p1(1),p1(2),p1(3),p1(4) ,q2(1),q2(2),q2(3),q2(4))
        u1 = pdot(w2(1),w2(2),w2(3),w2(4) ,q1(1),q1(2),q1(3),q1(4))
      ELSEIF(noga.eq.3) THEN
        DO 3 k=1,4
 3      v1(k) = q1(k) - ph(k)
        s  = pdot(p1(1),p1(2),p1(3),p1(4) ,v1(1),v1(2),v1(3),v1(4))
        s1 = pdot(p2(1),p2(2),p2(3),p2(4) ,q2(1),q2(2),q2(3),q2(4))
        u  = pdot(p1(1),p1(2),p1(3),p1(4) ,q2(1),q2(2),q2(3),q2(4))
        u1 = pdot(p2(1),p2(2),p2(3),p2(4) ,v1(1),v1(2),v1(3),v1(4))
      ELSEIF(noga.eq.4) THEN
        DO 4 k=1,4
 4      v2(k) = q2(k) + ph(k)
        s  = pdot(p1(1),p1(2),p1(3),p1(4) ,q1(1),q1(2),q1(3),q1(4))
        s1 = pdot(p2(1),p2(2),p2(3),p2(4) ,v2(1),v2(2),v2(3),v2(4))
        u  = pdot(p1(1),p1(2),p1(3),p1(4) ,v2(1),v2(2),v2(3),v2(4))
        u1 = pdot(p2(1),p2(2),p2(3),p2(4) ,q1(1),q1(2),q1(3),q1(4))
      ELSE
        write(6,*) 'DIX2: Wrong noga = ',noga
      ENDIF
      dix2   = (s**2+u**2+s1**2+u1**2)/(4*s**2)
      END

