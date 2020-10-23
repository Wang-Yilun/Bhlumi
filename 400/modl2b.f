      SUBROUTINE MODL2B(MODE)
!     ***********************
!======================================================================!
!                                                                      !
!    NEW version  of the EXPONENTIATED matrix element, type (B)        !
!    Started May 92.                                                   !
!    O(alf2) part of the matrix element implemented Oct. 93            !
!    and fully debugged July 1995                                      !
!======================================================================!
!     ************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0, ALFINV=137.03604D0)
      PARAMETER( ALFA= 1D0/ALFINV, ALF1 = 1D0/PI/ALFINV)
      SAVE
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
      DIMENSION sfu(100),bu10(100),bu11(100),hsu(100)
      DIMENSION sfl(100),bl10(100),bl11(100),hsl(100)
      DIMENSION bu11g(100)
      DIMENSION bl11g(100)
      DIMENSION ph1(4)
! ------------------ Inline functions ------------------
! Dot products
! Elements of single bremss. distribution
      CHI(x)= (1+(1-x)**2)/2
      XDELS(a,b,delt) =
     $ delt*(a**2+b**2)*(a**2+b**2)/((1-a)**2+(1-b)**2)/(a*b)
      XRHOR(a,b,psi) =(1-a)*(1-b)+a*b
     $  +2*sqrt(a*b*max(1-a,0d0)*max(1-b,0d0))*cos(psi)
      XCHIR(a,b,rho,zeta) = ( 1      +(1-zeta/max(1-a,zeta))**2
     $                       +rho**2 +(rho -(1-b)*zeta)**2 )/4
!********** alternative ***********
!*     XSFCM(a,b,aa,bb,y,z,delt) = a*b/(y*z)**2
!*    $ +delt/(y*z)**2 *(y**2+z**2)*(a**2+b**2)/((1-aa)**2+(1-bb)**2)
!*     dis1   = facpq* xchir(a,b,xrhor(a,b,psi1),zeta)
!*    $              * xsfcm(a,b,a,b,y,z,deltp)
!**********************************
      DATA icont /0/

!==================================================================!
      s = cmsene**2
      zeta = abs(tran)/s
!*****************************************************
!  Crude MC distribution (S-tilde factors omitted)
      CRUDE  =  s**2/tran**2
      facpq  =  s**2/(tranp*tranq)
      zeta   =  abs(tran)/s
      chi0   =  ( 1+(1-zeta)**2 )/2
      delta  =  amel**2/tran
      deltp  =  amel**2/tranp
      deltq  =  amel**2/tranq
      gam    =  2*alf1 *( dlog(1/delta)-1 )
      gamp   =  2*alf1 *( dlog(1/deltp)-1 )
      gamq   =  2*alf1 *( dlog(1/deltq)-1 )
      psip   = psi1
      psiq   = psi2
!==================================================================!
!    ///////////////////////////////////////
!    /////     UPPER + LOWER line      /////
!    ///////////////////////////////////////
! ------------------
      bt10u =0
      bt10l =0
      bt11u =0
      bt11l =0
      bt20u =0
      bt20l =0
      bt2ul =0
! DEBUG
      bt11ug =0
      bt11lg =0
! ------------------
! beta0
! ------------------
      bt00   =  chi0* s**2 /(tranp*tranq)
      bt01   =  bt00*(1+gam)
      bt02   =  bt00*(1+gam+0.5*gam**2)
! ------------------
! beta1 upper line
! ------------------
      DO i=1,nphot1
         a   = al1(i)
         b   = be1(i)
         y   = a +b*deltp
         z   = b +a*deltp
         v   = a +b -a*b
! soft factor precisely as in mass weight from the generator
         sfc  = a*b/(y*z)**2
         sfu(i) = sfc
! one photon bremss. distribution
         chir    = xchir(a,b,xrhor(a,b,psip),zeta)
         dis10   = facpq * chir *(1+xdels(a,b,deltp)) * sfc
         dis11g  = dis10 *(1+gam + gam/2*log((1-b)/(1-a))) ! with LL bug
         dis11   = dis10 *(1+gam + gam/4*log((1-b)/(1-a)))
         hsu(i) = dis10 / facpq
! beta1 O(alf1),O(alf2)
         bu10(i) =  dis10   -bt00*sfc
         bu11(i) =  dis11   -bt01*sfc
         bt10u   =  bt10u   +bu10(i)/sfc
         bt11u   =  bt11u   +bu11(i)/sfc
! DEBUG
         bu11g(i) =  dis11g   -bt01*sfc ! with LL bug
         bt11ug   =  bt11ug   +bu11g(i)/sfc
      ENDDO
! ------------------
! beta1 lower line
! ------------------
      DO i=1,nphot2
         a   = al2(i)
         b   = be2(i)
         y   = a +b*deltq
         z   = b +a*deltq
         v   = a +b -a*b
! soft factor precisely as in mass weight from the generator
         sfc  = a*b/(y*z)**2
         sfl(i) = sfc
! one photon bremss. distribution
         chir    = xchir(a,b,xrhor(a,b,psiq),zeta)
         dis10   = facpq * chir *(1+xdels(a,b,deltq)) * sfc
         dis11g  = dis10 *(1+gam + gam/2*log((1-b)/(1-a)))
         dis11   = dis10 *(1+gam + gam/4*log((1-b)/(1-a)))
         hsl(i) = dis10 / facpq
! beta1 O(alf1),O(alf2)
         bl10(i) =  dis10   -bt00*sfc
         bl11(i) =  dis11   -bt01*sfc
         bt10l   =  bt10l   +bl10(i)/sfc
         bt11l   =  bt11l   +bl11(i)/sfc
! DEBUG
         bl11g(i) =  dis11g   -bt01*sfc ! with LL bug
         bt11lg   =  bt11lg   +bl11g(i)/sfc
      ENDDO
!------------------
! beta2 upper line
!------------------
      DO i=1,nphot1
      DO j=i+1,nphot1
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
      hs1    = hsu(i)
      hs2    = hsu(j)
      gs1st  = sf1*chi(v1st)
      gs2st  = sf2*chi(v2st)
      hs1st  = xchir(a1st,b1st,xrhor(a1st,b1st,psip),zeta)
     $            *(1+xdels(a1st,b1st,deltp))* sf1
      hs2st  = xchir(a2st,b2st,xrhor(a2st,b2st,psip),zeta)
     $            *(1+xdels(a2st,b2st,deltp))* sf2
! The actual two photon distributions
      IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
        IF (a1.gt.b1) THEN
!          FIRST photon in INITIAL state, second in the final
           IF(v1.gt.v2) THEN
               d2uu  = facpq  *hs1st*gs2
           ELSE
               d2uu  = facpq  *gs1st*hs2
           ENDIF
        ELSE
!          SECOND photon in INITIAL state, first in the final
           IF(v1.gt.v2) THEN
               d2uu  = facpq  *hs1*gs2st
           ELSE
               d2uu  = facpq  *gs1*hs2st
           ENDIF
        ENDIF
      ELSE
!       SAME directions two photons
        IF(v1.gt.v2) THEN
            d2uu  = facpq/2 *( hs1*gs2st + hs1st*gs2)
        ELSE
            d2uu  = facpq/2 *( hs2*gs1st + hs2st*gs1)
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
      hs1    = hsl(i)
      hs2    = hsl(j)
      gs1st  = sf1*chi(v1st)
      gs2st  = sf2*chi(v2st)
      hs1st  = xchir(a1st,b1st,xrhor(a1st,b1st,psiq),zeta)
     $            *(1+xdels(a1st,b1st,deltq))* sf1
      hs2st  = xchir(a2st,b2st,xrhor(a2st,b2st,psiq),zeta)
     $            *(1+xdels(a2st,b2st,deltq))* sf2
! The actual two photon distributions
      IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
        IF (a1.gt.b1) THEN
!          FIRST photon in INITIAL state, second in the final
           IF(v1.gt.v2) THEN
               d2ll  = facpq  *hs1st*gs2
           ELSE
               d2ll  = facpq  *gs1st*hs2
           ENDIF
        ELSE
!          SECOND photon in INITIAL state, first in the final
           IF(v1.gt.v2) THEN
               d2ll  = facpq  *hs1*gs2st
           ELSE
               d2ll  = facpq  *gs1*hs2st
           ENDIF
        ENDIF
      ELSE
!       SAME directions two photons
        IF(v1.gt.v2) THEN
           d2ll= facpq/2 *( hs1*gs2st  + hs1st*gs2)
       ELSE
           d2ll= facpq/2 *( hs2*gs1st  + hs2st*gs1)
        ENDIF
      ENDIF
      bl20  = d2ll-bl10(i)*sfl(j) -bl10(j)*sfl(i) -bt00*sfl(j)*sfl(i)
      bt20l = bt20l + bl20/sfl(i)/sfl(j)
      ENDDO
      ENDDO
!----------------------------
! beta2 upper and lower line
!----------------------------
      DO i=1,nphot1
      DO j=1,nphot2
! Basic variables first photon
      a1  = al1(i)
      b1  = be1(i)
      y1  = a1 +b1*deltp
      z1  = b1 +a1*deltp
      v1  = a1 +b1-a1*b1
! Basic variables second photon
      a2  = al2(j)
      b2  = be2(j)
      y2  = a2 +b2*deltq
      z2  = b2 +a2*deltq
      v2  = a2+b2-a2*b2
! soft factors
      sf1  =  sfu(i)
      sf2  =  sfl(j)
      IF(v1.gt.v2) THEN
         d2ul= facpq  *hsu(i)       *sf2*chi(v2)
       ELSE
         d2ul= facpq  *sf1*chi(v1)  *hsl(j)
      ENDIF
      b2ul  = d2ul-bu10(i)*sfl(j)-sfu(i)*bl10(j)-bt00*sfu(i)*sfl(j)
      bt2ul = bt2ul + b2ul/sfu(i)/sfl(j)
      ENDDO
      ENDDO

!        **************************************
!        ****   Definitions of MC weights  ****
!        **************************************
!          ---------------------------------
!          !    UPPER line ONLY,  New      !
!          ---------------------------------
! Case of upper line only has separate entries in WTSET because it is 
! analysed in the separate run with separate plotting programs 
! with different entries in vvrho, differerent histogram names etc.
!
! All beta's:  TOTAL O(alf0),O(alf1),O(alf2)
      WTSET(130) =    bt00                /crude
      WTSET(131) =   (bt01+bt10u)         /crude
      WTSET(132) =   (bt02+bt11u +bt20u)  /crude
! Individual beta's in various orders.
! Special O(alf0) bet0 for tests only (without chi0) called WTZERO
      WTSET(134) =   ( s**2/(tranp*tranq) )/crude
! O(alf1)
      WTSET(135) =   bt01  /crude                  ! beta0
      WTSET(136) =   bt10u /crude                  ! beta1
! O(alf2)
      WTSET(137) =   bt02  /crude                  ! beta0
      WTSET(138) =   bt11u /crude                  ! beta1
      WTSET(139) =   bt20u /crude                  ! beta2
! DEBUG
      WTSET(232) =   (bt02+bt11ug +bt20u)  /crude  ! total (LL bug)
      WTSET(238) =   bt11u /crude                  ! beta1 (LL bug)
!          ---------------------------------
!          !  UPPER + LOWER line, New      !
!          ---------------------------------
! All beta's:  O(alf0),O(alf1),O(alf2)
      WTSET(140) =    bt00  /crude
      WTSET(141) =   (bt01+bt10u+bt10l) /crude
      WTSET(142) =   (bt02+bt11u+bt11l+bt20u+bt20l+bt2ul) /crude
! Individual beta's in various orders.
! O(alf1)
      WTSET(143) =    bt01   /crude                ! beta0
      WTSET(144) =   (bt10u+bt10l) /crude          ! beta1
      WTSET(145) =    bt10u  /crude                ! beta1 components
      WTSET(146) =    bt10l  /crude
! O(alf2)
      WTSET(147) =    bt02   /crude                ! beta0
      WTSET(148) =   (bt11u+bt11l)/crude           ! beta1
      WTSET(149) =   (bt20u+bt20l+bt2ul)/crude     ! beta2
      WTSET(150) =    bt11u  /crude                ! beta1 components
      WTSET(151) =    bt11l  /crude
      WTSET(152) =    bt2ul  /crude
      WTSET(153) =    bt20u  /crude                ! beta2 components
      WTSET(154) =    bt20l  /crude
! DEBUG
      WTSET(242) = (bt02+bt11ug+bt11lg+bt20u+bt20l+bt2ul) /crude
      WTSET(250) =    bt11ug  /crude  ! (LL bug)
      WTSET(251) =    bt11l   /crude  ! (LL bug)
!==================================================================
!==================================================================
!                  Non-exponentiated version                      !
!==================================================================
!==================================================================
      pdel  = del*bcud(p1,p2,phsu1)
      qdel  = del*bcud(q1,q2,phsu2)
      fyfsu = exp(-gamp*dlog(1/pdel)  +gamp/4 -alf1/2)
      fyfsl = exp(-gamq*dlog(1/qdel)  +gamq/4 -alf1/2)
! temporary test CONV distr.
      bcon0u    =0
      bcon1u    =0
      bcon2u    =0
! entire distributions
      dis1uC   =0
      dis0u    =0
      dis1u    =0
      dis2u    =0
! beta-like components
      bt00u    =0
      bt01u    =0
      bt02u    =0
      bt10u    =0
      bt11u    =0
      bt20u    =0
! DEBUG part
      dis2ug   = 0
      bt11ug   = 0
!------------------------------------
! Normal one-line virtual correction
      virt1   =  3d0/4*gam -1d0/2*alf1 
! Tweo lines with one YFS formfactor subtraction
! that is: 5/4gam = 2(3/4)gam -1/4gam
      virt1x = 5d0/4*gam-1d0/2*alf1
      virt2x = 0.5d0*(5d0/4*gam)**2
! Note that in order to shorten our formulas we omitt S-factors,
! They always cancel with S-factors from crude distributions
!    ///////////////////////////////////////
!    /////  UNEXP  UPPER line only     /////
!    ///////////////////////////////////////
      IF(nphot1.eq.0) THEN
!     ====================
! [00] No photons, note here gamp=gam
!------
! Clasical O(alf0,1)
        dis0uC  = facpq*chi0 
        dis1uC  = facpq*chi0 *(1 +virt1 +gamp*dlog(del))
! NCONV  
        bcon0u= facpq*chi0 
        bcon1u= facpq*chi0*(1 +gamp*dlog(del))
        bcon2u= facpq*chi0*(1 +gamp*dlog(del) +1./2*(gamp*dlog(del))**2)
! O(alf0,1,2) Entire distribution 
        dis0u = facpq*chi0
        dis1u = facpq*chi0*(1 +gamp*dlog(del))
     $         +facpq*chi0*virt1x
        dis2u = facpq*chi0*(1+gamp*dlog(del) +1./2*(gamp*dlog(del))**2)
     $         +facpq*chi0*(1+gamp*dlog(del))*virt1x
     $         +facpq*chi0*virt2x
!------
! Beta's for controll
        bt00u= dis0u
        bt01u= dis1u
        bt02u= dis2u
! DEBUG part
        dis2ug = dis2u
      ELSEIF(nphot1.eq.1) THEN
!     ========================
! [10] One upper line photon, note here gamp.ne.gam !!!
!------
! Basic distribution
        a   = al1(1)
        b   = be1(1)
        v   = a+b-a*b
! O(alf1) Entire distributions 
        hs1  = xchir(a,b,xrhor(a,b,psip),zeta)*(1+xdels(a,b,deltp))
!------
! Clasical O(alf1)
        dis1uC = facpq *hs1
! NCONV 
        bcon1u   = facpq*chi0
        bcon2u   = facpq*chi0 *(1 +gamp*dlog(del))
!------
! O(alf1,2) Entire distribution
        dis1u  = facpq*hs1
        dis2ug = facpq*hs1*(1 +gamp*dlog(del) +virt1x 
     $       -1/2d0*gam*log(1-a) -1/2d0*gam*log(1-b) )
        dis2u  = facpq*hs1*(1 +gamp*dlog(del) +virt1x 
     $       -1/4d0*gam*log(1-a) -3/4d0*gam*log(1-b) )
! Beta's for control
        bt01u  = facpq*chi0
        bt02u  = facpq*chi0 *(1 +gamp*dlog(del) +virt1x)
        bt10u  = dis1u  -bt01u
        bt11u  = dis2u  -bt02u
        bt11ug = dis2ug -bt02u ! with LL bug

      ELSEIF(nphot1.eq.2) THEN
!     ========================
! [20] Two upper line photon 
!------
! NCONV 
        bcon2u   = facpq*chi0
!------
! Basic distribution
        a1  = al1(1)
        b1  = be1(1)
        a2  = al1(2)
        b2  = be1(2)
        v1  = a1+b1-a1*b1
        v2  = a2+b2-a2*b2
        a1st = min(a1/(1-v2),1-zeta)
        b1st = min(b1/(1-v2),1-zeta)
        a2st = min(a2/(1-v1),1-zeta)
        b2st = min(b2/(1-v1),1-zeta)
        v1st = min(v1/(1-v2),1-zeta)
        v2st = min(v2/(1-v1),1-zeta)
!       H-functions exact and Leading Log
        hs1 =xchir(a1,b1,xrhor(a1,b1,psip),zeta)*(1+xdels(a1,b1,deltp))
        hs2 =xchir(a2,b2,xrhor(a2,b2,psip),zeta)*(1+xdels(a2,b2,deltp))
        hs1st  = xchir(a1st,b1st,xrhor(a1st,b1st,psip),zeta)
     $            *(1+xdels(a1st,b1st,deltp))
        hs2st  = xchir(a2st,b2st,xrhor(a2st,b2st,psip),zeta)
     $            *(1+xdels(a2st,b2st,deltp))
!       The actual two photon distributions
        IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
          IF (a1.gt.b1) THEN
!          FIRST photon in INITIAL state, second in the final
             IF(v1.gt.v2) d2u  = hs1st*chi(v2)
             IF(v2.ge.v1) d2u  = chi(v1st)*hs2
          ELSE
!          SECOND photon in INITIAL state, first in the final
             IF(v1.gt.v2) d2u  = hs1*chi(v2st)
             IF(v2.ge.v1) d2u  = chi(v1)*hs2st
          ENDIF
        ELSE
!       SAME directions two photons
          IF(v1.gt.v2) d2u  = 1./2 *( hs1*chi(v2st) + hs1st*chi(v2))
          IF(v2.ge.v1) d2u  = 1./2 *( hs2*chi(v1st) + hs2st*chi(v1))
        ENDIF
!------
! O(alf2) Entire distribution dis2u 
        dis2u = facpq*d2u
!------
! Beta's for control
        bt02u = facpq*chi0
        bt11u = facpq*( (hs1-chi0) +(hs2-chi0) )
        bt20u = dis2u -bt11u  -bt02u
! DEBUG part
        dis2ug = dis2u
        bt11ug = bt11u ! with LL bug
      ENDIF
!     =====
!        **************************************
!        ****   Definitions of MC weights  ****
!        ****    Non-exponentiated case    ****
!        **************************************
! Case of upper-line-only has separate entries in WTSET because it is 
! analysed in the separate run with separate plotting programs 
! with different entries in vvrho, differerent histogram names etc.
!
! YFS formfactor fyfsu has to be removed in the unexponentiated case!
!-----
! NCONV for relative comparisons
      WTSET(156) =   bcon1u      /fyfsu/crude ! NCONV O(alf1)
      WTSET(157) =   bcon2u      /fyfsu/crude ! NCONV O(alf2)
!          ---------------------------------
!          /////    UPPER line ONLY    /////
!          ---------------------------------
! TOTAL O(alf0),O(alf1),O(alf2)
      WTSET(160) =   bt00u       /fyfsu/crude
      WTSET(161) =   dis1u       /fyfsu/crude
      WTSET(162) =   dis2u       /fyfsu/crude
      WTSET(164) =   dis1uC      /fyfsu/crude
! Individual beta's in various orders.
! O(alf1)
      WTSET(165) =   bt01u  /fyfsu/crude            ! bt0
      WTSET(166) =   bt10u  /fyfsu/crude            ! bt1
! O(alf2)
      WTSET(167) =   bt02u  /fyfsu/crude            ! bt0
      WTSET(168) =   bt11u  /fyfsu/crude            ! bt1
      WTSET(169) =   bt20u  /fyfsu/crude            ! bt2
! Debug part
      WTSET(262) =   dis2ug   /fyfsu/crude    ! O(alf2) tot BUG
      WTSET(290) =   bt11ug   /fyfsu/crude    ! O(alf2) bt1 BUG

!    ////////////////////////////////////////////
!    /////  UNEXP   UPPER + LOWER line      /////
!    ////////////////////////////////////////////
! CONV for primary tests
      bcon0  =0
      bcon1  =0
      bcon2  =0
! Entire 0,1,2-photon distributions
      dis0   =0
      dis1   =0
      dis2   =0
! beta-like components
      bt00    =0
      bt01    =0
      bt02    =0
      bt10    =0
      bt11    =0
      bt20    =0
! ---------- salami on second order beta1
      bt11rg =0
      bt11uu =0
      bt11ll =0
      bt11ul =0
! ----------salami on second order beta2
      bt20u =0
      bt20l =0
      bt2ul =0
! DEBUG
      dis2g = 0
!
      virt1y = 3/2d0*gam -alf1
      virt2y = 1/2d0*(3/2d0*gam)**2 -3/2d0*alf1*gam
!=================================================
      IF(nphot1.eq.0.and.nphot2.eq.0) THEN
! [00] No photons
!------CONV 
        bcon0=facpq*chi0
        bcon1=bcon0*(1 -2*gam*dlog(1/del))
        bcon2=bcon0*(1 -2*gam*dlog(1/del)+0.5*(2*gam*dlog(1/del))**2 )
!------
! O(alf0,1,2) Entire distribution 
        dis0= facpq*chi0
        dis1= facpq*chi0*(1+2*gam*dlog(del))
     $       +facpq*chi0*virt1y
        dis2= facpq*chi0*(1+2*gam*dlog(del)+0.5*(2*gam*dlog(del))**2)
     $       +facpq*chi0*(1+2*gam*dlog(del))*virt1y
     $       +facpq*chi0* virt2y
!------
! Beta's for controll
        bt00= dis0
        bt01= dis1
        bt02= dis2
! DEBUG
        dis2g = dis2
!=================================================
      ELSEIF(nphot1.eq.1.and.nphot2.eq.0) THEN
! [10] One upper line photon 
!------CONV
        bcon1  = facpq*chi0
        bcon2  = bcon1  *(1 -(gamp+gamq)*dlog(1/del))
!------Basic distribution
        hs1 = (1+xdels(al1(1),be1(1),deltp))
     $        *xchir(al1(1),be1(1),xrhor(al1(1),be1(1),psip),zeta)
!------
! O(alf1,2) Entire distribution
! Note that 1/2 in the log(1-v) term = {1/2 from two tree diagrams} 
        dis1 = facpq*hs1
        dis2g= facpq*hs1*(1 +(gamp+gamq)*dlog(del) +virt1y
     $       -1/2d0*gam*log(1-al1(1)) -1/2d0*gam*log(1-be1(1)) )
        dis2 = facpq*hs1*(1 +(gamp+gamq)*dlog(del) +virt1y
     $       -1/4d0*gam*log(1-al1(1)) -3/4d0*gam*log(1-be1(1))
     $       +1/4d0*(gamp-gam) -(gamp-gam)*log(1-be1(1))  ! emulate 4.02a.exp
     $       )
! Beta's for control
        bt01 = facpq*chi0
        bt02 = facpq*chi0 *(1 +(gamp+gamq)*dlog(del) +virt1y)
        bt10 = dis1  -bt01
        bt11 = dis2  -bt02
! DEBUG part
        bt11g= dis2g -bt02 ! with LL bug
!/////////  SALAMI on bet11 ///////// 
! bet11 rewritten as bt10 +bt11rg +bt11uu +bt11ll +bt11ul
!        bt11 = facpq*(
!     $  (hs1-chi0)                            ! first order
!     $ +(hs1-chi0)*virt1y -hs1*gam/2*log(1-v) ! regular part + bug corr.
!     $ +(hs1-chi0)*gamp*dlog(del)             ! re_uper-im_uper
!     $ +(hs1-chi0)*gamq*dlog(del)             ! re_uper-im_lower
!     $)
        bt11rg = facpq*((hs1-chi0)*virt1y 
     $  +hs1*( -1/4d0*gam*log(1-al1(1)) -3/4d0*gam*log(1-be1(1)) )
     $       +1/4d0*(gamp-gam) -(gamp-gam)*log(1-be1(1))  ! emulate 4.02a.exp
     $       )
        bt11uu = facpq*(hs1-chi0)*gamp*dlog(del)
        bt11ul = facpq*(hs1-chi0)*gamq*dlog(del)

!=================================================
      ELSEIF(nphot1.eq.0.and.nphot2.eq.1) THEN
! [01] One lower line photon 
!------CONV 
        bcon1  = facpq*chi0
        bcon2  = bcon1  *(1 -(gamp+gamq)*dlog(1/del))
!------Basic distribution
        hs1 = (1+xdels(al2(1),be2(1),deltq))
     $        *xchir(al2(1),be2(1),xrhor(al2(1),be2(1),psiq),zeta)
!------
! O(alf1,2) Entire distribution
! Note that 1/2 in the log(1-v) term = {1/2 from two tree diagrams} 
        dis1 = facpq*hs1
        dis2g= facpq*hs1*(1 +(gamp+gamq)*dlog(del) +virt1y
     $       -1/2d0*gam*log(1-al2(1)) -1/2d0*gam*log(1-be2(1)) )
        dis2 = facpq*hs1*(1 +(gamp+gamq)*dlog(del) +virt1y
     $       -1/4d0*gam*log(1-al2(1)) -3/4d0*gam*log(1-be2(1)) 
     $       +1/4d0*(gamq-gam) -(gamq-gam)*log(1-be2(1))   ! emulate 4.02a.exp
     $       )
! Beta's for control
        bt01 = facpq*chi0
        bt02 = facpq*chi0 *(1 +(gamp+gamq)*dlog(del) +virt1y)
        bt10 = dis1  -bt01
        bt11g= dis2g -bt02 ! with LL bug
        bt11 = dis2  -bt02
!/////////  SALAMI on bet11 ///////// 
        bt11rg = facpq*((hs1-chi0)*virt1y 
     $   +hs1*(-1/4d0*gam*log(1-al2(1)) -3/4d0*gam*log(1-be2(1)) )
     $       +1/4d0*(gamq-gam) -(gamq-gam)*log(1-be2(1))   ! emulate 4.02a.exp
     $       )

        bt11ll = facpq*(hs1-chi0)*gamq*dlog(del)
        bt11ul = facpq*(hs1-chi0)*gamp*dlog(del)
!=================================================
      ELSEIF(nphot1.eq.2.and.nphot2.eq.0) THEN
! [20] Two upper line photons 
!------CONV
        bcon2  = facpq*chi0
!------
! Basic distribution
        a1  = al1(1)
        b1  = be1(1)
        a2  = al1(2)
        b2  = be1(2)
        v1  = a1+b1-a1*b1
        v2  = a2+b2-a2*b2
        a1st = min(a1/(1-v2),1-zeta)
        b1st = min(b1/(1-v2),1-zeta)
        a2st = min(a2/(1-v1),1-zeta)
        b2st = min(b2/(1-v1),1-zeta)
        v1st = min(v1/(1-v2),1-zeta)
        v2st = min(v2/(1-v1),1-zeta)
!       H-functions exact and Leading Log
        hs1 =xchir(a1,b1,xrhor(a1,b1,psip),zeta)*(1+xdels(a1,b1,deltp))
        hs2 =xchir(a2,b2,xrhor(a2,b2,psip),zeta)*(1+xdels(a2,b2,deltp))
        hs1st  = xchir(a1st,b1st,xrhor(a1st,b1st,psip),zeta)
     $            *(1+xdels(a1st,b1st,deltp))
        hs2st  = xchir(a2st,b2st,xrhor(a2st,b2st,psip),zeta)
     $            *(1+xdels(a2st,b2st,deltp))
!       The actual two photon distributions
        IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
          IF (a1.gt.b1) THEN
!          FIRST photon in INITIAL state, second in the final
             IF(v1.gt.v2) d2u  = hs1st*chi(v2)
             IF(v2.ge.v1) d2u  = chi(v1st)*hs2
          ELSE
!          SECOND photon in INITIAL state, first in the final
             IF(v1.gt.v2) d2u  = hs1*chi(v2st)
             IF(v2.ge.v1) d2u  = chi(v1)*hs2st
          ENDIF
        ELSE
!       SAME directions two photons
          IF(v1.gt.v2) d2u  = 1./2 *( hs1*chi(v2st) + hs1st*chi(v2))
          IF(v2.ge.v1) d2u  = 1./2 *( hs2*chi(v1st) + hs2st*chi(v1))
        ENDIF
!------
! O(alf2) Entire distribution dis2u 
        dis2 = facpq*d2u
!------
! Beta's for control
        bt02 = facpq*chi0
        bt11 = facpq*( (hs1-chi0) +(hs2-chi0) )
        bt20 = dis2 -bt11  -bt02
        bt20u= bt20
!/////////  SALAMI on bet11 ///////// 
        bt11uu = bt11
! DEBUG
        dis2g = dis2
!=================================================
      ELSEIF(nphot1.eq.0.and.nphot2.eq.2) THEN
! [02] Two lower line photons 
!------CONV
        bcon2  = facpq*chi0
!------
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
!       H-functions exact and Leading Log
        hs1 =xchir(a1,b1,xrhor(a1,b1,psiq),zeta)*(1+xdels(a1,b1,deltq))
        hs2 =xchir(a2,b2,xrhor(a2,b2,psiq),zeta)*(1+xdels(a2,b2,deltq))
        hs1st  = xchir(a1st,b1st,xrhor(a1st,b1st,psiq),zeta)
     $            *(1+xdels(a1st,b1st,deltq))
        hs2st  = xchir(a2st,b2st,xrhor(a2st,b2st,psiq),zeta)
     $            *(1+xdels(a2st,b2st,deltq))
!       The actual two photon distributions
        IF ( (a1-b1)*(a2-b2).le.0d0 ) THEN
!       OPPOSITE directions two photons
          IF (a1.gt.b1) THEN
!          FIRST photon in INITIAL state, second in the final
             IF(v1.gt.v2) d2u  = hs1st*chi(v2)
             IF(v2.ge.v1) d2u  = chi(v1st)*hs2
          ELSE
!          SECOND photon in INITIAL state, first in the final
             IF(v1.gt.v2) d2u  = hs1*chi(v2st)
             IF(v2.ge.v1) d2u  = chi(v1)*hs2st
          ENDIF
        ELSE
!       SAME directions two photons
          IF(v1.gt.v2) d2u  = 1./2 *( hs1*chi(v2st) + hs1st*chi(v2))
          IF(v2.ge.v1) d2u  = 1./2 *( hs2*chi(v1st) + hs2st*chi(v1))
        ENDIF
!------
! O(alf2) Entire distribution dis2u 
        dis2 = facpq*d2u
!------
! Beta's for control
        bt02 = facpq*chi0
        bt11 = facpq*( (hs1-chi0) +(hs2-chi0) )
        bt20 = dis2 -bt11  -bt02
        bt20l= bt20
c/////////  SALAMI on bet11 ///////// 
        bt11ll = bt11
! DEBUG
        dis2g = dis2
!=================================================
      ELSEIF(nphot1.eq.1.and.nphot2.eq.1) THEN
! [11] One upper and one lower line photons 
!------CONV
        bcon2  = facpq*chi0
!------
! O(alf2) Entire distribution dis2u 
      a1  = al1(1)
      b1  = be1(1)
      a2  = al2(1)
      b2  = be2(1)
      v1  = a1+b1-a1*b1
      v2  = a2+b2-a2*b2
      hs1 = xchir(a1,b1,xrhor(a1,b1,psip),zeta)*(1+xdels(a1,b1,deltp))
      hs2 = xchir(a2,b2,xrhor(a2,b2,psiq),zeta)*(1+xdels(a2,b2,deltq))
      IF(v1.gt.v2) THEN
        d2ul=  hs1     *chi(v2)
      ELSE
        d2ul=  chi(v1) * hs2
      ENDIF
        dis2 = facpq*d2ul
!------
! Beta's for control
        bt02  = facpq*chi0
        bt11  = facpq*( (hs1-chi0) +(hs2-chi0) )
        bt20  = dis2 -bt11 -bt02
        bt2ul = bt20
!/////////  SALAMI on bet11 ///////// 
        bt11ul = bt11
! DEBUG
        dis2g = dis2
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
      WTSET(170) =    dis0                /fyfs/crude
      WTSET(171) =    dis1                /fyfs/crude
      WTSET(172) =    dis2                /fyfs/crude
!***************************************************************
! Individual beta's in various orders.
! O(alf1)
      WTSET(173) =    bt01         /fyfs/crude       ! beta0
      WTSET(174) =    bt10         /fyfs/crude       ! beta1
! O(alf2)
      WTSET(175) =    bt11uu       /fyfs/crude       ! beta11 up-up
      WTSET(176) =    bt11ll       /fyfs/crude       ! beta11 low-low
      WTSET(177) =    bt02         /fyfs/crude       ! beta0
      WTSET(178) =    bt11         /fyfs/crude       ! beta1
      WTSET(179) =    bt20         /fyfs/crude       ! beta2
      WTSET(180) =    bt11rg       /fyfs/crude       ! beta11 regular 
      WTSET(181) =    bt11ul       /fyfs/crude       ! beta11 up-low
      WTSET(182) =    bt2ul        /fyfs/crude
      WTSET(183) =    bt20u        /fyfs/crude       ! beta2 components
      WTSET(184) =    bt20l        /fyfs/crude
! CONV distribution for tests
      WTSET(185) =   bcon0         /fyfs/crude ! NCONV O(alf0)
      WTSET(186) =   bcon1         /fyfs/crude ! NCONV O(alf1)
      WTSET(187) =   bcon2         /fyfs/crude ! NCONV O(alf2)
! DEBUG part
      WTSET(272) =    dis2g        /fyfs/crude ! Total with LL bug
!***************************************************************
      END
