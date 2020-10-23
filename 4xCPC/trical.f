************************************************************************
!      Collection of trigger (selection) routines
!      used during 1995/96 workshop, see also silicon.f
************************************************************************
************************************************************************
! =====> trical1 is  BERE1 CALO1 and CALO2
! =====> trical2 is  idealized LCAL of TH.6118 (1991, PLB268)
************************************************************************
************************************************************************

      SUBROUTINE trical1
     $     (th1n,th2n,th1w,th2w,delcon,dlthe,dlphi,zcalA,zcalB,zcalC)
!     ******************************************************************
! ----------------------------------------------------------------------
! LCAL type trigger with theta-rings and phi-sectors
! Descendant of TRIGAS0, see zcalA,
! ====> zcalA is for BERE1 <==========
! ====> zcalB is for CALO1 <==========
! ====> zcalC is for CALO2 <==========
! CALO1: Associate photon energy with electrons wthin delcon cone
! CALO2: Associate photon energy with electrons 
!        wthin dlthe*dlphi plaquette
! ----------------------------------------------------------------------
C     ****************************************** 
      IMPLICIT REAL*8(A-H,O-Z) 
      SAVE
      PARAMETER( pi = 3.1415926535897932D0 )
      DIMENSION zcalA(*),zcalB(*),zcalC(*)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT 
      LOGICAL LFORW,LFORN,LBACW,LBACN
      DATA icont /0/

      IF(icont.eq.0) THEN
        icont=icont+1
        WRITE(6,*) '*******************************'
        WRITE(6,*) '***     trical1 Trigger     ***'
        WRITE(6,*) '*******************************'
        WRITE(6,'(5A12/4F12.9,I12)')
     $  ' th1w',' th2w',' th1n',' th2n','nphi',
     $    th1w ,  th2w ,  th1n ,  th2n , nphi
      ENDIF

! wide/narrow sectors forward
      thetp=angfi(p2(3),dsqrt(p2(1)**2+p2(2)**2))
      phip =angfi(p2(1),p2(2))
      lforw= thetp .GT. th1w .AND. thetp .LT. th2w
      lforn= thetp .GT. th1n .AND. thetp .LT. th2n
! wide/narrow sectors backward
      thetq=angfi(-q2(3),dsqrt(q2(1)**2+q2(2)**2))
      phiq =angfi(-q2(1),-q2(2))
      lbacw= thetq .GT. th1w .AND. thetq .LT. th2w
      lbacn= thetq .GT. th1n .AND. thetq .LT. th2n
!--------------------------------------------------------
! BARE1: Clasical choice as in trigas1 of TH.6118
      sv1 = (p2(4)+q2(4))**2 -(p2(3)+q2(3))**2
     $     -(p2(2)+q2(2))**2 -(p2(1)+q2(1))**2 
      sv  = (p1(4)+q1(4))**2 -(p1(3)+q1(3))**2
     $     -(p1(2)+q1(2))**2 -(p1(1)+q1(1))**2  
      z1 = sv1/sv   
      xwide= 0D0
      xnarr= 0D0
      xmix1= 0D0
      xmix2= 0D0
      IF(lforw .AND. lbacw) xwide= z1
      IF(lforn .AND. lbacn) xnarr= z1
      IF(lforw .AND. lbacn) xmix1= z1
      IF(lforn .AND. lbacw) xmix2= z1
!--------------------------------------------------------
! CALO1: Associate photon energy with electrons wthin delcon cone
      e1=p2(4)
      e2=q2(4)
      DO k=1,NPHOT
        cosph1= (p2(1)*phot(k,1)+p2(2)*phot(k,2)+p2(3)*phot(k,3))
     $         /phot(k,4)/p2(4)
        angph1=acos(min(cosph1,1d0))
        IF(angph1.lt.delcon) e1=e1+phot(k,4)
        cosph2= (q2(1)*phot(k,1)+q2(2)*phot(k,2)+q2(3)*phot(k,3))
     $         /phot(k,4)/q2(4)
        angph2=acos(min(cosph2,1d0))
        IF(angph2.lt.delcon) e2=e2+phot(k,4)
      ENDDO
      ebeam = p1(4)
      z2 = e1*e2/ebeam**2
      z2ww=0d0
      z2nn=0d0
      z2wn=0d0
      z2nw=0d0
      IF(lforw .AND. lbacw) z2ww= z2
      IF(lforn .AND. lbacn) z2nn= z2
      IF(lforw .AND. lbacn) z2wn= z2
      IF(lforn .AND. lbacw) z2nw= z2
!--------------------------------------------------------
! CALO2: Associate photon energy with electrons 
! wthin dlthe*dlphi semi-ring
      e1=p2(4)
      e2=q2(4)
      DO k=1,NPHOT
! forward
        thet1=angfi( phot(k,3),dsqrt(phot(k,1)**2+phot(k,2)**2))
        phi1 =angfi( phot(k,1),phot(k,2))
        difphi1 = min(ABS(phi1 -phip),2*pi-ABS(phi1 -phip))
        IF( abs(thet1-thetp) .LT. dlthe .AND.
     $               difphi1 .LT. dlphi)  e1=e1+phot(k,4)
! backward
        thet2=angfi(-phot(k,3),dsqrt(phot(k,1)**2+phot(k,2)**2))
        phi2 =angfi(-phot(k,1),-phot(k,2))
        difphi2 = min(ABS(phi2 -phiq),2*pi-ABS(phi2 -phiq))
        IF( abs(thet2-thetq) .LT. dlthe .AND.
     $               difphi2 .LT. dlphi)  e2=e2+phot(k,4)
      ENDDO
      ebeam = p1(4)
      z3 = e1*e2/ebeam**2
      z3ww=0d0
      z3nn=0d0
      z3wn=0d0
      z3nw=0d0
      IF(lforw .AND. lbacw) z3ww= z3
      IF(lforn .AND. lbacn) z3nn= z3
      IF(lforw .AND. lbacn) z3wn= z3
      IF(lforn .AND. lbacw) z3nw= z3
!--------------------------------------------------------
! BARE1: Clasical choice as in trigas1 of TH.6118
      zcalA(1) = xwide
      zcalA(2) = xnarr
      zcalA(3) = xmix1
      zcalA(4) = xmix2
! CALO1 choice with cones around electrons
      zcalB(1) = z2ww
      zcalB(2) = z2nn
      zcalB(3) = z2wn
      zcalB(4) = z2nw
! CALO2 choice with rectanguler semi-rings
      zcalC(1) = z3ww
      zcalC(2) = z3nn
      zcalC(3) = z3wn
      zcalC(4) = z3nw

      END

************************************************************************
************************************************************************
************************************************************************
************************************************************************
************************************************************************
************************************************************************
************************************************************************

      SUBROUTINE trical2
     $   (nphi,th1n,th2n,th1w,th2w,zcal,zcal1)
!     **********************************************************
! LCAL Idealized exper. CALORIMETRIC trigger on dressed final electrons.
! Electrons and photons not distinguished!
! This was used initially in paper TH.6118.
! Input:  th1n,th2n,th1w,th2w  theta limits of the Narrow/Wide cone
!         nphi                 number of phi-sectors, nphi=1 is possible
! Output: vcal(i),i=1,2,3,4:  
!                each entry parametrizes total enery 'cought' 
!                in the corresponding N-N, W-W, N-W, W-N angular range. 
! The energy cut X>X_min is to be imposed outside this program.
!     ****************************************** 
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      PARAMETER( pi = 3.1415926535897932D0 )
      PARAMETER( nmx1 = 20)
      DIMENSION zcal(*),zcal1(*)
      LOGICAL langw,langn,lphi
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      DIMENSION pc(100,4)
      DIMENSION pcw1(nmx1,4),pcn1(nmx1,4)
      DIMENSION pcw2(nmx1,4),pcn2(nmx1,4)
      DATA icont /0/

      IF(icont.eq.0) THEN
        icont=icont+1
        WRITE(6,*) '*******************************'
        WRITE(6,*) '***     trical2 Trigger     ***'
        WRITE(6,*) '*******************************'
        WRITE(6,'(5A12/4F12.9,I12)')
     $  ' th1w',' th2w',' th1n',' th2n','nphi',
     $    th1w ,  th2w ,  th1n ,  th2n , nphi
      ENDIF

! Beam energy
      ene = p1(4)
! Final electrons and photons not distinguished
      DO k=1,4
        pc(1,k)=p2(k)
        pc(2,k)=q2(k)
      ENDDO
      DO i=1,nphot
      DO k=1,4
        pc(2+i,k)=phot(i,k)
      ENDDO
      ENDDO
      np = nphot+2
      DO i=1,nmx1
      DO k=1,4
        pcw1(i,k)=0d0
        pcn1(i,k)=0d0
        pcw2(i,k)=0d0
        pcn2(i,k)=0d0       
      ENDDO
      ENDDO
!
! Collecting energies in calorimeter sectors
!
      DO I=1,NP
! wide/narrow sectors forward
      thet=angfi(pc(i,3),dsqrt(pc(i,1)**2+pc(i,2)**2))
      phi =angfi(pc(i,1),pc(i,2))
      langw= thet .GT. th1w .AND. thet .LT. th2w
      langn= thet .GT. th1n .AND. thet .LT. th2n
      DO jphi=1,nphi
        phi1 = (jphi-1)*(2d0*pi/nphi)
        phi2 =    jphi *(2d0*pi/nphi)
        lphi= phi .GT. phi1 .AND. phi .LT. phi2
        IF( langw .AND. lphi ) THEN
          DO k=1,4
            pcw1(jphi,k)=pcw1(jphi,k)+ pc(i,k)
          ENDDO
        ENDIF
        IF( langn .AND. lphi ) THEN
          DO k=1,4
            pcn1(jphi,k)=pcn1(jphi,k)+ pc(i,k)
          ENDDO
        ENDIF
      ENDDO
! wide/narrow sectors backward
      thet = angfi(-pc(i,3),dsqrt(pc(i,1)**2+pc(i,2)**2))
      phi  = angfi(-pc(i,1),-pc(i,2))
      langw= thet .GT. th1w .AND. thet .LT. th2w
      langn= thet .GT. th1n .AND. thet .LT. th2n
      DO jphi=1,nphi
        phi1 = (jphi-1)*(2d0*pi/nphi)
        phi2 =    jphi *(2d0*pi/nphi)
        lphi= phi .GT. phi1 .AND. phi .LT. phi2
        IF( langw .AND. lphi) THEN
          DO k=1,4
            pcw2(jphi,k)=pcw2(jphi,k)+ pc(i,k)
          ENDDO
        ENDIF
        IF( langn .AND. lphi) THEN
          DO k=1,4
            pcn2(jphi,k)=pcn2(jphi,k)+ pc(i,k)
          ENDDO
        ENDIF  
      ENDDO
      ENDDO
     
! at least one coincidences in a pair of opposite calorimetric blocks
      xwide= -3d0
      xnarr= -3d0
      xmix1= -3d0
      xmix2= -3d0
      yw1= 0d0
      yw2= 0d0
      yn1= 0d0
      yn2= 0d0
      DO jphi = 1,nphi
! minimum of energies deposed in the opposite blocks,
! we will require energies in opposite block above certain minimum
        xeww   = DMIN1(pcw1(jphi,4),pcw2(jphi,4))/ene
        xenn   = DMIN1(pcn1(jphi,4),pcn2(jphi,4))/ene
        xewn   = DMIN1(pcw1(jphi,4),pcn2(jphi,4))/ene
        xenw   = DMIN1(pcn1(jphi,4),pcw2(jphi,4))/ene
! maximum over pairs of blocks, 
! we will ask at least one pair of blocks excited above certain min.
        xwide = DMAX1(xwide,xeww)
        xnarr = DMAX1(xnarr,xenn)
        xmix1 = DMAX1(xmix1,xewn)
        xmix2 = DMAX1(xmix2,xenw)
        yw1 = yw1 +pcw1(jphi,4)/ene
        yw2 = yw2 +pcw2(jphi,4)/ene
        yn1 = yn1 +pcn1(jphi,4)/ene
        yn2 = yn2 +pcn2(jphi,4)/ene
      ENDDO
! max(x1,x2) type clasical choice as in trigas1 of TH.6118
      zcal(1) = xwide
      zcal(2) = xnarr
      zcal(3) = xmix1
      zcal(4) = xmix2
! No phi sectors and z-variable of the LL type
      zcal1(1) = yw1*yw2
      zcal1(2) = yn1*yn2
      zcal1(3) = yw1*yn2
      zcal1(4) = yn1*yw2

      END 

************************************************************************
************************************************************************
************************************************************************
