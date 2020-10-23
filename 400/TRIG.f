************************************************************************
************************************************************************
*****    Collection of LEP lumi event selection subrograms (triggers)
************************************************************************
!        1990-93, see also trical.f (even older trigas.f)
! =====> trical1 is  BERE1 CALO1 and CALO2, Figs.12-13 in workshop95
! =====> trical2 is  simplifies LCAL of TH.6118 (PLB268 1991)
************************************************************************
!        1994-96, see also silicon.f
! =====> trisic2  is  SICAL semirealistic ALEPH (Phys.lett. 95)
! =====> trisic2w is  SICAL2 simplified SICAL, Fig.14 in workshop95
! =====> triosi2W is  semirealistic OPAL silicon
************************************************************************
************************************************************************
!               !!! NEW!!! NEW!!! NEW!!! NEW!!!
!     Separate namespace for all subroutines in TRIG.f
!     and common block in TRIG.h is achieved using prefix TRIG_
!     (In this way TRIG.* package can be safely linked within bigger f77 code.)
!     TRIG_UniSical is configurable trigger for silicon LEP detectors.
!     All cut-off parameters reside TRIG.h can be set/reset using setter routines.
!     Setters are used in c++ wrappers.
!     TRIG_UniSical can be also used to emulate triggers for FCAL, BGO and SAT detectors.
!     TRIG_UniElcal configurable version of trigas1 of TH.6118,
!     and trical2 of workshop96. (So far not used much.)
!     TRIG_TRIGAS2, TRIG_TRISIC2W, TRIG_trical1 and TRIG_trical2
!     are left in original form for backward compatibility.
!     Version of TRIG_TRIGAS2 is slightly modified by Patrick Janot
!     ANGFI function is incorporated as TRIG_ANGFI.
************************************************************************
************************************************************************

      SUBROUTINE TRIG_UniSical(lWW,lNN,lNW,lWN)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Universal trigger covering all LEP cases
! Output lWW,lNN,lNW,lWN are 0,1 for rejected/accepted event (c++ convention)
      IMPLICIT NONE
      SAVE
      INTEGER lWW,lNN,lNW,lWN
!
      INCLUDE 'TRIG.h'
!
      INTEGER NPHOT
      DOUBLE PRECISION P1,Q1,P2,Q2,PHOT
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT

      DOUBLE PRECISION  XCL1(3),XCL2(3)           ! results of clustering
      DOUBLE PRECISION  EneF,EneB,PhiF,PhiB,TheF,TheB,DELPHI,EBEAM
!
      DOUBLE PRECISION  xAngFW,xAngFN,xAngBW,xAngBN,xEcutMin,xEcutMax,xEcutSum,xEcutProd,xAcopl,xAcoli
      INTEGER           lCut
!
      DOUBLE PRECISION  PI
      PARAMETER( PI = 3.1415926535897932D0 )
!
      INTEGER icont
      DATA icont /0/
!
      icont=icont+1
      IF(icont.eq.1) THEN
        WRITE(6,*) '*******************************'
        WRITE(6,*) '***  TRIG_UniSical Trigger   ***'
        WRITE(6,*) '*******************************'
! Defasult values as in SICAL of ALEPH
        TMIND = 0.024d0
        TMAXD = 0.058d0
        NPHI = 32  ! detector granularity
        NTHE = 16  ! detector granularity
        Npad = 16  ! clustering
        Nseg = 2   ! clustering
        WRITE(6,'(6A12/2F12.9,4I12)')
     $  ' TMIND',' TMAXD','nphi','nthe','npad','nseg',
     $    TMIND ,  TMIND , nphi , nthe , npad , nseg
! Full fidutial angular range
        th1w = TMIND ! wide minimum
        th2w = TMAXD ! wide maximum
        th1n = TMIND ! narrow minimum
        th2n = TMAXD ! narrow maximum
! Initialization of cutoff parameters to "no cuts"
        EcutMin = 0d0
        EcutMax = 0d0
        EcutSum = 0d0
        EcutProd= 0d0
        Acopl   = 1d9
        Acoli   = 1d9
      ENDIF

! Clustering procedure
      CALL TRIG_CLUSTE(TMIND,TMAXD,NPHI,NTHE,XCL1,XCL2,NPAD,NSEG)
! Theta acceptance, wide and narrow fiducial regions
      TheF =XCL1(3)
      TheB =XCL2(3)
      xAngFW = (th1W-TheF)*(TheF-th2W)
      xAngFN = (th1N-TheF)*(TheF-th2N)
      xAngBW = (th1W-TheB)*(TheB-th2W)
      xAngBN = (th1N-TheB)*(TheB-th2N)

! Delta-phi cut, Acoplanarity
      PhiF =XCL1(2)
      PhiB =XCL2(2)
      DELPHI=ABS(PhiF-PhiB)
ccc   Lacopl= DELPHI.GT.(PI-Acopl) .AND. DELPHI.LT.(PI+Acopl)
      xAcopl = ((PI-Acopl) -DELPHI)*(DELPHI -(PI+Acopl))

! Acolinearity cut, in theta diff.
      xAcoli  = Acoli -ABS(TheF-TheB)

! three kinds of energy cuts
ccc   Lecut= ECL1/Ebeam.GE.Ecut .AND. ECL2/Ebeam.GE.Ecut
ccc   Lscut= (ECL1+ECL2)/(2*Ebeam).GE.Scut
      EneF =XCL1(1)
      EneB =XCL2(1)
      EBEAM=P1(4)
      xEcutMin  = MIN(EneF,EneB)/EBEAM - EcutMin
      xEcutMax  = MAX(EneF,EneB)/EBEAM - EcutMax
      xEcutSum  = (EneF+EneB)/(2*EBEAM)- EcutSum
      xEcutProd = EneF*EneB/EBEAM/EBEAM- EcutProd

      lCut =0
      IF( xAcopl.GT.0 .AND. xAcoli.GT.0 .AND. xEcutMin.GT.0 .AND. xEcutMax.GT.0
     $    .AND. xEcutSum.GT.0 .AND. xEcutProd.GT.0 ) lCut=1
!!!      IF( xEcutProd.GT.0 ) lCut=1
      lWW =0
      lNN =0
      lNW =0
      lWN =0
      IF((xAngFW.GT.0).AND.(xAngBW.GT.0).AND.(lCut.EQ.1)) lWW =1  !WW
      IF((xAngFN.GT.0).AND.(xAngBN.GT.0).AND.(lCut.EQ.1)) lNN =1  !NN
      IF((xAngFW.GT.0).AND.(xAngBN.GT.0).AND.(lCut.EQ.1)) lWN =1  !WN
      IF((xAngFN.GT.0).AND.(xAngBW.GT.0).AND.(lCut.EQ.1)) lNW =1  !NW
!
      END ! of TRIG_UniSical

      SUBROUTINE TRIG_SetGrid(yTMIND,yTMAXD,yNPHI,yNTHE)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION     yTMIND,yTMAXD
      INTEGER              yNPHI,yNTHE
!
      INCLUDE 'TRIG.h'
!
      TMIND = yTMIND ! Fidutial minimum angle
      TMAXD = yTMAXD ! Fidutial maximum angle
      NPHI  = yNPHI  ! No of sectors in phi
      NTHE  = yNTHE  ! No of sectors in theta
      END ! TRIG_SetGrid


      SUBROUTINE TRIG_SetClust(yNpad,yNseg)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      SAVE
      INTEGER              yNpad,yNseg
!
      INCLUDE 'TRIG.h'
!
!     Cluster size params
      Npad = yNpad ! (+-3) No. of pads in theta taken to make a cluster
      Nseg = yNseg ! (+-1) No. of segments in phi taken to make a cluster
      END ! TRIG_SetClust


      SUBROUTINE TRIG_SetAng(yth1w,yth2w, yth1n,yth2n)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION yth1w,yth2w, yth1n,yth2n
!
      INCLUDE 'TRIG.h'

      th1w = yth1w
      th2w = yth2w
      th1n = yth1n
      th2n = yth2n
      END ! of TRIG_SetAng

      SUBROUTINE TRIG_SetAcol(yAcopl,yAcoli)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION yAcopl,yAcoli
!
      INCLUDE 'TRIG.h'
!
      Acopl = yAcopl
      Acoli = yAcoli

      END ! TRIG_SetAcol

      SUBROUTINE TRIG_SetEcut(yEcutMin,yEcutMax,yEcutSum,yEcutProd)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION    yEcutMin,yEcutMax,yEcutSum,yEcutProd
!
      INCLUDE 'TRIG.h'
!
      EcutMin = yEcutMin
      EcutMax = yEcutMax
      EcutSum = yEcutSum
      EcutProd= yEcutProd
!
      END ! TRIG_SetEcut


      SUBROUTINE TRIG_UniElCal(lWW,lNN,lNW,lWN)
!     **********************************************************
! LCAL Idealized exper. CALORIMETRIC trigger on dressed final electrons.
! Formely trigas1 of TH.6118, and trical2 of workshop96
! Electrons and photons not distinguished!
! This was used initially in paper TH.6118.
!     ******************************************
      IMPLICIT NONE
      SAVE
      INTEGER lWW,lNN,lNW,lWN
!
      INCLUDE 'TRIG.h'
!
      INTEGER NPHOT
      DOUBLE PRECISION P1,Q1,P2,Q2,PHOT
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT

      INTEGER nmx1
      PARAMETER( nmx1 = 20)
      DOUBLE PRECISION pc(100,4)
      DOUBLE PRECISION pcw1(nmx1,4),pcn1(nmx1,4)
      DOUBLE PRECISION pcw2(nmx1,4),pcn2(nmx1,4)

      DOUBLE PRECISION phi, phi1, phi2, thet, xenn, xenw, xewn, xeww, Ebeam
      DOUBLE PRECISION xnarr, xwide, yn1, yn2, yw1, yw2, xmix1, xmix2
      DOUBLE PRECISION TRIG_angfi
      LOGICAL langw,langn,lphi

      DOUBLE PRECISION pi
      PARAMETER( pi = 3.1415926535897932D0 )

      INTEGER jphi, k, np, i, icont
      DATA icont /0/

      IF(icont.eq.0) THEN
        icont=icont+1
        nphi = 6
        th1w = 0.043d0  ! ALEPH 1990-92
        th2w = 0.125d0
        th1n = 0.057d0
        th2n = 0.107d0
        EcutMin = 0.44d0  ! ALEPH 1990-92 ??
        WRITE(6,*) '*******************************'
        WRITE(6,*) '*** TRIG_UniElCal Trigger   ***'
        WRITE(6,*) '*******************************'
        WRITE(6,'(5A12/4F12.9,I12)')
     $  ' th1w',' th2w',' th1n',' th2n','nphi',
     $    th1w ,  th2w ,  th1n ,  th2n , nphi
      ENDIF

! Beam energy
ccc      ene = p1(4)
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
      thet=TRIG_ANGFI(pc(i,3),dsqrt(pc(i,1)**2+pc(i,2)**2))
      phi =TRIG_ANGFI(pc(i,1),pc(i,2))
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
      thet = TRIG_ANGFI(-pc(i,3),dsqrt(pc(i,1)**2+pc(i,2)**2))
      phi  = TRIG_ANGFI(-pc(i,1),-pc(i,2))
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
        xeww   = DMIN1(pcw1(jphi,4),pcw2(jphi,4))
        xenn   = DMIN1(pcn1(jphi,4),pcn2(jphi,4))
        xewn   = DMIN1(pcw1(jphi,4),pcn2(jphi,4))
        xenw   = DMIN1(pcn1(jphi,4),pcw2(jphi,4))
! maximum over pairs of blocks,
! we will ask at least one pair of blocks excited above certain min.
        xwide = DMAX1(xwide,xeww)
        xnarr = DMAX1(xnarr,xenn)
        xmix1 = DMAX1(xmix1,xewn)
        xmix2 = DMAX1(xmix2,xenw)
        yw1 = yw1 +pcw1(jphi,4)
        yw2 = yw2 +pcw2(jphi,4)
        yn1 = yn1 +pcn1(jphi,4)
        yn2 = yn2 +pcn2(jphi,4)
      ENDDO

      lWW =0
      lNN =0
      lNW =0
      lWN =0
      Ebeam = P1(4)
      IF( xwide/Ebeam .GT. EcutMin) lWW =1  !WW
      IF( xnarr/Ebeam .GT. EcutMin) lNN =1  !NN
      IF( xmix1/Ebeam .GT. EcutMin) lWN =1  !WN
      IF( xmix2/Ebeam .GT. EcutMin) lNW =1  !NW

      END


************************************************************************
************************************************************************
************************************************************************
************************************************************************
************************************************************************


      SUBROUTINE TRIG_trical1(th1n,th2n,th1w,th2w,delcon,dlthe,dlphi,zcalA,zcalB,zcalC)
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
      thetp=TRIG_ANGFI(p2(3),dsqrt(p2(1)**2+p2(2)**2))
      phip =TRIG_ANGFI(p2(1),p2(2))
      lforw= thetp .GT. th1w .AND. thetp .LT. th2w
      lforn= thetp .GT. th1n .AND. thetp .LT. th2n
! wide/narrow sectors backward
      thetq=TRIG_ANGFI(-q2(3),dsqrt(q2(1)**2+q2(2)**2))
      phiq =TRIG_ANGFI(-q2(1),-q2(2))
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
        thet1=TRIG_ANGFI( phot(k,3),dsqrt(phot(k,1)**2+phot(k,2)**2))
        phi1 =TRIG_ANGFI( phot(k,1),phot(k,2))
        difphi1 = min(ABS(phi1 -phip),2*pi-ABS(phi1 -phip))
        IF( abs(thet1-thetp) .LT. dlthe .AND.
     $               difphi1 .LT. dlphi)  e1=e1+phot(k,4)
! backward
        thet2=TRIG_ANGFI(-phot(k,3),dsqrt(phot(k,1)**2+phot(k,2)**2))
        phi2 =TRIG_ANGFI(-phot(k,1),-phot(k,2))
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


      SUBROUTINE TRIG_trical2(nphi,th1n,th2n,th1w,th2w,zcal,zcal1)
!     **********************************************************
! LCAL Idealized exper. CALORIMETRIC trigger on dressed final electrons.
! Formely trigas1 of TH.6118
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
      thet=TRIG_ANGFI(pc(i,3),dsqrt(pc(i,1)**2+pc(i,2)**2))
      phi =TRIG_ANGFI(pc(i,1),pc(i,2))
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
      thet = TRIG_ANGFI(-pc(i,3),dsqrt(pc(i,1)**2+pc(i,2)**2))
      phi  = TRIG_ANGFI(-pc(i,1),-pc(i,2))
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
      SUBROUTINE TRIG_TRISIC2W (Nphi,Nthe,Tmind,Tmaxd,Nseg,Npad,Z1,Z2,Z3)
!     *****************************************************************
!---------------------------------------------------------------------!
!     Version of trisic routine for SICAL2 trigger of WORKSHOP 95     !
!     Simplified, customized by S. Jadach 03.04.95                    !
!     Developed by Wieslaw Placzek: 30.09.1993                        !
!     Last update by S. Jadach :    03.04.95                          !
!---------------------------------------------------------------------!
! Input:
!    COMMON / MOMSET / Four momenta and number of photons
!    TMIND,TMAXD theta limits of the detector (radians)
!    NPHI:  total number of phi-sectors
!    NTHE:  total number of theta-sectors
!    NPAD,  (=3) Number of pads in theta taken to make a cluster
!    NSEG,  (=1) Number of segments in phi taken to make a cluster
! Output: Z1, Z2, Z3 energy cut parameters
!         all of them implicitly include angular cuts!!!!
!---------------------------------------------------------------------!
!     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0 )
      SAVE
!---------------------------------------------------------------------!
      DIMENSION Z1(*),Z2(*),Z3(*)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      DIMENSION XCL1(3),XCL2(3)
      LOGICAL Lthet(4)
      LOGICAL langw1,langn1,langw2,langn2
!
      INTEGER icont
      DATA icont /0/

      Ebeam = P1(4)
      Icont=Icont+1
      IF(Icont.eq.1) THEN
        PAD = (TMAXD-TMIND)/NTHE
        TH1W = TMIND + PAD
        TH2W = TMAXD - PAD
        TH1N = TMIND + 2*PAD
        TH2N = TMAXD - 4*PAD
        WRITE(6,*) '*******************************'
        WRITE(6,*) '***       SICAL Trigger     ***'
        WRITE(6,*) '*******************************'
        WRITE(6,'(6A12/6F12.9)')
     $  'TMIND','TMAXD',' TH1W',' TH2W',' TH1N',' TH2N',
     $   TMIND , TMAXD ,  TH1W ,  TH2W ,  TH1N ,  TH2N
        WRITE(6,'(2A12/2I12)')
     $  'NPHI','NTHE',
     $   NPHI , NTHE
      ENDIF
!
! Clustering procedure
      CALL TRIG_CLUSTE(TMIND,TMAXD,NPHI,NTHE,XCL1,XCL2,Npad,Nseg)
      ECL1  =XCL1(1)
      ECL2  =XCL2(1)
! Theta acceptance, wide and narrow fiducial regions
      THECL1=XCL1(3)
      THECL2=XCL2(3)
      LANGW1= THECL1 .GT. TH1W .AND. THECL1 .LT. TH2W
      LANGN1= THECL1 .GT. TH1N .AND. THECL1 .LT. TH2N
      LANGW2= THECL2 .GT. TH1W .AND. THECL2 .LT. TH2W
      LANGN2= THECL2 .GT. TH1N .AND. THECL2 .LT. TH2N
! Theta and phi Acceptance  for various triggers
      Lthet(1) = LANGW1 .AND. LANGW2 ! WW
      Lthet(2) = LANGN1 .AND. LANGN2 ! NN
      Lthet(3) = LANGW1 .AND. LANGN2 ! WN
      Lthet(4) = LANGN1 .AND. LANGW2 ! NW
!
!--- Final definitions
!
! Preferred by theorist variable   z=s'/s
      Ebeam=P1(4)
      Ztrue  = ecl1*ecl2/Ebeam**2
!---
      DO i=1,4
        Z1(i)= -5
        IF(Lthet(i))  Z1(i) = Ztrue
      ENDDO
!
      END ! TRISIC2W


      SUBROUTINE TRIG_CLUSTE(THMIN,THMAX,NPHI,NTHE,XCL1,XCL2,NPAD,NSEG)
!     ************************************************************
! Idealized procedure for constructing calorimetric clusters.
! Electrons and photons not distinguished!
! It provides exerimental values of energy and angles for largest
! clusters in both sides of the detector.
! Input:  THMIN, THMAX,  min. and max. polar angles of the trigger
!         NPHI,       number of azimuthal sections of the trigger (<101)
!         NTHE,       number of polar     sections of the trigger (<101)
!         NPAD,  Number of pads in theta taken to make a cluster
!         NSEG,  Number of segments in phi taken to make a cluster
! Output: XCL1(1),  energy of largest cluster on the e+ side
!         XCL1(2),  azimuthal angle of largest cluster on the e+ side
!         XCL1(3),  polar angle of largest cluster on the e+ side
!         XCL2(K), (K=1,3) as the above for the e- side
! NOTE:   Azimuthal angles are calculated for +z axis along e+ beam!
!         While polar angles are in the range (THMIN,THMAX), i.e. they
!         are evaluated for each detector according to corresponding beam.
! --- Developed by Wieslaw Placzek: Aug./Sept. 1993
! --- Last update:                  21.09.1993
!     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0 )
      PARAMETER( DWAPI = 2D0*PI )
! NTMX: max. number of clusters
      PARAMETER( NCMX = 10 )
      SAVE
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      DIMENSION PP(4,100),XDET(200,200),XCL1(3),XCL2(3)
!
      DO 6 J=1,3
      XCL1(J) = 0
      XCL2(J) = 0
    6 CONTINUE
      ENC1=0D0
      PHC1=0D0
      THC1=0D0
      ENC2=0D0
      PHC2=0D0
      THC2=0D0

      DO 10 J=1,NPHI
      DO 10 I=1,NTHE
   10 XDET(I,J)=0D0
! All final particle into one array
      DO 30 K=1,4
      PP(K,1)=P2(K)
   30 PP(K,2)=Q2(K)
      DO 35 I=1,NPHOT
      DO 35 K=1,4
   35 PP(K,2+I)=PHOT(I,K)
      NP=2+NPHOT
! Fill calorimetric cells on the e+ side
      EDE1=0D0
      DO 40 I=1,NP
      IF( (pp(3,i)**2+pp(1,i)**2+pp(2,i)**2) .GT. 0d0) THEN ! avoid 1/zero
      THET=TRIG_ANGFI(PP(3,I),DSQRT(PP(1,I)**2+PP(2,I)**2))
      IF (THET.GT.THMIN .AND. THET.LT.THMAX) THEN
        PHI =TRIG_ANGFI(PP(1,I),PP(2,I))
        JPHI=PHI/DWAPI*DBLE(NPHI) + 1D0
        JTHE=(THET-THMIN)/(THMAX-THMIN)*DBLE(NTHE) + 1D0
        XDET(JTHE,JPHI)=XDET(JTHE,JPHI)+PP(4,I)
        EDE1=EDE1+PP(4,I)
      ENDIF
      ENDIF
   40 CONTINUE
! Make clusters on the e+ side
      IF (EDE1.GT.0D0) THEN
        ESUM1=0D0
        ENC1=0D0
        DO 50 IC=1,NCMX
        CALL TRIG_ONECLU(THMIN,THMAX,NPHI,NTHE,XDET,ENEC,PHIC,THEC,Npad,Nseg)
        IF (ENEC.GT.ENC1) THEN
          ENC1=ENEC
          PHC1=PHIC
          THC1=THEC
        ENDIF
        ESUM1=ESUM1+ENEC
        IF (ESUM1.GE.0.99D0*EDE1) GOTO 51
   50   CONTINUE
   51   CONTINUE
      ELSE
        ENC1=0D0
        PHC1=0D0
        THC1=0D0
      ENDIF
! Fill calorimetric cells on the e- side
      DO 20 J=1,NPHI
      DO 20 I=1,NTHE
   20 XDET(I,J)=0D0
      EDE2=0D0
      DO 60 I=1,NP
      IF( (pp(3,i)**2+pp(1,i)**2+pp(2,i)**2) .GT. 0d0) THEN ! avoid 1/zero
      THET=TRIG_ANGFI(-PP(3,I),DSQRT(PP(1,I)**2+PP(2,I)**2))
      IF (THET.GT.THMIN .AND. THET.LT.THMAX) THEN
        PHI =TRIG_ANGFI(-PP(1,I),-PP(2,I))
        JPHI=0.5D0*PHI/PI*DBLE(NPHI) + 1D0
        JTHE=(THET-THMIN)/(THMAX-THMIN)*DBLE(NTHE) + 1D0
        XDET(JTHE,JPHI)=XDET(JTHE,JPHI)+PP(4,I)
        EDE2=EDE2+PP(4,I)
      ENDIF
      ENDIF
   60 CONTINUE
! Make clusters on the e- side
      IF (EDE2.GT.0D0) THEN
        ESUM2=0D0
        ENC2=0D0
        DO 70 IC=1,NCMX
        CALL TRIG_ONECLU(THMIN,THMAX,NPHI,NTHE,XDET,ENEC,PHIC,THEC,Npad,Nseg)
        IF (ENEC.GT.ENC2) THEN
          ENC2=ENEC
          PHC2=PHIC
          THC2=THEC
        ENDIF
        ESUM2=ESUM2+ENEC
        IF (ESUM2.GE.0.99D0*EDE2) GOTO 71
   70   CONTINUE
   71   CONTINUE
      ELSE
        ENC2=0D0
        PHC2=0D0
        THC2=0D0
      ENDIF
      XCL1(1)=ENC1
      XCL1(2)=PHC1
      XCL1(3)=THC1
      XCL2(1)=ENC2
      XCL2(2)=DMOD(PHC2+PI,DWAPI)
      XCL2(3)=THC2
      END

      SUBROUTINE TRIG_ONECLU(THMIN,THMAX,NPHI,NTHE,XDET,ENEC,PHIC,THEC,NPAD,NSEG)
!     ************************************************************
! Idealized procedure for constructing one calorimetric cluster
! Electrons and photons not distinguished!
! It provides exerimental values of energy and angles of one cluster.
! Input:  THMIN, THMAX,  min. and max. polar angles of the trigger
!         NPHI,       number of azimuthal sections of the trigger
!         NTHE,       number of polar     sections of the trigger
!         XDET(NTHE,NPHI),  array of filled calorimetric cells
! Output: ENEC,       energy of cluster
!         PHIC,       azimuthal angle of cluster
!         THEC,       polar     angle of cluster
! --- Developed by Wieslaw Placzek: Aug./Sept. 1993
! --- Last update:                  30.09.1993
!     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0 )
      PARAMETER( DWAPI = 2D0*PI )
      SAVE
      DIMENSION XDET(200,200)

! Find a cell with maximum energy
      EX=0D0
      DO 10 J=1,NPHI
      DO 10 I=1,NTHE
      EC=XDET(I,J)
      IF (EC.GT.EX) THEN
        EX=EC
        IX=I
        JX=J
      ENDIF
   10 CONTINUE
      IF (EX.GT.0D0) THEN
        IMIN=IX-NPAD
        IMAX=IX+NPAD
        IF (IMIN.LT.1) IMIN=1
        IF (IMAX.GT.NTHE) IMAX=NTHE
        ENEC=0D0
        EPHI=0D0
        ETHE=0D0
        PHCE=2D0*PI/DBLE(NPHI)
        THCE=(THMAX-THMIN)/DBLE(NTHE)
! Make a cluster
        DO 20 J=JX-Nseg,JX+Nseg
!!!!+++        DO 20 J=JX-1,JX+1
        DO 20 IC=IMIN,IMAX
        JC=J
        IF (JC.LT.1) JC=NPHI+J
        IF (JC.GT.NPHI) JC=JC-NPHI
        ECEL=XDET(IC,JC)
        ENEC=ENEC+ECEL
        EPHI=EPHI+ECEL*(DBLE(J)-0.5D0)*PHCE
        ETHE=ETHE+ECEL*(THMIN+(DBLE(IC)-0.5D0)*THCE)
        XDET(IC,JC)=0D0
   20   CONTINUE
        PHIC=EPHI/ENEC
        IF (PHIC.LT.0D0) THEN
          PHIC=DWAPI+PHIC
        ELSEIF (PHIC.GT.DWAPI) THEN
          PHIC=PHIC-DWAPI
        ENDIF
        THEC=ETHE/ENEC
      ELSE
        ENEC=0D0
        PHIC=0D0
        THEC=0D0
      ENDIF
! ==== End of clustering package ============================
      END

****************************************************************
*                   ALEPH LCAL
* Narrow range approximately   3.3-6.3 deg.
* Energy default cut
*     min(Ecl_1+Ecl_2)/E_beam  > 0.44 and
*     (Ecl_1+Ecl_2)/(2*E_beam) > 0.6
****************************************************************
* Calling sequence:
*     CALL TRIGAS2(0.057D0,0.107D0,0.043D0,0.125D0,16,XWIDE,XNARR,XMIX1,XMIX2)
*     IF ( XMIX1 .LT. 0.6D0 ) GOTO 999
****************************************************************
      SUBROUTINE TRIG_TRIGAS2(TH1N,TH2N,TH1W,TH2W,NPHI,XWIDE,XNARR,XMIX1,XMIX2)
C     **********************************************************
C MODIFIED TRIGAS1 FOR ALEPH GEOMETRICAL CUTS
C Idealized exper. CALORIMETRIC trigger on dressed final electrons.
C Electrons and photons not distinguished!
C     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0 )
      PARAMETER( DWAPI = 2D0*PI )
      LOGICAL LANGW,LANGN,LPHI
      LOGICAL LNARROW,LWIDE,TRIG_KEEPIT
      LOGICAL LAWIDE, LANARR, LAMIX1, LAMIX2
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      REAL *4 X4,Y4
      DIMENSION PC(100,4)
      DIMENSION PCW1(4),PCN1(4)
      DIMENSION PCW2(4),PCN2(4)
      DATA ICONT /0/

      IF(icont.eq.0) THEN
        icont=icont+1
        WRITE(6,*) '*******************************'
        WRITE(6,*) '***     ELCAL  Trigger     ***'
        WRITE(6,*) '*******************************'
        WRITE(6,'(5A12/4F12.9,I12)')
     $  ' TH1W',' TH2W',' TH1N',' TH2N','NPHI',
     $    TH1W ,  TH2W ,  TH1N ,  TH2N , NPHI
      ENDIF
C Beam energy
      ENE = P1(4)
C Final electrons and photons not distinguished
      DO 10 K=1,4
      PC(1,K)=P2(K)
   10 PC(2,K)=Q2(K)
      DO 20 I=1,NPHOT
      DO 20 K=1,4
   20 PC(2+I,K)=PHOT(I,K)
      NP = NPHOT+2
      DO 40 K=1,4
      PCW1(K)=0D0
      PCN1(K)=0D0
      PCW2(K)=0D0
   40 PCN2(K)=0D0
C
C Collecting energies in calorimeter sectors
C
      Z0 = 2800.D0
      DO 100 I=1,NP
C wide/narrow sectors forward
C Staszek's angels
      THETA=TRIG_ANGFI(ABS(PC(I,3)), DSQRT(PC(I,1)**2+PC(I,2)**2))
      PHI  =TRIG_ANGFI(PC(I,1),PC(I,2))

C Bolek's angels
      X     = Z0*PC(I,1)/ABS(PC(I,3))
      Y     = Z0*PC(I,2)/ABS(PC(I,3))

      LNARROW = .FALSE.
      LWIDE   = .TRUE.
      IF( TRIG_KEEPIT(X,Y,5) ) LNARROW = .TRUE.
      CALL TRIG_LOOSEC(THETA,PHI,LWIDE)
      LANGW = .NOT. LWIDE
      LANGN = LNARROW
C**      X4 = X
C**      Y4 = Y
C**      IF(LANGN) CALL HF2(100,X4,Y4,1.)
C**      IF(LANGW) CALL HF2(200,X4,Y4,1.)
C**      TYPE *,I,PC(I,1),PC(I,2),THETA,THET,PHIA,PHI,X,Y,
C**     &       LNARROW,LWIDE,langw,LANGN

      IF(PC(I,3) .GT. 0D0) THEN

      IF(LANGW) THEN
        DO 50 K=1,4
   50   PCW1(K)=PCW1(K)+ PC(I,K)
      ENDIF
      IF(LANGN) THEN
        DO 51 K=1,4
   51   PCN1(K)=PCN1(K)+ PC(I,K)
      ENDIF

      ELSE

C wide/narrow sectors backward
      IF(LANGW) THEN
        DO 70 K=1,4
   70   PCW2(K)=PCW2(K)+ PC(I,K)
      ENDIF
      IF(LANGN) THEN
        DO 71 K=1,4
   71   PCN2(K)=PCN2(K)+ PC(I,K)
      ENDIF

      ENDIF
C**      TYPE *,PCW1(4),PCW2(4),PCN1(4),PCN2(4)
  100 CONTINUE

C at least one coincidences in a pair of opposite calorimetric blocks
      XWIDE= 0D0
      XNARR= 0D0
      XMIX1= 0D0
      XMIX2= 0D0
c[[[[[[[[[ Replaced by S.J.
c      IF(PCW1(4)/ENE .GT. 0.44D0  .AND.
c     &   PCW2(4)/ENE .GT. 0.44D0  .AND.
c     &  (PCW1(4)+PCW2(4))/(2D0*ENE) .GT. 0.6D0) XWIDE = PCW1(4)+PCW2(4)
c      IF(PCN1(4)/ENE .GT. 0.44D0  .AND.
c     &   PCN2(4)/ENE .GT. 0.44D0  .AND.
c     &  (PCN1(4)+PCN2(4))/(2D0*ENE) .GT. 0.6D0) XNARR = PCN1(4)+PCN2(4)
c      IF(PCW1(4)/ENE .GT. 0.44D0  .AND.
c     &   PCN2(4)/ENE .GT. 0.44D0  .AND.
c     &  (PCW1(4)+PCN2(4))/(2D0*ENE) .GT. 0.6D0) XMIX1 = PCW1(4)+PCN2(4)
c      IF(PCN1(4)/ENE .GT. 0.44D0  .AND.
c     &   PCW2(4)/ENE .GT. 0.44D0  .AND.
c     &  (PCN1(4)+PCW2(4))/(2D0*ENE) .GT. 0.6D0) XMIX2 = PCN1(4)+PCW2(4)
C
C Added by PJ
C
      PHI1W = TRIG_ANGFI(PCW1(1),PCW1(2))
      PHI2W = TRIG_ANGFI(PCW2(1),PCW2(2))
      PHI1N = TRIG_ANGFI(PCN1(1),PCN1(2))
      PHI2N = TRIG_ANGFI(PCN2(1),PCN2(2))
      DPWIDE =ABS(PHI2W-PHI1W)
      DPNARR =ABS(PHI2N-PHI1N)
      DPMIX1 =ABS(PHI2N-PHI1W)
      DPMIX2 =ABS(PHI2W-PHI1N)
      Acopl = 0.17453293D0
      LAWIDE = DPWIDE.GT.(PI-Acopl) .AND. DPWIDE.LT.(PI+Acopl)
      LANARR = DPNARR.GT.(PI-Acopl) .AND. DPNARR.LT.(PI+Acopl)
      LAMIX1 = DPMIX1.GT.(PI-Acopl) .AND. DPMIX1.LT.(PI+Acopl)
      LAMIX2 = DPMIX2.GT.(PI-Acopl) .AND. DPMIX2.LT.(PI+Acopl)
!
! Replaced by S.J.
!
      IF(PCW1(4)/ENE .GT. 0.44D0  .AND. PCW2(4)/ENE .GT. 0.44D0 .AND.
     $   LAWIDE )  XWIDE = (PCW1(4)+PCW2(4))/(2*Ene)
      IF(PCN1(4)/ENE .GT. 0.44D0  .AND. PCN2(4)/ENE .GT. 0.44D0 .AND.
     &   LANARR )  XNARR = (PCN1(4)+PCN2(4))/(2*Ene)
      IF(PCW1(4)/ENE .GT. 0.44D0  .AND. PCN2(4)/ENE .GT. 0.44D0 .AND.
     $   LAMIX1 )  XMIX1 = (PCW1(4)+PCN2(4))/(2*Ene)
      IF(PCN1(4)/ENE .GT. 0.44D0  .AND. PCW2(4)/ENE .GT. 0.44D0 .AND.
     $   LAMIX2 )  XMIX2 = (PCN1(4)+PCW2(4))/(2*Ene)
c
      END

C***************** ADDED BY BOLEK 20 - DEC - 1990 *****************

      LOGICAL FUNCTION TRIG_KEEPIT(X,Y,METH)
C----------------------------------------------------------------------
C!  -
C!
C!   Author   :- John Renner Hansen    23-JAN-1990
C!
C!   Inputs:Method number:  METHOD = (5,6)
C!          X,Y of particle
C!        -
C!
C!   Outputs: KEEPIT
C!        -
C!
C!   Libraries required:
C!
C!   Description
C!   ===========
C!
C?
      IMPLICIT REAL*8(A-H,O-Z)
      REAL *8 PADS/29.75D0/
      REAL *8 X,Y,XP,YP
      REAL *8 LXP5(9)/4.5D0,4.5D0,3.5D0,2.5D0,5*1.5D0/
      REAL *8 HXP(9)/10.5D0,9.5D0,9.5D0,9.5D0,8.5D0,7.5D0
     &              ,6.5D0,5.5D0,3.5D0/
      REAL *8 LYP(10)/0D0,3.D0,4.D0,5.D0,6.D0,7.D0
     &               ,8.D0,9.D0,10.D0,11.D0/
      REAL *8 LXP6(9)/5.5D0,5.5D0,4.5D0,2.5D0,5*1.5D0/
      INTEGER*4 METH
      LOGICAL LOG1,LOG2,LOG3,LOG4
C!======================================================================
      TRIG_KEEPIT = .FALSE.
C
C     DEFINE ACTIVE REGION
      XP  =  ABS(X/PADS)
      YP  =  ABS(Y/PADS)
C
      IF(METH.EQ.5) THEN
      DO I = 1,9
        LOG1 = XP.GT.LXP5(I)
        LOG2 = XP.LT.HXP(I)
        LOG3 = YP.GE.LYP(I)
        LOG4 = YP.LE.LYP(I+1)
        IF(LOG1.AND.LOG2.AND.LOG3.AND.LOG4) THEN
          TRIG_KEEPIT = .TRUE.
          GOTO 999
        ENDIF
      ENDDO
      ELSE
      DO I = 1,9
        LOG1 = XP.GT.LXP6(I)
        LOG2 = XP.LT.HXP(I)
        LOG3 = YP.GT.LYP(I)
        LOG4 = YP.LT.LYP(I+1)
        IF(LOG1.AND.LOG2.AND.LOG3.AND.LOG4) THEN
          TRIG_KEEPIT = .TRUE.
          GOTO 999
        ENDIF
      ENDDO
      ENDIF
  999 RETURN
      END

c********************************************************************
      SUBROUTINE TRIG_LOOSEC(THETA,PHI,REJFL)
C----------------------------------------------------------------------
C!  - Bhabha selection method 5 on loose side
C!
C!   Author   :- Peter H. Hansen       27-NOV-1990
C!
C!   Inputs: THETA,PHI of cluster centroid in global system
C!   Output: REJFL = .TRUE. if rejected
C!                 = .FALSE. if accepted
C?
C!======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DXN(4),DYN(4),DZN(4),OMXN(4),OMYN(4),OMZN(4)
      DIMENSION DIST(4)
      LOGICAL REJFL
C-----------------------------------------------------------------------
      REJFL = .TRUE.
*      CALL VZERO(DXN,2*4)
*      CALL VZERO(DYN,2*4)
*      CALL VZERO(DZN,2*4)
*      CALL VZERO(OMXN,2*4)
*      CALL VZERO(OMYN,2*4)
*      CALL VZERO(OMZN,2*4)
*      CALL VZERO(DIST,2*2)
      DO k=1,4
      DXN(k) =0d0
      DYN(k) =0d0
      DZN(k) =0d0
      OMXN(k)=0d0
      OMYN(k)=0d0
      OMZN(k)=0d0
      DIST(k)=0d0
      ENDDO
C
C Sign z and module number
      IF(THETA.GT.1.D0) THEN
        SIG  =-1.D0
        MODU = 2
      ELSE
        SIG  = 1.D0
        MODU = 4
      ENDIF
      IF(COS(PHI).LT.0.D0) MODU=MODU-1
C
C X and Y in global system
C
      X = SIG*280.D0*COS(PHI)*TAN(THETA)
      Y = SIG*280.D0*SIN(PHI)*TAN(THETA)
C
C Local system
      ZL = SIG*17.5D0
      XLOC = X - DXN(MODU) - OMZN(MODU)*Y + OMYN(MODU)*ZL
      YLOC = Y - DYN(MODU) - OMXN(MODU)*ZL + OMZN(MODU)*X
C
C Fold into first quadrant
      XLOC = ABS(XLOC)
      YLOC = ABS(YLOC)
C
C Find distance to edges (as in LCLUTW, but database hardwired)
      DIST(1) = XLOC - 1.9D0
      DIST(2) = 24.5D0
      IF(YLOC.GT.8.4D0+1.2D0) THEN
        DIS = XLOC-(17.5D0-YLOC)/0.75D0
        IF(DIS.LT.DIST(1)) DIST(1)=DIS
        DIST(2) = YLOC-(17.5D0-XLOC*0.75D0)
      ELSE
        DIST(1) = XLOC-11.9D0
      ENDIF
C
C Now make the cut on the inner boundary
C and on the outer boundary
      IF(THETA.GT.1.D0) THETA = 3.14159D0-THETA
      IF(DIST(1).GT.1.D0  .AND.
     &   DIST(2).GT.0.75D0  .AND.
     &   THETA.LT.0.125D0)  REJFL = .FALSE.
C**      TYPE *,DIST(1),DIST(2),THETA,REJFL
C
  999 CONTINUE
      END



      FUNCTION TRIG_ANGFI(X,Y)
C     *******************
* CALCULATES ANGLE IN (0,2*PI) RANGE OUT OF X-Y
*     ***********************
      IMPLICIT REAL*8(A-H,O-Z)
      DATA PI /3.1415926535897932D0/

      IF(ABS(Y).LT.ABS(X)) THEN
        THE=ATAN(ABS(Y/X))
        IF(X.LE.0D0) THE=PI-THE
      ELSE
        THE=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      IF(Y.LT.0D0) THE=2D0*PI-THE
      ANGFI=THE
      END




