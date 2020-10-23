************************************************************************
************************************************************************
************************************************************************
! =====> trisic2  is  SICAL semirealistic ALEPH (Phys.lett. 95)
! =====> trisic2w is  SICAL2, simplified SICAL for workshop95
! =====> triosi2W is  semirealistic OPAL silicon
! See also trical.f where
! trical1  is  BERE1 CALO1 and CALO2
! trical2  is  simplified LCAL of TH.6118 (PLB268 1991)
************************************************************************
************************************************************************
************************************************************************
      SUBROUTINE TRISIC2W (Nphi,Nthe,Tmind,Tmaxd,Nseg,Npad,Z1,Z2,Z3)
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

      Ebeam = P1(4)
      IF(Icont.eq.0) THEN
        Icont=Icont+1
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
      CALL CLUSTE(TMIND,TMAXD,NPHI,NTHE,XCL1,XCL2,Npad,Nseg)
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
      END


      SUBROUTINE TRIOSIW(Nphi,Nthe,Tmind,Tmaxd,Z1,Z2,Z3,Z4)
!     *****************************************************************
! Input:  TMIND,TMAXD theta limits of the detector
!         NPAD        cluster polar     half width, default 16?
!         NSEG        cluster azimuthal half width, default  2
! Output: LWIDE,LNARR,LMIX1,LMIX2, each LOGICAL variable L parametrizes
!         acceptance condition in the corresponding N-N, W-W, N-W, W-N
! NPHI:  number of phi-sectors
! NTHE:  number of theta-sectors
!---------------------------------------------------------------------!
!
!--- Tailored for SiW by M.Dallavalle  Oct.94
!
! Detector size:
! Detector Z distance (mm) from I.P., minimum and maximum radial 
! coverage (mm)
!     PARAMETER (RMIN=62.008D0  , RMAX=142.008D0, ZDET= 2460.225D0 )
!
!   TMIND,TMAXD theta limits of the detector
!     TMIND = RMIN/ZDET   !=0.025204199
!     TMAXD = RMAX/ZDET   !=0.057721555
!
! NTHE:  number of detector theta-pads
! NPHI:  number of detector phi-sectors
!     NTHE = 32
!     NPHI = 32
! Cluster size:
!         NPAD        cluster polar     half width
!         NSEG        cluster azimuthal half width
!  NPAD varies in the range 3 to 16.
!  NSEG varies in the range 2 to 16.
!     NPAD = 16
!     NSEG =  2
! Energy of largest cluster on each side >= ECUT *Ebeam
! Sum of the energies of these clusters  >= SCUT *2*Ebeam
! Limits for Delta-phi cuts: PHIMIN < ABS(PHI1-PHI2) < PHIMAX
!---------------------------------------------------------------------!
!     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0 )
      SAVE
!---------------------------------------------------------------------!
!                     OPAL OSiW
!---------------------------------------------------------------------!
!CMD--   the following PARAMETER statements are for SiW
!**************** ALEPH for comparison************
!***      PARAMETER ( Ecut = 0.43865902, Scut =0.60315615 )
!***      PARAMETER ( Phimin = PI-0.52359878,  Phimax=PI+0.52359878)
!*************************************************
      PARAMETER ( Ecut = 0.5D0, Scut = 0.75D0 )
!**   PARAMETER ( PHIMIN = PI-0.2D0, PHIMAX=PI+0.2D0)
      PARAMETER ( Acopl = 0.2D0)
!**   PARAMETER ( Accut = 25.0D0/2460.225D0 )    !=0.010161672
!
      DIMENSION Z1(*),Z2(*),Z3(*),Z4(*)
!
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      DIMENSION XCL1(3),XCL2(3)
      LOGICAL Lthet(4),Lsolid
      LOGICAL LANGW1,LANGN1,LANGW2,LANGN2,Lacopl
      LOGICAL Lecut,Lscut,Lacoll(4)
      DATA Icont /0/

      Ebeam = P1(4)
      IF(Icont.eq.0) THEN
        Icont=Icont+1
        PAD = (TMAXD-TMIND)/NTHE
!*********** ALEPH for comparison ***********
!***        TH1W = TMIND + PAD
!***        TH2W = TMAXD - PAD
!***        TH1N = TMIND + 2*PAD
!***        TH2N = TMAXD - 4*PAD
!********************************************
        TH1W = TMIND + 2*PAD
        TH2W = TMAXD - 2*PAD
        TH1N = TMIND + 6*PAD
        TH2N = TMAXD - 6*PAD
        WRITE(6,*) '*******************************'
        WRITE(6,*) '***       OSiW  Trigger     ***'
        WRITE(6,*) '*******************************'
        WRITE(6,'(6A12/6F12.9)')
     $  'TMIND','TMAXD',' TH1W',' TH2W',' TH1N',' TH2N',
     $   TMIND , TMAXD ,  TH1W ,  TH2W ,  TH1N ,  TH2N
        WRITE(6,'(2A12/2I12)')
     $  'NPHI','NTHE',
     $   NPHI , NTHE
      ENDIF
! Cluster size
      Npad = 16
      Nseg = 2
! Clustering procedure
      CALL CLUSTE(TMIND,TMAXD,NPHI,NTHE,XCL1,XCL2,NPAD,NSEG)
      ECL1  =XCL1(1)
      ECL2  =XCL2(1)
! Energy cuts
      Lecut = (ECL1/Ebeam).GE.ECUT .AND. (ECL2/Ebeam).GE.ECUT 
      Lscut = ((ECL1+ECL2)/Ebeam).GE.(2*SCUT)
! Acoplanarity Delta-phi cut
      PHICL1=XCL1(2)
      PHICL2=XCL2(2)
      DELPHI=ABS(PHICL1-PHICL2)
      Lacopl= DELPHI.GT.(PI-Acopl) .AND. DELPHI.LT.(PI+Acopl)
! Acollinearity Delta-theta cut
      THECL1 =XCL1(3)
      THECL2 =XCL2(3)
      DELthe =ABS(THEcl1-THEcl2)
! Acollinearity cuts (radians)
      Accut1 = 0.005       ! Sharp 
      Accut2 = Pi          ! None
      Accut3 = 0.010161672 ! Standard
      Lacoll(1)= DELthe .LT. Accut1
      Lacoll(2)= DELthe .LT. Accut2
      Lacoll(3)= DELthe .LT. Accut3
      Lsolid = Lacopl .AND. Lacoll(3)
! Limits of wide and narrow fiducial regions
      LANGW1= THECL1.GT.TH1W.AND.THECL1.LT.TH2W
      LANGN1= THECL1.GT.TH1N.AND.THECL1.LT.TH2N
      LANGW2= THECL2.GT.TH1W.AND.THECL2.LT.TH2W
      LANGN2= THECL2.GT.TH1N.AND.THECL2.LT.TH2N
! Acceptance conditions for various triggers
      Lthet(1) = LangW1 .AND. LangW2 ! WW
      Lthet(2) = LangN1 .AND. LangN2 ! NN
      Lthet(3) = LangW1 .AND. LangN2 ! WN
      Lthet(4) = LangN1 .AND. LangW2 ! NW
!
! Preferred variable   z=s'/s
      Ztrue  = ecl1*ecl2/Ebeam**2
! OSiW basic cut variable
      Zsum   = (ecl1+ecl2)/(2*Ebeam)
! Not used
      Zmin   = MIN(ecl1,ecl2)/Ebeam
! Angular trigger included in all z's below
      DO i=1,4
        Z1(i)= 5
! Energy cut versus asymetricity
! True cut is for Zsum=Scut=0.75, 1-Zsum=0.25
        IF(Lthet(i) .AND. Lacopl .AND. Lacoll(3))  Z1(i) = 1-Zsum
! Energy cut versus acollinearity
        Z2(i)= 5
        IF(Lthet(3) .AND. Lacopl .AND. Lacoll(i))  Z2(i) = 1-Zsum
! Theorist choice for energy cut
        Z3(i)= 5
        IF(Lthet(i) .AND. Lsolid)  Z3(i) = 1-Ztrue
! Acollinearity cut versus Asymetricity
        Z4(i)= 5
        IF(Lthet(i) .AND. Lsolid)  Z4(i) = DELthe
      ENDDO
!
      END


      SUBROUTINE TRISIC2 (Nphi,Nthe,Tmind,Tmaxd,Z1,Z2,Z3)
!     *****************************************************************
!---------------------------------------------------------------------!
!     Version of TRISIC routine customized by S. Jadach 04.03.94      !
!     Developed by Wieslaw Placzek: 30.09.1993                        !
!     Last update by S. Jadach :    29.11.94                          !
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
!Short descroption of the trigger by Placzek.
!(programmed according to Bolek's description of the trigger).
!1) Angular range:
!   theta_min = 24 mrad (milirads)
!   theta_max = 58 mrad
!2) Division into pads  (calorimeter cells):
!   w theta (polar angle)   -   16 sektorow
!   w phi (azimuthal angle) -   32 sektory
! ----------------------------------------
!   altogether 16*32 pady
!3) trigger ranges in theta:
!   wide:   pady 2 - 15
!   narrow: pady 3 - 12
!4) Acceptance conditions:
!   E1_c, E2_c - cluster energies on both sides;
!   (i)    E1_c > 20 GeV,
!   (ii)   E2_c > 20 GeV,
!   (iii)  E1_c + E2_c > 55 GeV
!   (iv)   150 stopni < abs(phi_1 - phi_2) < 210 stopni
!   where phi_1, phi_2 - azimuthal angles of the clusters on both sides;
!Construction algorithm for cluster:
! - find pad with the biggest energy
! - add energies of the neighbour pads as fillows:
!   take Nseg pads on each side in the phi disrection
!   and Npad pads on each side in theta direction,
!   altogether  (1+2*Nseg)*(1+2*Npad) pads, except the case near 
!   the edge -- in this case we have less pads.
!   This algorith is repeated until all energy is exhausted
!   and we find several clusters on each side.
! - Total energies E_i of the clusters and 
!   the central positions (theta_i, phi_i) of the cluster 
!   (weighted with energy) is determined and additional cuts 
!   are applied on them.
!---------------------------------------------------------------------!
!     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0 )
      SAVE
!---------------------------------------------------------------------!
!                  ALEPH - SICAL
!Czesc,
!Jezeli bys chcial byc bardzo scisly to 
!  scut jest 0.6 (55 GeV) a ecut! 0.45 (20 Gev) 
! (tak bylo w analizie danych 92 roku. 
! Ruszanie cutami wyglada rozsadnie
!                            Powodzenia
!                                    Bolek
!---------------------------------------------------------------------!
! Energy of largest cluster on each side >= ECUT*Ebeam
! Sum of the energies of these clusters  >= SCUT*2*Ebeam
!**   PARAMETER ( Ecut = 20D0/91.187D0*2d0, Scut = 55D0/91.187D0 )
      PARAMETER ( Ecut = 0.43865902,        Scut = 0.60315615 )
! Limits for Delta-phi cuts: PHIMIN < ABS(PHI1-PHI2) < PHIMAX
! Here, Phimin,max = pi +- 0.52359878
!**   PARAMETER ( Phimin = 150D0*PI/180D0, Phimax=210D0*PI/180D0)
!**   PARAMETER ( Phimin = PI-0.52359878,  Phimax=PI+0.52359878)
      PARAMETER ( Acopl = 0.52359878 )
      DIMENSION Z1(*),Z2(*),Z3(*)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      DIMENSION XCL1(3),XCL2(3)
      LOGICAL Lthet(4), Lsize(4)
      LOGICAL LANGW1,LANGN1,LANGW2,LANGN2,Lacopl
      LOGICAL Lecut,Lscut
      DATA Icont /0/

      Ebeam = P1(4)
      IF(Icont.eq.0) THEN
        Icont=Icont+1
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
! Note important trick in this loop over cluster sizes:
! Standared choice of the cluster size is the LAST one k=3
      DO k=1,3
! Standard Cluster size
      IF(k .EQ. 1) THEN
        Npad = 0
        Nseg = 0
      ELSEIF(k .EQ. 2) THEN
        Npad = 8
        Nseg = 16
      ELSE
        Npad = 3
        Nseg = 1
      ENDIF
! Clustering procedure
      CALL CLUSTE(TMIND,TMAXD,NPHI,NTHE,XCL1,XCL2,Npad,Nseg)
      ECL1  =XCL1(1)
      ECL2  =XCL2(1)
! Energy cuts
      Lecut= ECL1/Ebeam.GE.Ecut .AND. ECL2/Ebeam.GE.Ecut 
      Lscut= (ECL1+ECL2)/(2*Ebeam).GE.Scut
! Delta-phi cut, Acoplanarity
      PHICL1=XCL1(2)
      PHICL2=XCL2(2)
      DELPHI=ABS(PHICL1-PHICL2)
      Lacopl= DELPHI.GT.(PI-Acopl) .AND. DELPHI.LT.(PI+Acopl)
! Theta acceptance, wide and narrow fiducial regions
      THECL1=XCL1(3)
      THECL2=XCL2(3)
      LANGW1= THECL1.GT.TH1W.AND.THECL1.LT.TH2W
      LANGN1= THECL1.GT.TH1N.AND.THECL1.LT.TH2N
      LANGW2= THECL2.GT.TH1W.AND.THECL2.LT.TH2W
      LANGN2= THECL2.GT.TH1N.AND.THECL2.LT.TH2N
! Theta and phi Acceptance  for various triggers          
      Lthet(1) = LANGW1.AND.LANGW2 ! WW
      Lthet(2) = LANGN1.AND.LANGN2 ! NN
      Lthet(3) = LANGW1.AND.LANGN2 ! WN
      Lthet(4) = LANGN1.AND.LANGW2 ! NW
! Size of the cluster dependence encoded in this variable
      Lsize(k) = Lthet(3) .AND.  Lacopl .AND.  Lscut
      ENDDO
!
!--- Final definitions
!
! Preferred by theorist variable   z=s'/s
      Ztrue  = ecl1*ecl2/Ebeam**2
! Additional SICAL cut on:
      Zsum   = (ecl1+ecl2)/(2*Ebeam)
! Basic SICAL cut on:
      Zmin   = MIN(ecl1,ecl2)/Ebeam
!---
      DO i=1,4
! Energy cut versus Asymetricity
! Aleph standard cut is for Zmin=Ecut=0.43, 1-Zmin=0.57
        Z1(i)= 5
        IF(Lthet(i) .AND. Lacopl .AND. Lscut)  Z1(i) = 1-Zmin
! Energy cut versus Cluster Size, 
! ALeph Standard size for  Lsize(3)
        Z2(i)= 5
        IF(Lsize(i))  Z2(i) = 1-Zmin
! Theoretically clean cut
        Z3(i)= 5
        IF(Lthet(i) .AND.  Lacopl)  Z3(i) = 1-Ztrue
      ENDDO
!
      END



      SUBROUTINE CLUSTE(THMIN,THMAX,NPHI,NTHE,XCL1,XCL2,NPAD,NSEG)
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
      DIMENSION PP(4,100),XDET(100,100),XCL1(3),XCL2(3)

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
      THET=ANGFI(PP(3,I),DSQRT(PP(1,I)**2+PP(2,I)**2))
      IF (THET.GT.THMIN .AND. THET.LT.THMAX) THEN
        PHI =ANGFI(PP(1,I),PP(2,I))
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
        CALL ONECLU(THMIN,THMAX,NPHI,NTHE,XDET,ENEC,PHIC,THEC,Npad,Nseg)
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
      THET=ANGFI(-PP(3,I),DSQRT(PP(1,I)**2+PP(2,I)**2))
      IF (THET.GT.THMIN .AND. THET.LT.THMAX) THEN
        PHI =ANGFI(-PP(1,I),-PP(2,I))
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
        CALL ONECLU(THMIN,THMAX,NPHI,NTHE,XDET,ENEC,PHIC,THEC,Npad,Nseg)
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

      SUBROUTINE 
     $ ONECLU(THMIN,THMAX,NPHI,NTHE,XDET,ENEC,PHIC,THEC,NPAD,NSEG)
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
      DIMENSION XDET(100,100)

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


