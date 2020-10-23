      PROGRAM MAIN
!     ***********************************
! To execute: make prod2-dvi
!     ***************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib / b(50000)
      COMMON / INOUT  / NINP,NOUT
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!---------------------------------------------------------------------- 
      Tesnam    = 'PlotProd2'
      TeXfile   = 'PlotProd2.tex'
      CALL GLIMIT(50000)
      NINP=  5
      NOUT= 16
      OPEN( NOUT, file='output-'//Tesnam)
      CALL GOUTPU(NOUT)
!--------------------------------------------------------
!----    initialize GPLOT
!--------------------------------------------------------
      CALL GPLINT(0)
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
! Stored histograms and corresponding histograms
!--------------------------------------------------------
      lendan = 0
! BHLUMI
      lendan = lendan+1
!>      Dname(lendan)  = '../prod2/prod2.data.2146M'  ! from 4.x-cpc
!>      Hname(lendan)  = '../prod2/prod2.hst.2146M'   ! from 4.x-cpc
!>      Dname(lendan)  = '../prod2/prod2.data.2037M'  ! April 96
!>      Hname(lendan)  = '../prod2/prod2.hst.2037M'   ! April 96
      Dname(lendan)  = '../prod2/prod2.data'        ! Current
      Hname(lendan)  = '../prod2/bhl.hst'           ! Current

!==========================================================
      CALL prod2
!==========================================================
! ------------dumping histogram for control -------------------
      NOUTH=20
      OPEN(NOUTH,file='dump.hst')
      CALL GRFILE(NOUTH,DNAME,'N')
      CALL GROUT( 0,ICY,' ')
      CALL GREND(DNAME)
! ------------THE END OF HISTO WRITING -------------------------
      CALL GPLEND
      CLOSE(NOUT)
      END



      SUBROUTINE prod2
C     *****************
! This is test BHLUMI-(OLDBIS+LUMLOG) for SICAL type triggers
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  

!----------------------------------------------------------------------
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G
      PARAMETER(
     $BXOPE =  '(//1X,15(5H=====)    )',
     $BXTXT =  '(1X,1H=,                  A48,25X,    1H=)',
     $BXL1I =  '(1X,1H=,I17,                 16X, A20,A12,A7, 1X,1H=)',
     $BXL1F =  '(1X,1H=,F17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2F =  '(1X,1H=,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXL1G =  '(1X,1H=,G17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2G =  '(1X,1H=,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXCLO =  '(1X,15(5H=====)/   )'    )
!----------------------------------------------------------------------
      SAVE
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!---------------------------------------------------------------------- 
      COMMON / INOUT  / NINP,NOUT     
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI


      WRITE(NOUT,*) '==============================================' 
      WRITE(NOUT,*) '=============  prod2     =====================' 
      WRITE(NOUT,*) '==============================================' 


!==================================================================
!     Restore data and histo (1)
!==================================================================
      jene=1
      WRITE(   6,BXOPE)
      WRITE(   6,BXTXT) Dname(jene)
      WRITE(   6,BXTXT) Hname(jene)
!-- Reading input-data file used for MC
      CALL ReaDat(Dname(jene))
      WRITE(   6,BXL1F) cmsene, 'total CMS energy  ','<<<---','=='
!-- Restore histograms from j-th  directory
      NINPH=10
      OPEN(NINPH,file=Hname(jene))
      CALL GRFILE(NINPH,' ',' ')
      CALL GRIN(0,0,0)
      CLOSE(NINPH)
!-- End of restoring

C==========================================================
C==========================================================
C==========================================================
! Introductory plots BLUM2, OLDBIS, LUMLOG separately
! and interesting differences for BLUM2
      CALL FGCOR1
!
C==========================================================
C==========================================================
C==========================================================
      END



      SUBROUTINE FGCOR1
C     *****************
C     Basic plots BLUMI
C     ********************************* 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI     =  3.1415926535897932D0) 
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G 
      PARAMETER(
     $BXOPE =  '(//1X,15(5H=====)    )',
     $BXTXT =  '(1X,1H=,                  A48,25X,    1H=)',
     $BXL1I =  '(1X,1H=,I17,                 16X, A20,A12,A7, 1X,1H=)',
     $BXL1F =  '(1X,1H=,F17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2F =  '(1X,1H=,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXL1G =  '(1X,1H=,G17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2G =  '(1X,1H=,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXCLO =  '(1X,15(5H=====)/   )'    )   
      SAVE
      COMMON / INOUT  / NINP,NOUT
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI
      DIMENSION BorW(10),BorN(10)
      LOGICAL gexist


! Born for OPSiW, wide and narrow
      THosi1 = 0.025204199d0
      THosi2 = 0.057721555d0
      Nthe = 32
      PAD = (THosi2-THosi1)/Nthe
      BorW(1)  = BORNB(CMSENE,THosi1 +2*PAD, THosi2 -2*PAD)
      BorN(1)  = BORNB(CMSENE,THosi1 +6*PAD, THosi2 -6*PAD)
! Born for SICAL, wide and narrow
      THsic1 = .024d0
      THsic2 = .058d0
      Nthe = 16
      PAD = (THsic2-THsic1)/Nthe
      BorW(2)  = BORNB(CMSENE,THsic1   +PAD, THsic2   -PAD)
      BorN(2)  = BORNB(CMSENE,THsic1 +2*PAD, THsic2 -4*PAD)
! Born for BARE1, wide and narrow
      THbar1 = 0.024d0
      THbar2 = 0.058d0
      Nthe = 16
      PAD = (THbar2-THbar1)/Nthe
      BorW(3)  = BORNB(CMSENE,THbar1       , THbar2       )
      BorN(3)  = BORNB(CMSENE,THbar1 +1*PAD, THbar2 -1*PAD)
! Born for SICAL2 and CALO2, wide and narrow
      THbar1 = 0.024d0
      THbar2 = 0.058d0
      Nthe = 16
      PAD = (THbar2-THbar1)/Nthe
      BorW(4)  = BORNB(CMSENE,THbar1 +1*PAD, THbar2 -1*PAD)
      BorN(4)  = BORNB(CMSENE,THbar1 +2*PAD, THbar2 -4*PAD)
!
      WRITE(NOUT,BXL1I) KEYTRI,     'Type of Trigger    ','KEYTRI','  '
      WRITE(NOUT,BXL1F) BorW(4),    'Born Wide    [nb]  ','BORNW ','  '
      WRITE(NOUT,BXL1F) BorN(4),    'Born Narrow  [nb]  ','BORNN ','  '
!=========================================================
C --------- BHLUMI -------
      KeyGen = 3
      JDA = 2000 +10000*KEYGEN
      kda = JDA +10000000
!
      DO itr=1,3
! O(alf2)expB real exp.
        CALL CUMHIS(KeyGen,JDA+200+itr,kda+200+itr)
        CALL CUMHIS(KeyGen,JDA+210+itr,kda+210+itr)
        CALL CUMHIS(KeyGen,JDA+220+itr,kda+220+itr)
        CALL CUMHIS(KeyGen,JDA+230+itr,kda+230+itr)
!
        CALL CUMHIS(KeyGen,JDA+300+itr,kda+300+itr)
        CALL CUMHIS(KeyGen,JDA+310+itr,kda+310+itr)
        CALL CUMHIS(KeyGen,JDA+320+itr,kda+320+itr)
!
        CALL CUMHIS(KeyGen,JDA+340+itr,kda+340+itr)
        CALL CUMHIS(KeyGen,JDA+350+itr,kda+350+itr)
        CALL CUMHIS(KeyGen,JDA+360+itr,kda+360+itr)
!
        CALL CUMHIS(KeyGen,JDA+370+itr,kda+370+itr)
        CALL CUMHIS(KeyGen,JDA+380+itr,kda+380+itr)
        CALL CUMHIS(KeyGen,JDA+390+itr,kda+390+itr)
      ENDDO
      YMAX =  .0
      YMIN = -.15
! Energy cut versus Asymetricity
      CALL PL3MC('[200] OPSiW: Bhlum2eB, WW, NN, NW, 1-Umin $' 
     $,kda+201,kda+202,kda+203,YMIN,YMAX,-BorW(1),-BorN(1),-BorN(1))
! THeorist Energy cut versus  Asymetricity
      CALL PL3MC('[210] OPSiW: Bhlum2eB, WW, NN, NW, 1-Zmin $' 
     $,kda+211,kda+212,kda+213,YMIN,YMAX,-BorW(1),-BorN(1),-BorN(1))
! Energy cut versus Acollinearity
      CALL PL3MC('[220] OPSiW: Bhlum2eB, Acol. Sharp,None,Stnd,1-Zsum$' 
     $,kda+221,kda+222,kda+223,YMIN,YMAX,-BorN(1),-BorN(1),-BorN(1))
! Acollinearity versus  Asymetricity
      CALL PL3MC('[230] OPSiW: Bhlum2eB, WW, NN, NW, 1-Zsum$' 
     $,kda+231,kda+232,kda+233,YMIN,YMAX,-BorW(1),-BorN(1),-BorN(1))

      YMAX =  .0
      YMIN = -.15
! Energy cut versus Asymetricity
      CALL PL3MC('[300] SICAL: Bhlum2eB, WW, NN, WN, 1-Umin $' 
     $,kda+301,kda+302,kda+303,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
! THeorist Energy cut versus  Asymetricity
      CALL PL3MC('[310] SICAL: Bhlum2eB, WW, NN, NW, 1-Zmin $' 
     $,kda+311,kda+312,kda+313,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
! Energy cut versus cluster size
      YMAX = -.03
      YMIN = -.07
      CALL PL3MC('[320] SICAL: 2eB; Clust. Smal, Huge, Stnd, 1-Zmin$' 
     $,kda+321,kda+322,kda+323,YMIN,YMAX,-BorN(2),-BorN(2),-BorN(2))

      YMAX =  .0
      YMIN = -.15
! Energy cut versus Asymetricity
      CALL PL3MC('[340] BARE1: Bhlum2eB, WW, NN, WN, 1-Zmin $' 
     $,kda+341,kda+342,kda+343,YMIN,YMAX,-BorW(3),-BorN(3),-BorN(3))
      CALL PL3MC('[350] CALO2: Bhlum2eB, WW, NN, WN, 1-Zmin $' 
     $,kda+351,kda+352,kda+353,YMIN,YMAX,-BorW(4),-BorN(4),-BorN(4))
      CALL PL3MC('[360] SICAL2: Bhlum2eB, WW, NN, WN, 1-Zmin $' 
     $,kda+361,kda+362,kda+363,YMIN,YMAX,-BorW(4),-BorN(4),-BorN(4))

! Energy cut versus Asymetricity
      CALL PL3MC('[370]  SICAL: 2eB FPZ, WW, NN, WN, 1-Zmin $' 
     $,kda+371,kda+372,kda+373,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
      CALL PL3MC('[380]  CALO2: 2eB FPZ, WW, NN, WN, 1-Zmin $' 
     $,kda+381,kda+382,kda+383,YMIN,YMAX,-BorW(4),-BorN(4),-BorN(4))
      CALL PL3MC('[390] SICAL2: 2eB FPZ, WW, NN, WN, 1-Zmin $' 
     $,kda+391,kda+392,kda+393,YMIN,YMAX,-BorW(4),-BorN(4),-BorN(4))


*===================================================
*               Differences
*===================================================
      JDU = 4000 +10000*KEYGEN
      kdu = JDU +10000000
      DO itr=1,3
! Differences
        CALL CUMHIS(KeyGen,JDU+300+itr,kdu+300+itr)
        CALL CUMHIS(KeyGen,JDU+200+itr,kdu+200+itr)

        CALL CUMHIS(KeyGen,JDU+260+itr,kdu+260+itr)
        CALL CUMHIS(KeyGen,JDU+270+itr,kdu+270+itr)
        CALL CUMHIS(KeyGen,JDU+400+itr,kdu+400+itr)

        CALL CUMHIS(KeyGen,JDU+310+itr,kdu+310+itr)
        CALL CUMHIS(KeyGen,JDU+320+itr,kdu+320+itr)
        CALL CUMHIS(KeyGen,JDU+330+itr,kdu+330+itr)
        CALL CUMHIS(KeyGen,JDU+340+itr,kdu+340+itr)

        CALL CUMHIS(KeyGen,JDU+350+itr,kdu+350+itr)
        CALL CUMHIS(KeyGen,JDU+360+itr,kdu+360+itr)
        CALL CUMHIS(KeyGen,JDU+370+itr,kdu+370+itr)
        CALL CUMHIS(KeyGen,JDU+380+itr,kdu+380+itr)

        CALL CUMHIS(KeyGen,JDU+410+itr,kdu+410+itr)
        CALL CUMHIS(KeyGen,JDU+420+itr,kdu+420+itr)
        CALL CUMHIS(KeyGen,JDU+430+itr,kdu+430+itr)
        CALL CUMHIS(KeyGen,JDU+440+itr,kdu+440+itr)

        CALL CUMHIS(KeyGen,JDU+600+itr,kdu+600+itr)
        CALL CUMHIS(KeyGen,JDU+610+itr,kdu+610+itr)
        CALL CUMHIS(KeyGen,JDU+700+itr,kdu+700+itr)
        CALL CUMHIS(KeyGen,JDU+710+itr,kdu+710+itr)
      ENDDO

      YMAX =  .003
      YMIN = -.003
! BHLUMI  4.03 - 2.02  difference [Bhlumi4-Bhlumi2 'Correction']
      CALL PL3MC('[300] SICAL 2eB-1eA: CORR. PL 95, WW,NN,NW, 1-Zmin$' 
     $,kdu+301,kdu+302,kdu+303,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
      CALL PL3MC('[200] OPSiW 2eB-1eA: CORRECTION,  WW,NN,NW, 1-Zsum$' 
     $,kdu+201,kdu+202,kdu+203,YMIN,YMAX,BorW(1),BorN(1),BorN(1))

! Specials on CUT-OFF variations for
! O(alf2)eB, Energy cut versus Acollinearity
      CALL PL3MC('[260] OPSiW 2eB-1eA: CORR. 1-Zsum: Acol=5,1000,10$' 
     $,kdu+261,kdu+262,kdu+263,YMIN,YMAX,BorW(1),BorN(1),BorN(1))
! O(alf2)eB, Acollinearity versus Asimetricity
      CALL PL3MC('[270] OPSiW 2eB-1eA: CORR. Acol, for WW,NN,NW$' 
     $,kdu+271,kdu+272,kdu+273,YMIN,YMAX,BorW(1),BorN(1),BorN(1))
! O(alf2)eB, Energy cut versus Cluster size
      CALL PL3MC('[400] SICAL 2eA-1eA: CORR.1-Zmin: ClSiz=1x1,all,7x3$' 
     $,kdu+401,kdu+402,kdu+403,YMIN,YMAX,BorW(2),BorN(2),BorN(2))

! Switching-off exponentiation
! O(alf2)eB-O(alf2)eA, Energy cut versus Asimetricity
      YMAX =  .003
      YMIN = -.003
      CALL PL3MC('[310] SICAL 2eB-2B: expon.on/off PL95: 1-Z:WW,NN,NW$' 
     $,kdu+311,kdu+312,kdu+313,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
      CALL PL3MC('[320] SICAL2 2eB-2B: expon.on/off: 1-Z:WW,NN,NW$' 
     $,kdu+321,kdu+322,kdu+323,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
      CALL PL3MC('[320] CALO2  2eB-2B: expon.on/off: 1-Z:WW,NN,NW$' 
     $,kdu+331,kdu+332,kdu+333,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
      CALL PL3MC('[340] CALO2  2eB-2B: expon.on/off: 1-Z:WW,NN,NW$' 
     $,kdu+341,kdu+342,kdu+343,YMIN,YMAX,BorW(3),BorN(3),BorN(3))

! BHLUMI  4.03-2.02 difference  [Bhlumi4-Bhlumi2 'Correction']
      CALL PL3MC('[350] SICAL 2eB-1eA: CORR. PL95, WW, NN, NW, 1-Z$' 
     $,kdu+351,kdu+352,kdu+353,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
      CALL PL3MC('[360] SICAL2 2eB-1eA: CORRECTION, WW, NN, NW, 1-Z$' 
     $,kdu+361,kdu+362,kdu+363,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
      CALL PL3MC('[370]  CALO2 2eB-1eA: CORRECTION, WW, NN, NW, 1-Z$' 
     $,kdu+371,kdu+372,kdu+373,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
      CALL PL3MC('[380]  BARE1 2eB-1eA: CORRECTION, WW, NN, NW, 1-Z$' 
     $,kdu+381,kdu+382,kdu+383,YMIN,YMAX,BorW(3),BorN(3),BorN(3))

!  (B)-(A)  EXP   Second Order      O(alf2)eB-O(alf2)eA
      CALL PL3MC('[410] SICAL 2eB-2eA: New-Old: 1-Z: WW,NN,NW$' 
     $,kdu+411,kdu+412,kdu+413,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
!  (B)-(A)  UNEXP Second Order      O(alf2)B-O(alf2)A
      CALL PL3MC('[420] SICAL 2B-2A: New-Old: 1-Z: WW,NN,NW$' 
     $,kdu+421,kdu+422,kdu+423,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
!  (B) Second-First Order  EXP      O(alf2)eB-O(alf1)eB
      CALL PL3MC('[430] SICAL 2eB-1eB:  WW, NN, NW, 1-Z$' 
     $,kdu+431,kdu+432,kdu+433,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
!  (B) Second-First Order  UNEXP    O(alf2)B -O(alf1)B
      CALL PL3MC('[440] SICAL 2B-1B:  WW, NN, NW, 1-Z$' 
     $,kdu+441,kdu+442,kdu+443,YMIN,YMAX,BorW(4),BorN(4),BorN(4))

!=========================================================
! LL-bug correction
      YMAX =  .003
      YMIN = -.003
      IF(gexist(kdu+600+1)) THEN
      CALL PL3MC('[600] OPSiW LL-bug CORRECTION, WW, NN, NW, 1-Zsum$' 
     $,kdu+601,kdu+602,kdu+603,YMIN,YMAX,BorW(1),BorN(1),BorN(1))
      CALL PL3MC('[700] SICAL LL-bug CORRECTION, WW, NN, NW, 1-Zmin$' 
     $,kdu+701,kdu+702,kdu+703,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
!
      CALL PL3MC('[610] OPSiW LL-bug CORRECTION, WW, NN, NW, 1-Z$' 
     $,kdu+611,kdu+612,kdu+613,YMIN,YMAX,BorW(1),BorN(1),BorN(1))
      CALL PL3MC('[710] SICAL LL-bug CORRECTION, WW, NN, NW, 1-Z$' 
     $,kdu+711,kdu+712,kdu+713,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
      ENDIF

!=========================================================
      END
