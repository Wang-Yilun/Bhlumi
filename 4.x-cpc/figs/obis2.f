      PROGRAM MAIN
!     ***********************************
! To execute: make obis2-plot
!     ***************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib / b(50000)
      COMMON / INOUT  / NINP,NOUT
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!---------------------------------------------------------------------- 
      Tesnam    = 'obis2'
      TeXfile   = 'obis2.tex'
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
! OLDBIS
      lendan = lendan+1
!>      Dname(lendan)  = '../obis2/obis2.data.2118M'  ! April 96
!>      Hname(lendan)  = '../obis2/obis2.hst.2118M'   ! April 96
      Dname(lendan)  = '../obis2/obis2.data'    ! Current
      Hname(lendan)  = '../obis2/bhl.hst'       ! Current

!==========================================================
      CALL obis2
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



      SUBROUTINE obis2
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
      WRITE(NOUT,*) '=============  plofi2    =====================' 
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
! Introductory plots OLDBIS separately
      CALL FGCOR1
C==========================================================
C==========================================================
C==========================================================
      END



      SUBROUTINE FGCOR1
C     *****************
C Basic plots BLUM2, OLDBIS, LUMLOG separately
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
! Born for Sical, wide and narrow
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
      WRITE(NOUT,BXL1F) BorW(1),    'Born Wide    [nb]  ','BORNW ','  ' 
      WRITE(NOUT,BXL1F) BorN(1),    'Born Narrow  [nb]  ','BORNN ','  ' 
!
!=========================================================
!=========================================================
C --------- OLDBIS -------
      KeyGen = 1
C O(alf)exp correction pure bremss.
      JDA = 2000 +10000*KEYGEN
      kda = JDA +10000000
      DO itr=1,3
        CALL CUMHIS(KeyGen,JDA+200+itr,kda+200+itr)
        CALL CUMHIS(KeyGen,JDA+210+itr,kda+210+itr)
        CALL CUMHIS(KeyGen,JDA+220+itr,kda+220+itr)
        CALL CUMHIS(KeyGen,JDA+230+itr,kda+230+itr)

        CALL CUMHIS(KeyGen,JDA+300+itr,kda+300+itr)
        CALL CUMHIS(KeyGen,JDA+310+itr,kda+310+itr)
        CALL CUMHIS(KeyGen,JDA+320+itr,kda+320+itr)

        CALL CUMHIS(KeyGen,JDA+340+itr,kda+340+itr)
        CALL CUMHIS(KeyGen,JDA+350+itr,kda+350+itr)
        CALL CUMHIS(KeyGen,JDA+360+itr,kda+360+itr)
      ENDDO
!
      YMAX =  .0
      YMIN = -.15
! Energy cut versus Asymetricity
      CALL PL3MC('[200] OPSiW: OLDbis, WW, NN, NW, 1-Zsum $' 
     $,kda+201,kda+202,kda+203,YMIN,YMAX,-BorW(1),-BorN(1),-BorN(1))
! THeorist Energy cut versus  Asymetricity
      CALL PL3MC('[210] OPSiW: OLDbis, WW, NN, NW, 1-Z $' 
     $,kda+211,kda+212,kda+213,YMIN,YMAX,-BorW(1),-BorN(1),-BorN(1))


! Energy cut versus Asymetricity
      CALL PL3MC('[300] SICAL: OLDbis, WW, NN, WN, 1-Zmin $' 
     $,kda+301,kda+302,kda+303,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
! THeorist Energy cut versus  Asymetricity
      CALL PL3MC('[310] SICAL:  OLDbis, WW, NN, NW, 1-Z $' 
     $,kda+311,kda+312,kda+313,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
! Energy cut versus cluster size
      YMAX = -.04
      YMIN = -.08
      CALL PL3MC('[320] SICAL: OLDBIS, Clust.Smal,Huge,Stnd, 1-Zmin$' 
     $,kda+321,kda+322,kda+323,YMIN,YMAX,-BorN(2),-BorN(2),-BorN(2))

! Energy cut versus Acollinearity
      YMAX =  .0
      YMIN = -.15
      CALL PL3MC('[220] OPSiW: OLDBIS, Acol.Sharp,None,Stnd, 1-Zsum$' 
     $,kda+221,kda+222,kda+223,YMIN,YMAX,-BorN(1),-BorN(1),-BorN(1))
! Acollinearity versus  Asymetricity
      CALL PL3MC('[230] OPSiW: OLDBIS, WW, NN, NW, , Acol. [rad] $' 
     $,kda+231,kda+232,kda+233,YMIN,YMAX,-BorW(1),-BorN(1),-BorN(1))

! Energy cut versus  Asymetricity
      CALL PL3MC('[340] BARE1:  OLDbis, WW, NN, NW, 1-Z $' 
     $,kda+341,kda+312,kda+343,YMIN,YMAX,-BorW(3),-BorN(3),-BorN(3))
      CALL PL3MC('[350] CALO2:  OLDbis, WW, NN, NW, 1-Z $' 
     $,kda+351,kda+352,kda+353,YMIN,YMAX,-BorW(4),-BorN(4),-BorN(4))
      CALL PL3MC('[360] SICAL2:  OLDbis, WW, NN, NW, 1-Z $' 
     $,kda+361,kda+362,kda+363,YMIN,YMAX,-BorW(4),-BorN(4),-BorN(4))
!=========================================================
      END
