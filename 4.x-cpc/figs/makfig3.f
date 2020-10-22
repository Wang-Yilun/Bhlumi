      PROGRAM MAIN
!     ***********************************
! This is series of tests in which we vary dummy parameters in M.C.
! and we look whether distributions has changed.
! Note that all runs are from May95 with program with LL bug in Mat. Elm.
! this should not matter for these tests which are rather on crude level MC.
!---------------------
! To execute: make makfig3-dvi
!             make makfig3-ps
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib / b(50000)
      COMMON / INOUT  / NINP,NOUT
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!----------------------------------------------------------------------
      CALL GLIMIT(50000)
      NINP=  5
      NOUT= 16
      Tesnam    = 'fig-MC-params;'
      TeXfile   = 'fig-MC-params.tex'
      OPEN( NOUT, file='output-'//Tesnam)
      CALL GOUTPU(NOUT)
!--------------------------------------------------------
! Stored histograms and corresponding histograms
!--------------------------------------------------------
      lendan = 0
! Results with changed parameters
      lendan = lendan+1
!-------
! Deafult which may not work
      Dname(lendan)  = '../prod3/prod3.data'
      Hname(lendan)  = '../prod3/bhl.hst'
!-------
! epsilon changed
c>      TeXfile   = 'prod3.2083M.eps.e-5.tex'
c>      Dname(lendan)  = '../prod3/prod3.data.eps.e-5'
c>      Hname(lendan)  = '../prod3/prod3.hst.2083M.eps.e-5'
!-------
! keywgt=0
c>      TeXfile   = 'prod3.2090M.keywgt.eq.0.tex'
c>      Dname(lendan)  = '../prod3/prod3.data.keywgt.eq.0'
c>      Hname(lendan)  = '../prod3/prod3.hst.2090M.keywgt.eq.0'
! keywgt=1
c>      TeXfile   = 'prod3.1973M.keywgt.eq.1.tex'
c>      Dname(lendan)  = '../prod3/prod3.data.keywgt.eq.1'
c>      Hname(lendan)  = '../prod3/prod3.hst.1973M.keywgt.eq.1'
! keywgt=1
c>      TeXfile   = 'prod3.2028M.keywgt.eq.2.tex'
c>      Dname(lendan)  = '../prod3/prod3.data.keywgt.eq.2'
c>      Hname(lendan)  = '../prod3/prod3.hst.2028M.keywgt.eq.2'
!
! Reference results !!!  from May95 with bug !!!
! But bug for these tests is irrelevant
      lendan = lendan+1
      Dname(lendan)  = '../prod3/prod2.data.1995' ! ? recreated data
      Hname(lendan)  = '../prod3/prod2.hst.1995'  ! from May95 with bug

!--------------------------------------------------------
!----    initialize GPLOT
!--------------------------------------------------------
      CALL GPLINT(0)
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)

! =====================================
        CALL prod3
! =====================================

!--------------------------------------------------------
! ------------dumping histogram for control -------------------
!--------------------------------------------------------
      NOUTH=20
      OPEN(NOUTH,file='dump.hst')
      CALL GRFILE(NOUTH,DNAME,'N')
      CALL GROUT( 0,ICY,' ')
      CALL GREND(DNAME)
! ------------THE END OF HISTO WRITING -------------------------
      CALL GPLEND
      CLOSE(NOUT)
      END


      SUBROUTINE prod3
!     *****************
! BHLUMI-(OLDBIS+LUMLOG) and other corrections/components
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!----------------------------------------------------------------------
      SAVE
!----------------------------------------------------------------------
      COMMON / INOUT  / NINP,NOUT
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
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!----------------------------------------------------------------------
      WRITE(NOUT,*) '=============================================='
      WRITE(NOUT,*) '=============  prod3     ====================='
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
!==================================================================
!     Restore data and histo (2)
!==================================================================
      jene=2
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
!==================================================================
!==================================================================

! ==========================================================
! ==========================================================
! ==========================================================
! Introductory plots BLUM2, OLDBIS, LUMLOG separately
! and interesting differences for BLUM2
      CALL pfi3a
!
      END


      SUBROUTINE pfi3a
!     *****************
! Basic plots BLUM2,
!     *********************************
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
!
      WRITE(NOUT,BXL1I) KEYTRI,     'Type of Trigger    ','KEYTRI','  '
      WRITE(NOUT,BXL1F) BorW(1),    'Born Wide    [nb]  ','BORNW ','  '
      WRITE(NOUT,BXL1F) BorN(1),    'Born Narrow  [nb]  ','BORNN ','  '
! =========================================================
! --------- BHLUMI -------
      KeyGen = 3
      KeyGe2 = 3 +1000000
      JDA = 2000 +10000*KEYGEN
      JDU = 4000 +10000*KEYGEN
      KDA = 2000 +10000*KEYGEN +1000000
      KDU = 4000 +10000*KEYGEN +1000000
      LDA = 2000 +10000*KEYGEN +3000000
      LDU = 4000 +10000*KEYGEN +3000000
!
      DO itr=1,3
! O(alf2)expB real exp.
        CALL CUMHIS(KeyGen,JDA+200+itr,JDA+400+itr)
        CALL CUMHIS(KeyGen,JDA+300+itr,JDA+500+itr)
! O(alf2)expB, Theorists 1-Z
        CALL CUMHIS(KeyGen,JDA+360+itr,JDA+560+itr)
        CALL Gprint(JDA+300+itr)
!  Switching-off exponentiation,  Theorists 1-Z
!##     CALL CUMHIS(KeyGen,JDU+340+itr,JDU+540+itr)
      ENDDO
      DO itr=1,3
! O(alf2)expB real exp.
        CALL Gprint(KDA+200+itr)
        CALL CUMHIS(KeyGe2,KDA+200+itr,KDA+400+itr)
        CALL CUMHIS(KeyGe2,KDA+300+itr,KDA+500+itr)
! O(alf2)expB, Theorists
        CALL CUMHIS(KeyGe2,KDA+360+itr,KDA+560+itr)
        CALL Gprint(KDA+300+itr)
!  Switching-off exponentiation,  Theorists 1-Z
!##     CALL CUMHIS(KeyGen,KDU+340+itr,KDU+540+itr)
      ENDDO
      YMAX =  .0
      YMIN = -.15
! Energy cut versus Asymetricity
      CALL PL3MC('(200) OPSiW: Bhlum2eB, WW, NN, NW, 1-Zsum $'
     $,JDA+401,JDA+402,JDA+403,YMIN,YMAX,-BorW(1),-BorN(1),-BorN(1))
      CALL PL3MC('(300)SICAL: Bhlum2eB, WW, NN, WN, 1-Zmin $'
     $,JDA+501,JDA+502,JDA+503,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
      CALL PL3MC('(360)SICAL: Bhlum2eB, WW, NN, NW, 1-Z $'
     $,JDA+561,JDA+562,JDA+563,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
! Energy cut versus Asymetricity
      CALL PL3MC('(200) OPSiW: Bhlum2eB, WW, NN, NW, 1-Zsum $'
     $,KDA+401,KDA+402,KDA+403,YMIN,YMAX,-BorW(1),-BorN(1),-BorN(1))
      CALL PL3MC('(300)SICAL: Bhlum2eB, WW, NN, WN, 1-Zmin $'
     $,KDA+501,KDA+502,KDA+503,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
      CALL PL3MC('(360)SICAL: Bhlum2eB, WW, NN, NW, 1-Z $'
     $,KDA+561,KDA+562,KDA+563,YMIN,YMAX,-BorW(2),-BorN(2),-BorN(2))
! Take differnces
      DO itr=1,3
       CALL Gopera(JDA+400+itr,'-',KDA+400+itr,LDA+400+itr,1d0,1d0)
       CALL Gopera(JDA+500+itr,'-',KDA+500+itr,LDA+500+itr,1d0,1d0)
       CALL Gopera(JDA+560+itr,'-',KDA+560+itr,LDA+560+itr,1d0,1d0)
!##    CALL Gopera(JDU+540+itr,'-',KDU+540+itr,LDU+540+itr,1d0,1d0)
      ENDDO
      YMAX =  .0020
      YMIN = -.0020
! Plot differences
      CALL PL3MC('(200) OPSiW Diff: Bhlum2eB, WW, NN, NW, 1-Zsum $'
     $,LDA+401,LDA+402,LDA+403,YMIN,YMAX,BorW(1),BorN(1),BorN(1))
      CALL PL3MC('(300) SICAL Diff: Bhlum2eB, WW, NN, WN, 1-Zmin $'
     $,LDA+501,LDA+502,LDA+503,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
      CALL PL3MC('(360) SICAL Diff: Bhlum2eB, WW, NN, NW, 1-Z$'
     $,LDA+561,LDA+562,LDA+563,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
!  Switching-off exponentiation,  Theorists 1-Z
!##  CALL PL3MC('(340) SICAL, exponen. on/off , WW, NN, NW, 1-Z$'
!##  $,LDU+541,LDU+542,LDU+563,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
      DO itr=1,3
       CALL GDELET(JDA+400+itr)
       CALL GDELET(JDA+500+itr)
       CALL GDELET(JDA+560+itr)
       CALL GDELET(KDA+400+itr)
       CALL GDELET(KDA+500+itr)
       CALL GDELET(KDA+560+itr)
      ENDDO
!=========================================================
      END

