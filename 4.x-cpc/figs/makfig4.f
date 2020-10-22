      PROGRAM MAIN
! ***********************************
! To execute: make makfig4-dvi
!             make makfig4-ps
! ***************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib / b(50000)
      COMMON / INOUT  / NINP,NOUT
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!---------------------------------------------------------------------- 
      Tesnam    = 'fig-MisO3'
      TeXfile   = 'fig-MisO3.tex'
      CALL GLIMIT(50000)
      NINP=  5
      NOUT= 16
      OPEN( NOUT, file='output-'//Tesnam)
      CALL GOUTPU(NOUT)
!--------------------------------------------------------
! Stored histograms and corresponding histograms
!--------------------------------------------------------
      lendan = 0
! LUMLOG
      lendan = lendan+1
!>      Dname(lendan)  = '../llog2/llog2.data.537M'  ! July96
!>      Hname(lendan)  = '../llog2/llog2.hst.537M'   ! July96
      Dname(lendan)  = '../llog2/llog2.data'     ! Current
      Hname(lendan)  = '../llog2/bhl.hst'        ! Current


!==========================================================
      CALL llog2
!==========================================================
! ------------dumping histogram for control -------------------
      NOUTH=20
      OPEN(NOUTH,file='dump.hst')
      CALL GRFILE(NOUTH,DNAME,'N')
      CALL GROUT( 0,ICY,' ')
      CALL GREND(DNAME)

      CLOSE(NOUT)
      END



      SUBROUTINE llog2
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
! Introductory plots LUMLOG separately
      CALL FGCOR1
C==========================================================
C==========================================================
C==========================================================
      END



      SUBROUTINE FGCOR1
!     *****************
! Basic plots BLUM2, OLDBIS, LUMLOG separately
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
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!----------------------------------------------------------------------
! Parameters for tables/figures
      DIMENSION    idl(6)
      CHARACTER*16 capt(7)
      CHARACTER*8  fmt(3),fmtx,fmty
      LOGICAL gexist
!---------------------------------------------------------------------- 
      DIMENSION BorW(10),BorN(10)
!---------------------------------------------------------------------- 
! Mark plots for plots
      CHARACTER*32 star,diamond,circle,times,disc,plus,box,dot
      PARAMETER (diamond ='\\makebox(0,0){\\LARGE $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\LARGE $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (times   ='\\makebox(0,0){\\LARGE $\\times$}')
      PARAMETER (disc    ='\\circle*{20}')
      PARAMETER (plus    ='\\makebox(0,0){\\LARGE $+$}')
      PARAMETER (box     ='\\makebox(0,0){\\Large $\\Box$}')
      PARAMETER (dot     ='\\circle*{10}')
!----------------------------------------------------------------------
      CHARACTER*64 labsical(60)
      DATA labsical /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put(600,1150){\\makebox(0,0)[t]{\\huge SICAL}}',
     $'%%%%%%%%',
     $' \\put(400,850){\\begin{picture}( 300,400)',
!*** $'     \\put(  0, 80){\\makebox(0,0){\\LARGE $\\times$}}',
     $'     \\put(  0, 80){\\circle*{10}}',
     $'     \\put(100, 80){\\makebox(0,0)[l]{\\LARGE Wide-Wide}} ',
!
     $'     \\put(  0,  0){\\circle{30}}',
     $'     \\put(100,  0){\\makebox(0,0)[l]{\\LARGE Wide-Narrow}}',
     $' \\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 600,  40){\\makebox(0,0)[b]{\\huge  ',
     $'      $ 1-z_{\\min}^{^{\\rm{SICAL}}} $}}',
     $' \\put(  50,1050){\\makebox(0,0)[l]{\\huge ${',
     $'  \\sigma _{_{\\rm{Miss.}}} \\over\\sigma_{_{\\rm{Born}}} }$}}',
     $' \\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
      CHARACTER*64 labBARE1(60)
      DATA labBARE1 /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put(600,1150){\\makebox(0,0)[t]{\\huge BARE1}}',
     $'%%%%%%%%',
     $' \\put(400,850){\\begin{picture}( 300,400)',
     $'     \\put(  0, 80){\\circle*{10}}',
     $'     \\put(100, 80){\\makebox(0,0)[l]{\\LARGE Wide-Wide}} ',
!
     $'     \\put(  0,  0){\\circle{30}}',
     $'     \\put(100,  0){\\makebox(0,0)[l]{\\LARGE Wide-Narrow}}',
     $' \\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 600,  40){\\makebox(0,0)[b]{\\huge  ',
     $'     $ 1-z_{\\min}^{^{\\rm{BARE1}}} $}}',
     $' \\put(  50,1050){\\makebox(0,0)[l]{\\huge ${',
     $'  \\sigma _{_{\\rm{Miss.}}} \\over\\sigma_{_{\\rm{Born}}} }$}}',
     $' \\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
      CHARACTER*64 labCALO2(60)
      DATA labCALO2 /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put(600,1150){\\makebox(0,0)[t]{\\huge CALO2}}',
     $'%%%%%%%%',
     $' \\put(400,850){\\begin{picture}( 300,400)',
     $'     \\put(  0, 80){\\circle*{10}}',
     $'     \\put(100, 80){\\makebox(0,0)[l]{\\LARGE Wide-Wide}} ',
!
     $'     \\put(  0,  0){\\circle{30}}',
     $'     \\put(100,  0){\\makebox(0,0)[l]{\\LARGE Wide-Narrow}}',
     $' \\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 600,  40){\\makebox(0,0)[b]{\\huge  ',
     $'     $ 1-z_{\\min}^{^{\\rm{CALO2}}} $}}',
     $' \\put(  50,1050){\\makebox(0,0)[l]{\\huge ${',
     $'  \\sigma _{_{\\rm{Miss.}}} \\over\\sigma_{_{\\rm{Born}}} }$}}',
     $' \\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
      CHARACTER*64 labSICAL2(60)
      DATA labSICAL2 /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put(600,1150){\\makebox(0,0)[t]{\\huge SICAL2}}',
     $'%%%%%%%%',
     $' \\put(400,850){\\begin{picture}( 300,400)',
     $'     \\put(  0, 80){\\circle*{10}}',
     $'     \\put(100, 80){\\makebox(0,0)[l]{\\LARGE Wide-Wide}} ',
!
     $'     \\put(  0,  0){\\circle{30}}',
     $'     \\put(100,  0){\\makebox(0,0)[l]{\\LARGE Wide-Narrow}}',
     $' \\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 600,  40){\\makebox(0,0)[b]{\\huge  ',
     $'     $ 1-z_{\\min}^{^{\\rm{SICAL2}}} $}}',
     $' \\put(  50,1050){\\makebox(0,0)[l]{\\huge ${',
     $'  \\sigma _{_{\\rm{Miss.}}} \\over\\sigma_{_{\\rm{Born}}} }$}}',
     $' \\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      CHARACTER*64 labSICALu(60)
      DATA labSICALu /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put(600,1150){\\makebox(0,0)[t]{\\huge SICAL}}',
     $'%%%%%%%%',
     $' \\put(500,800){\\begin{picture}( 300,400)',
     $'   \\put(100,160){\\makebox(0,0)[l]{\\LARGE No Exponentiation}}',
!
     $'     \\put(  0, 80){\\circle*{10}}',
     $'     \\put(100, 80){\\makebox(0,0)[l]{\\LARGE Wide-Wide}} ',
!
     $'     \\put(  0,  0){\\circle{30}}',
     $'     \\put(100,  0){\\makebox(0,0)[l]{\\LARGE Wide-Narrow}}',
     $' \\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 600,  40){\\makebox(0,0)[b]{\\huge  ',
     $'      $ 1-z_{\\min}^{^{\\rm{SICAL}}} $}}',
     $' \\put( -20, 950){\\makebox(0,0)[r]{\\LARGE ${',
     $'  \\sigma _{_{\\rm{Miss.}}} ',
     $'  \\over\\sigma_{_{\\rm{Born}}} }$}}',
     $' \\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
      CHARACTER*64 labBARE1u(60)
      DATA labBARE1u /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put(600,1150){\\makebox(0,0)[t]{\\huge BARE1}}',
     $'%%%%%%%%',
     $' \\put(500,800){\\begin{picture}( 300,400)',
     $'   \\put(100,160){\\makebox(0,0)[l]{\\LARGE No Exponentiation}}',
!
     $'     \\put(  0, 80){\\circle*{10}}',
     $'     \\put(100, 80){\\makebox(0,0)[l]{\\LARGE Wide-Wide}} ',
!
     $'     \\put(  0,  0){\\circle{30}}',
     $'     \\put(100,  0){\\makebox(0,0)[l]{\\LARGE Wide-Narrow}}',
     $' \\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 600,  40){\\makebox(0,0)[b]{\\huge  ',
     $'     $ 1-z_{\\min}^{^{\\rm{BARE1}}} $}}',
     $' \\put( -20, 950){\\makebox(0,0)[r]{\\LARGE ${',
     $'  \\sigma _{_{\\rm{Miss.}}} ',
     $'  \\over\\sigma_{_{\\rm{Born}}} }$}}',
     $' \\end{picture}}',
     $'% end-of-label'/
!----------------------------------------------------------------------
!----------------------------------------------------------------------

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

!=========================================================
C --------- LUMLOG -------
      KeyGen = 2
      JDA = 2000  +10000*KEYGEN
      kda = JDA +10000000
      JDU = 4000  +10000*KEYGEN
      kdu = JDU +10000000
      DO itr=1,3
        CALL CUMHIS(KeyGen,jdu+410+itr,kdu+410+itr) ! SICAL2 BhabhaWG
        CALL CUMHIS(KeyGen,jdu+420+itr,kdu+420+itr) ! CALO2  BhabhaWG
        CALL CUMHIS(KeyGen,jdu+430+itr,kdu+430+itr) ! BARE1  BhabhaWG
        CALL CUMHIS(KeyGen,jdu+440+itr,kdu+440+itr) ! SICAL  PL95
        CALL CUMHIS(KeyGen,jdu+450+itr,kdu+450+itr) ! SICAL  PL95
        CALL CUMHIS(KeyGen,jdu+460+itr,kdu+460+itr) ! BARE1  BhabhaWG
      ENDDO


! Missing third order
!--------------------------------------------------------
!----    initialize GPLOT
!--------------------------------------------------------
      CALL GPLINT(0)
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      YMAX =  .0011
      YMIN = -.0011
      CALL PL3MC('[410] SICAL2, Missing O(alf3),  WW, NN, NW, 1-Z$' 
     $     ,kdu+411,kdu+412,kdu+413,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
      CALL PL3MC('[420]  CALO2, Missing O(alf3),  WW, NN, NW, 1-Z$' 
     $     ,kdu+421,kdu+422,kdu+423,YMIN,YMAX,BorW(4),BorN(4),BorN(4))
      CALL PL3MC('[430]  BARE1, Missing O(alf3),  WW, NN, NW, 1-Z$' 
     $     ,kdu+431,kdu+432,kdu+433,YMIN,YMAX,BorW(3),BorN(3),BorN(3))
      CALL PL3MC('[440]  SICAL, Missing O(alf3),  WW, NN, NW, 1-Z$' 
     $     ,kdu+441,kdu+442,kdu+443,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
      YMAX =  .002
      YMIN = -.002
      CALL PL3MC('[450]  SICAL, Mis. O(alf3) UNEXP,  WW, NN, NW, 1-Z$' 
     $     ,kdu+451,kdu+452,kdu+453,YMIN,YMAX,BorW(2),BorN(2),BorN(2))
      CALL PL3MC('[460]  BARE1, Mis. O(alf3) UNEXP,  WW, NN, NW, 1-Z$' 
     $     ,kdu+461,kdu+462,kdu+463,YMIN,YMAX,BorW(3),BorN(3),BorN(3))
!------------THE END OF HISTO WRITING -------------------------
      CALL GPLEND
!--------------------------------------------------------
!=======================================================================
!=======================================================================
!=======================================================================
      ymin=-0.0011
      ymax= 0.0011
      fmtx='f10.2'
      fmty='f10.4'
!**************** NORMALIZATION *****************
      DO itr=1,3
! [410] SICAL2, Missing O(alf3)
      fact = 1/BorN(4)
      IF(itr .EQ. 1) fact = 1/BorW(4)
      CALL gopera(kdu+410+itr,'+',kdu+410+itr,kdu+410+itr,0d0,fact)
      CALL gmimax(kdu+410+itr,ymin,ymax)
      CALL gidopt(kdu+410+itr,'ERRO')
! [420] CALO2, Missing O(alf3)
      fact = 1/BorN(4)
      IF(itr .EQ. 1) fact = 1/BorW(4)
      CALL gopera(kdu+420+itr,'+',kdu+420+itr,kdu+420+itr,0d0,fact)
      CALL gmimax(kdu+420+itr,ymin,ymax)
      CALL gidopt(kdu+420+itr,'ERRO')
! [430] BARE1, Missing O(alf3)
      fact = 1/BorN(3)
      IF(itr .EQ. 1) fact = 1/BorW(3)
      CALL gopera(kdu+430+itr,'+',kdu+430+itr,kdu+430+itr,0d0,fact)
      CALL gmimax(kdu+430+itr,ymin,ymax)
      CALL gidopt(kdu+430+itr,'ERRO')
! [440] SICAL, Missing O(alf3)
      fact = 1/BorN(2)
      IF(itr .EQ. 1) fact = 1/BorW(2)
      CALL gopera(kdu+440+itr,'+',kdu+440+itr,kdu+440+itr,0d0,fact)
      CALL gmimax(kdu+440+itr,ymin,ymax)
      CALL gidopt(kdu+440+itr,'ERRO')
      ENDDO
      ymin=-0.0025
      ymax= 0.0025
      DO itr=1,3
! [450] SICAL, Missing O(alf3) in UNEXP
      fact = 1/BorN(2)
      IF(itr .EQ. 1) fact = 1/BorW(2)
      CALL gopera(kdu+450+itr,'+',kdu+450+itr,kdu+450+itr,0d0,fact)
      CALL gmimax(kdu+450+itr,ymin,ymax)
      CALL gidopt(kdu+450+itr,'ERRO')
! [460] BARE1, Missing O(alf3)
      fact = 1/BorN(3)
      IF(itr .EQ. 1) fact = 1/BorW(3)
      CALL gopera(kdu+460+itr,'+',kdu+460+itr,kdu+460+itr,0d0,fact)
      CALL gmimax(kdu+460+itr,ymin,ymax)
      CALL gidopt(kdu+460+itr,'ERRO')
      ENDDO
!*************************************************
!=========================================================
!  Separate PLOTS for translation into eps files
!=========================================================
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './MisOgam3-SICAL.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      CALL gplot2(kdu+443,' ','*',dot     ,fmtx,fmty) ! N-W
      CALL gplot2(kdu+441,'S','*',circle  ,fmtx,fmty) ! W-W
      CALL gplabel(labSICAL)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './MisOgam3-BARE1.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      CALL gplot2(kdu+433,' ','*',dot     ,fmtx,fmty) ! N-W
      CALL gplot2(kdu+431,'S','*',circle  ,fmtx,fmty) ! W-W
      CALL gplabel(labBARE1)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './MisOgam3-CALO2.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      CALL gplot2(kdu+423,' ','*',dot     ,fmtx,fmty) ! N-W
      CALL gplot2(kdu+421,'S','*',circle  ,fmtx,fmty) ! W-W
      CALL gplabel(labCALO2)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './MisOgam3-SICAL2.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      CALL gplot2(kdu+413,' ','*',dot     ,fmtx,fmty) ! N-W
      CALL gplot2(kdu+411,'S','*',circle  ,fmtx,fmty) ! W-W
      CALL gplabel(labSICAL2)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './MisOgam3u-SICAL.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      CALL gplot2(kdu+453,' ','*',dot     ,fmtx,fmty) ! N-W
      CALL gplot2(kdu+451,'S','*',circle  ,fmtx,fmty) ! W-W
      CALL gplabel(labSICALu)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
! initialize GPLOT
      TeXfile   = './MisOgam3u-BARE1.txp'
      CALL GPLINT(2) 
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------
      CALL gplot2(kdu+463,' ','*',dot     ,fmtx,fmty) ! N-W
      CALL gplot2(kdu+461,'S','*',circle  ,fmtx,fmty) ! W-W
      CALL gplabel(labBARE1u)
!--------------------------------------------------------
! The end of histos/plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
      END


