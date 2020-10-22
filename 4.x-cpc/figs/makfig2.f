      PROGRAM MAIN
!     ***********************************
! To execute:  make makfig2-dvi
!              make makfig2-ps
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
      Tesnam    = 'makfig2'
      OPEN( NOUT, file='output-'//Tesnam)
      CALL GOUTPU(NOUT)
!--------------------------------------------------------
! Stored histograms and corresponding histograms 
!--------------------------------------------------------
      lendan = 0
      lendan = lendan+1
!>      Dname(lendan)  = '../prod2/prod2.data.2037M'  ! April 96
!>      Hname(lendan)  = '../prod2/prod2.hst.2037M'   ! April 96
      Dname(lendan)  = '../prod2/prod2.data'        ! Current
      Hname(lendan)  = '../prod2/bhl.hst'           ! Current
      lendan = lendan+1
!>      Dname(lendan)  = '../obis2/obis2.data.2118M'  ! April 96
!>      Hname(lendan)  = '../obis2/obis2.hst.2118M'   ! April 96
      Dname(lendan)  = '../obis2/obis2.data'        ! Current
      Hname(lendan)  = '../obis2/bhl.hst'           ! Current
      lendan = lendan+1
!>      Dname(lendan)  = '../llog2/llog2.data.2114M'  ! April 96
!>      Hname(lendan)  = '../llog2/llog2.hst.2114M'   ! April 96
      Dname(lendan)  = '../llog2/llog2.data'        ! Current
      Hname(lendan)  = '../llog2/bhl.hst'           ! Current
!==========================================================
! Table for workshop
      CALL wshop
!==========================================================
!--------------------------------------------------------
!  dumping histogram for control
!--------------------------------------------------------
      NOUTH=20
      OPEN(NOUTH,file='dump.hst')
      CALL grfile(nouth,' ','N')
      CALL grout( 0,ICY,' ')
      CALL grend(tname)
      CLOSE(nout)
      END




      SUBROUTINE wshop
!     ************************
! tables for wshop proceedings
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
!-----------------------------------------------------------------------
! Communicates with MAIN
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!-----------------------------------------------------------------------
! Communicates with READAT
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI
!-----------------------------------------------------------------------
      DIMENSION cut(100)
      DIMENSION bin1(100),err1(100)
      DIMENSION bin2(100),err2(100)
      DIMENSION bin3(100),err3(100)
      CHARACTER*80 title
      DIMENSION BorW(10),BorN(10)
!-------------------------
! Parameters for tables
      DIMENSION    idl(5)
      CHARACTER*16 capt(6)
      CHARACTER*8  fmt(3), fmtx,fmty
!-----------------------------------
      CHARACTER*64 cpsical1(50)
      DATA cpsical1 /
     $'\\setcounter{figure}{3}',
     $'For the ALEPH SICAL detector we plot',
     $'$c_{4-2} = ({\\rm BHLUMI}.4 - {\\rm BHLUMI}.2)/{\\rm Born}$',
     $'as a function of the energy cut $1-U_{\\min}$, where',
     $'$U_{\\min}=\\min(E^{\\rm{cl}}_1,E^{\\rm{cl}}_2)/E_{\\rm beam}$.',
     $'The standard value   $(1-U_{min})^{CUT}=0.561341$ ',
     $'is marked with the vertical line.',
     $'Three curves plotted with small dots, open circles and ',
     $'big dots represent angular cuts W-W, N-N and N-W, ',
     $'respectively, where W and N denote',
     $'wide or narrow angular ranges on one side of the detector.',
     $'The wide (W) angular range is',
     $'$\\theta_A+\\Delta < \\theta_1^{\\rm{cl}} < \\theta_B-\\Delta$,',
     $'and the narrow (N) angular range is',
     $'$\\theta_A+2\\Delta<\\theta_2^{\\rm{cl}} < \\theta_B-4\\Delta$,',
     $'where $\\theta_A= 0.024$, $\\theta_B= 0.058$',
     $'and $\\Delta=(\\theta_B-\\theta_A)/16$.',
     $'The other cuts are:',
     $'(a) auxiliary energy cuts',
     $'$Y_{\\min} = 0.60315$, $Z_{\\min} = 0,$',
     $'(b) acoplanarity cut  $\\Delta\\phi_{\\max} = 0.52359$,',
     $'see Fig. 2 for cut-off definitions.',
     $'% end-of-caption'/
!-----------------------------------
      CHARACTER*64 cpsical2(50)
      DATA cpsical2 /
     $'The difference',
     $'$d_{3}=({\\rm BHLUMI}.4x-{\\rm OLDBIS-LUMLOG})/{\\rm Born}$',
     $'for the ALEPH SICAL detector',
     $'as a function of the energy cut $1-U_{\min}$.',
     $'Cuts are the same as in Fig. 3.',
     $'% end-of-caption'/
!-----------------------------------
      CHARACTER*64 cpsical3(50)
      DATA cpsical3 /
     $'Difference of the cross section calculated',
     $'with the ${\\cal O}(\\alpha^2)_{\\rm prag}$ matrix element',
     $'with YFS exponentiation and without  exponentiation',
     $'for the SICAL trigger. The difference ',
     $'divided by the Born cross section is plotted',
     $'as a function of the energy cut $1-U_{\\min}$.',
     $'The other cuts are the same as in Fig. 3.',
     $'% end-of-caption'/
!-------------------------
! Mark plots for plots
      CHARACTER*32 star,diamond,circle,times,disc,dot
      PARAMETER (diamond ='\\makebox(0,0){\\LARGE $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\LARGE $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (times   ='\\makebox(0,0){\\LARGE $\\times$}')
      PARAMETER (disc    ='\\circle*{20}')
      PARAMETER (dot     ='\\circle*{10}')
!-----------------------------------
      CHARACTER*64 labex1(40)
      DATA labex1 /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put( 600, 1170){\\makebox(0,0)[t]{\\Large ',
     $'  ${\\rm{BHLUMI 4.03}-\\rm{BHLUMI 2.01}\\over\\rm{Born} }$ }}',
     $'%%%%%%%%',
     $' \\put(200,120){\\begin{picture}( 300,400)',
!
     $' \\put(  0,120){\\circle*{10}}',
     $' \\put( 50,120){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^3)_{NNL}$}}',
     $' \\put(280,120){\\makebox(0,0)[l]{\\large WW}}',
!
     $' \\put(  0, 60){\\circle*{20}}',
     $' \\put( 50, 60){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^3)_{NNL}$}}',
     $' \\put(280, 60){\\makebox(0,0)[l]{\\large NW}}',
!
     $' \\put(  0,  0){\\circle{30}}',
     $' \\put( 50,  0){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280,  0){\\makebox(0,0)[l]{\\large NN}}',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 900,  40){\\makebox(0,0)[b]{\\LARGE $ 1-U_{\\min} $}}',
     $' \\put( 673,0){\\line(0,1){1200}}',
     $' \\multiput(0,900)(10,0){120}{\\circle*{  2}}',
     $' \\multiput(0,300)(10,0){120}{\\circle*{  2}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!-----------------------------------
      CHARACTER*64 labex2(40)
      DATA labex2 /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put( 600, 1170){\\makebox(0,0)[t]{\\Large ',
     $'  ${\\rm{BHLUMI 4.03} ',
     $'   -\\rm{(OLDBIS+LUMLOG)}\\over\\rm{Born} }$ }}',
     $'%%%%%%%%',
     $' \\put(200,120){\\begin{picture}( 300,400)',
!
     $' \\put(  0,120){\\circle*{10}}',
     $' \\put( 50,120){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^3)_{NNL}$}}',
     $' \\put(280,120){\\makebox(0,0)[l]{\\large WW}}',
!
     $' \\put(  0, 60){\\circle*{20}}',
     $' \\put( 50, 60){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^3)_{NNL}$}}',
     $' \\put(280, 60){\\makebox(0,0)[l]{\\large NW}}',
!
     $' \\put(  0,  0){\\circle{30}}',
     $' \\put( 50,  0){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280,  0){\\makebox(0,0)[l]{\\large NN}}',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 900,  40){\\makebox(0,0)[b]{\\LARGE $ 1-U_{\\min} $}}',
     $' \\put( 673,0){\\line(0,1){1200}}',
     $' \\multiput(0,900)(10,0){120}{\\circle*{  2}}',
     $' \\multiput(0,300)(10,0){120}{\\circle*{  2}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!-----------------------------------
      CHARACTER*64 labex3(50)
      DATA labex3 /
     $' \\put(300,250){\\begin{picture}( 1200,1200)',
     $' \\put( 600, 1170){\\makebox(0,0)[t]{\\Large ',
     $'  ${ {\\cal O}(\\alpha^2)^{\\rm exp}_{\\rm prag} ',
     $'   - {\\cal O}(\\alpha^2){\\rm prag}  \\over\\rm{Born} }$ }}',
     $'%%%%%%%%',
     $' \\put(200,120){\\begin{picture}( 300,400)',
!
     $' \\put(  0,120){\\circle*{10}}',
     $' \\put( 50,120){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^3)_{NNL}$}}',
     $' \\put(280,120){\\makebox(0,0)[l]{\\large WW}}',
!
     $' \\put(  0, 60){\\circle*{20}}',
     $' \\put( 50, 60){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^3)_{NNL}$}}',
     $' \\put(280, 60){\\makebox(0,0)[l]{\\large NW}}',
!
     $' \\put(  0,  0){\\circle{30}}',
     $' \\put( 50,  0){\\makebox(0,0)[l]{\\large ',
     $'            ${\\cal O}(\\alpha^2)_{exp}$}}',
     $' \\put(280,  0){\\makebox(0,0)[l]{\\large NN}}',
     $'\\end{picture}}',
     $'%%%%%%%%',
     $' \\put( 900,  40){\\makebox(0,0)[b]{\\LARGE $ 1-U_{\\min} $}}',
     $' \\multiput(0,900)(10,0){120}{\\circle*{  2}}',
     $' \\multiput(0,300)(10,0){120}{\\circle*{  2}}',
     $'\\end{picture}}',
     $'% end-of-label'/
!==================================================================

! ------------------------------------------------
! ----------------- BHLUMI -----------------------
! ------------------------------------------------
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
      JBHL = 30000 +2000
      kbhl = JBHL  +10000000
      LBHL = 30000 +4000
      mbhl = LBHL  +10000000
!------- normalize in nanobarns
      KeyGen = 3
      DO itr=1,3
! BHLUMI O(alf2)exp
        CALL cumhis(KeyGen,JBHL+300+itr,kbhl+300+itr) ! SICAL
! BHLUMI  4.03 - 2.02  difference [Bhlumi4-Bhlumi2 'Correction']
        CALL cumhis(KeyGen,LBHL+300+itr,mbhl+300+itr) ! SICAL
! (B) Switching-off exponentiation O(alf2)eB-O(alf2)B (Phys. Lett. 95)
        CALL cumhis(KeyGen,LBHL+310+itr,mbhl+310+itr) ! SICAL

        CALL cumhis(KeyGen,JBHL+340+itr,kbhl+340+itr) ! BARE1
        CALL cumhis(KeyGen,JBHL+350+itr,kbhl+350+itr) ! CALO2
        CALL cumhis(KeyGen,JBHL+360+itr,kbhl+360+itr) ! SICAL2
      ENDDO
! ------------------------------------------------
! ----------------- OLDBIS -----------------------
! ------------------------------------------------
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
      JBIS= 10000 +2000
      kbis= JBis +10000000
!------- normalize in nanobarns
      KeyGen = 1
      DO itr=1,3
        CALL cumhis(KeyGen,JBIS+300+itr,kbis+300+itr) ! SICAL
        CALL cumhis(KeyGen,JBIS+340+itr,kbis+340+itr) ! BARE1
        CALL cumhis(KeyGen,JBIS+350+itr,kbis+350+itr) ! CALO2
        CALL cumhis(KeyGen,JBIS+360+itr,kbis+360+itr) ! SICAL2
      ENDDO
! ------------------------------------------------
! ----------------- LUMLOG -----------------------
! ------------------------------------------------
      jene=3
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
      JLOG= 20000 +4000
      klog= JLOG +10000000
!------- normalize in nanobarns
      KeyGen = 2
      DO itr=1,3
! O(alf3-alf1) !!!
        CALL cumhis(KeyGen,JLOG+300+itr,klog+300+itr) ! SICAL
        CALL cumhis(KeyGen,JLOG+340+itr,klog+340+itr) ! BARE1
        CALL cumhis(KeyGen,JLOG+350+itr,klog+350+itr) ! CALO2
        CALL cumhis(KeyGen,JLOG+360+itr,klog+360+itr) ! SICAL2
      ENDDO
!==================================================================
!==================================================================
!             B-(O+L)
!==================================================================
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
!-------------------
      jhyb= 50000 +2000
      kdif = jhyb +10000000
      DO itr=1,3
      fact = 1/BorN(2)
      IF(itr .EQ. 1) fact = 1/BorW(2)
! (B-(O+L))/Born
      CALL gopera(kbis+300+itr,'+',klog+300+itr,jhyb+300+itr,1d0,1d0)
      CALL gopera(kbhl+300+itr,'-',jhyb+300+itr,kdif+300+itr,fact,fact)
! BHLUMI  4.03 - 2.02  difference [Bhlumi4-Bhlumi2 'Correction']
      CALL gopera(mbhl+300+itr,'+',mbhl+300+itr,mbhl+300+itr,fact,0d0)
! (B) Switching-off exponentiation O(alf2)eB-O(alf2)B (Phys. Lett. 95)
      CALL gopera(mbhl+310+itr,'+',mbhl+310+itr,mbhl+310+itr,fact,0d0)
      ENDDO
!=========================================================
!  PLOTS PLOTS PLOTS PLOTS PLOTS
!=========================================================
!--------------------------------------------------------
      ymin=-3d-3
      ymax= 3d-3
      fmtx='f10.2'
      fmty='f10.3'
!--------------------------------------------------------
!----    initialize GPLOT
!--------------------------------------------------------
ccc      CALL GPLINT( 2)
      CALL GPLINT( 0)
      NOUFIG=11
      TeXfile   = 'fig-95-38.tex'
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!=========================================================
!=========================================================
!-----------------------------------------------!
!   Corr 4x-2x      SICAL WW, NN, NW            !
!-----------------------------------------------!
      CALL gpltit(' PL95 Correction $')
      DO itr=1,3
      CALL gidopt(mbhl+300+itr,'ERRO')
      CALL gmimax(mbhl+300+itr,ymin,ymax)
      ENDDO
!-------------------------------------------------------
      CALL gplcapt(cpsical1)
      CALL gplot2(mbhl+301,' ','*',dot    ,fmtx,fmty) !WW
      CALL gplot2(mbhl+302,'S','*',circle ,fmtx,fmty) !NN
      CALL gplot2(mbhl+303,'S','*',disc   ,fmtx,fmty) !NW
      CALL gplabel(labex1)
!=========================================================
!-----------------------------------------------!
!   (B-(O+L))/Born  SICAL WW, NN, NW            !
!-----------------------------------------------!
      CALL gpltit(' PL95 B-(O+L) $')
      DO itr=1,3
      CALL gidopt(kdif+300+itr,'ERRO')
      CALL gmimax(kdif+300+itr,ymin,ymax)
      ENDDO
!-------------------------------------------------------
      CALL gplcapt(cpsical2)
      CALL gplot2(kdif+301,' ','*',dot    ,fmtx,fmty)
      CALL gplot2(kdif+302,'S','*',circle ,fmtx,fmty)
      CALL gplot2(kdif+303,'S','*',disc   ,fmtx,fmty)
      CALL gplabel(labex2)
!=========================================================
!-----------------------------------------------!
!   (B-(O+L))/Born  SICAL WW, NN, NW            !
!-----------------------------------------------!
      CALL gpltit(' PL95 B-(O+L) $')
      DO itr=1,3
      CALL gidopt(mbhl+310+itr,'ERRO')
      CALL gmimax(mbhl+310+itr,ymin,ymax)
      ENDDO
!-------------------------------------------------------
      CALL gplcapt(cpsical3)
      CALL gplot2(mbhl+311,' ','*',dot    ,fmtx,fmty)
      CALL gplot2(mbhl+312,'S','*',circle ,fmtx,fmty)
      CALL gplot2(mbhl+313,'S','*',disc   ,fmtx,fmty)
      CALL gplabel(labex3)
!--------------------------------------------------------
!   The end of plots writing
      CALL gplend
!=========================================================
      WRITE(   6,BXCLO)
      END

      SUBROUTINE subtra(id1,id2,id3)
!     ******************************
! subrats id2 from id1 and divides by id2
! Errors are only divided, not combined as in gopera!!!
!     *********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION bin1(100),bin2(100)
      DIMENSION err1(100),err2(100)

      CALL gunpak(id1,bin1,' ',0)
      CALL gunpak(id1,err1,'ERRO',0)
      CALL gunpak(id2,bin2,' ',0)
      DO i=1,100
         bin1(i)=(bin1(i)-bin2(i))/bin2(i)
         err1(i)=          err1(i)/bin2(i)
      ENDDO
      CALL gpak (id3,bin1)
      CALL gpake(id3,err1)
      END
