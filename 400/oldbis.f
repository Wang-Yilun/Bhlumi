      SUBROUTINE OLDBIS(MODE,XPAR,NPAR)
C     *********************************
C
C           **************************************************    
C           *       **********************************       *
C           *       *      *******************       *       *
C           *       *      *                 *       *       *
C           *       *      *   O L D B I S   *       *       *
C           *       *      *                 *       *       *
C           *       *      *******************       *       *
C           *       **********************************       *
C           **************************************************    
C
C======================================================================
C
C OLDBIS is an improved version of OLDBAB Monte Carlo of
C F. Berends and R. Kleiss Nucl. Phys. B228 (1983) 537.    
C
C OLDBIS was extensively used and partly described in 
C Phys. Lett. B253 (1991) 469, TH-5888.
C
C OLDBIS represents the true QED first order calculation 
C for small angle Bhabha ( < 10degr.) to within  0.02% TECHNICAL 
C precision, see Phys. Lett. B253 (1991) 469  for more details.
C
C For those who are too young to rememeber: 
C OLDBAB was the standard and unique Monte Carlo for luminosity 
C calculations in all PETRA and PEP experiments!
C
C In the following we anwer or commentent on the following questions:
C   (1) What are essential advantages of OLDBIS with respect to OLDBAB
C   (2) What corrections were done with respect to the original source
C   (3) How to use the program?
C======================================================================
C
C Answer to question (1) :
C ======================== 
C
C The most important advantage of OLDBIS with respect to OLDBAB is that
C OLDBIS has well established TECHMICAL precision of 0.02% in a very
C wide range of cut-off parameters, 
C see Phys. Lett. B253(1991)469,TH-5888.
C
C OLDBAB calculates first order QED correction to small angle Bhabha.
C The basic and burning question in the beginning of 1990 was: 
C to what TECHNICAL precision this
C first order correction from OLDBAB (and other programs like BABAMC)
C represent the true and unique QED answer!
C One could only guess that it is about 1%  but no firm statement with 
C solid justification could be found at the time
C in any published papers on the above guess!
C (N.B. 1% precision was enough for PETRA/PEP.)
C In Phys. Lett. B253(1991)469, we have proved that, after some 
C modifications (transforming OLDBAB into OLDBIS), the fantastic
C TECHNICAL precision 0.02% was reached! To this precision OLDBIS
C represents the TRUE first order QED calculation in the 
C low angle Bhabha.
C 
C Answer to question (2) :
C ========================        
C 
C We have started with source code of OLDBAB taken from RADCOR package 
C installed on CERN public disk by R. Kleiss.
C The algorithm and most of source code in present OLDBIS is identical 
C to that of OLDBAB.
C New corrections are of two types: 
C (a) cosmetic ones which concern mainly output/input formats, choice 
C of random number generator etc.
C (b) essential ones which modify
C algorithm, implementation of QED matrix element and usage of program.
C 
C Full list of corrections/modifications:
C 
C (a) Midifications of input/output and other of the cosmetic type:
C -----------------------------------------------------------------
C
C (==( 000 )==) 
C The histograming routines are removed.
C
C (==( 001 )==) 
C The input/output structure was aligned with other programs like
C BHLUMI, KORALZ, KORALB, LESKOF etc.
C In the CALL OLDBIS(MODE,XPAR,NPAR) initialisation, generation
C and postgeneration phases are determined by MODE= -1,0,1 parameter.
C
C (==( 002 )==) 
C All input is read through XPAR and NPAR in the obligatory
C initialization for MODE=-1 and immetiately printed out.
C
C (==( 003 )==) 
C Each generated event (for MODE=0) is encoded in /MOMBAB/.
C
C (==( 004 )==) 
C Calling OLDBIS with MODE=1 or 2 provides 
C cross-section resulting from Monte Carlo integration (necessary
C for histogram normalization).  
C For MODE=1 no output is printed.
C
C (==( 005 )==)
C For MODE=2 in addition to information encoded in XPAR and NPAR 
C certain output is printed. 
C Note that routine BABINF printing output in the same format 
C as in original OLDBAB is kept active.
C
C (==( 006 )==)
C All weight accounting is not done 'by hand' but rather using 
C special routine GMONIT from KORALZ.
C
C (==( 007 )==)
C The modern random number generator VARRAN replaces the original one.
C
C (==( 008 )==)
C Output four-momenta QP,QM,QK are now in GeV units.
C
C (==( 009 )==)
C The other side photon may go outside angular range, user has
C to reject it by himself if he wants.
C
C (b) Modifications of the algorithm and QED matrix element: 
C ----------------------------------------------------------
C
C (==( 100 )==)  
C Detailed insight into MC algorithm of OLDBAB reveals that variable 
C SIGS is a dummy variable i.e. none of the calculated x-sections 
C and distributions depends (within stat. err.) on SIGS. 
C For this to be true SIGS has to be positive, however!
C For very small k0 SIGS becomes negative and the results from 
C the OLDBAB do not represent true first order QED anymore.
C We have corrected this i.e. defined SIGS which is always positive. 
C
C (==( 101 )==)  
C The above corr. (100) is still not enough to avoid problems at the 
C very small k0=XK0 limit which has to be taken in order to check that 
C the program represent first order QED at the 0.02% technical 
C precision level. One has to alow negative weights i.e. introduction
C of weighted events is necessary. This option is implemented
C and the program provides weighted events for switch KEYWGT=1.
C The weight WTM is located in /MOMBAB/ and later transfered 
C to /WGTALL/.  We recommend to use OLDBIS with weighted events.
C
C (==( 102 )==)  
C The essential observation made in Phys. Lett. B253 (1991) 469,TH-5888
C and expoited also in other papers (TH-5995,TH-6118) is that certain
C class of anoying QED corrections, so called up-down interferences,
C is completely unimportant at low angles ( <10 degr.).
C To obtain this result we had to rewrite completely the QED soft and
C hard matrix element in OLDBAB. 
C New routine VIRSEL calculating soft and virtual corrections is added.
C The user has at his disposal four types of matrix element with
C various components switched on/off. 
C Each of them is represented by the separate model weight XWT(10:13). 
C Only one of them 
C                      XWT(10+KEYSIN) 
C is choosen as a model weight 
C for eventual rejection, where KEYSIN=0,1,2,3 is the input parameter:
C               Table 1, WTSET entries for OLDBIS
C----------------------------------------------------------------------
C  Entry        Corrections present/absent in matrix element
C----------------------------------------------------------------------
C               up-down int.   vac.pol.   Z-exch.  s-chan.phot.  
C               -------------------------------------------------------
C  XWT(10)        yes          yes         yes        yes     
C  XWT(11)         no           no          no         no       
C  XWT(12)        yes           no          no         no   principal)
C  XWT(13)        yes           no          no        yes   
C  -------------------------------------------------------------------
C The backward compatibility is kept, XWT(10) represents original
C OLDBAB matrix element (with new vacuum polarization and Z-width).
C The difference XWT(12)-XWT(11) accounts for pure up-down interferece.
C
C (==( 103 )==)  
C The archaic vacuum polarization subprogram REPI (for KEYSIN=0) 
C is replaced by the modern version of the same name.
C
C (==( 104 )==)
C The Z-gamma inerference correction includes now Z width 
C (for KEYSIN=0).
C
C (==( 105 )==) 
C No immediate rejection for events in which one fermion is in 
C theta-trigger but another one (due to photon emission) is out.
C Such an event goes through but has zero weight.
C
C (==( 106 )==) 
C Symmetrisation QP <=> QM can be supressed by setting KEYMIR switch
C KEYMIR = 0. This option was usefull in TH-5888 comparisons
C with semianalytical calculations.
C In the following common block TRAN is transfer to fermion without
C photon emission BEFORE symmetrization
C     COMMON / TRANSR / TRAN,TRANP,TRANQ
C
C Answer to question (3) :
C ======================== 
C
C The complete description of the usage of the program can be found
C in Long-write-up of BHLUMI 2.01, CERN preprint TH-6230.
C Here we only summarize on Input/Output parameters in
C                 CALL BHLUMI(MODE,XPAR,NPAR)
C The user may use directly OLDBIS as well
C                 CALL OLDBIS(MODE,XPAR,NPAR)
C
C IF( MODE =-1 ) THEN
C ===================
C   
C Initialisation is performed, all input parameters are transferred
C through XPAR and NPAR, see table below:
C     Table 2, Input parameters of OLDBIS
C----------------------------------------------------------------------
C  Entry    Variable   Meaning
C----------------------------------------------------------------------
C  NPAR( 1)  KEYOPT =1000*KEYGEN +10*KEYWGT +KEYRND 
C                   General option switch 
C            KEYGEN =1 for this sub-generator
C            KEYRND =1,2 type of random number generator RANMAR,RANECU
C            KEYWGT =0,1 for constant, variable weight WTM
C  NPAR( 2)  KEYRAD =10*KEYMIR +KEYSIN, QED switch defining WTM weight
C                   The meaning of KEYSIN is summarized in Table 1,
C            KEYMIR =0 photon emitted only from both fermions 
C            KEYMIR =1 photon emitted only from QM fermion line (tests)
C  XPAR( 1)  CMSENE Total center mass energy [GeV]
C  XPAR( 2)   THMIN Minimum theta angle for electron [degr]
C  XPAR( 3)   THMAX Maximum theta angle for electron [degr]
C  XPAR( 4)     XK0 k0 parameter, infrared cut on real photon energy
C                   is k0*CMSENE/2, recomended range 0.000001<k0<0.001
C  XPAR( 5)   XKMAX maximum real photon energy is XKMAX*CMSENE/2
C  XPAR( 6)   XKMIN minimum real photon energy is XKMIN*CMSENE/2,
C                   one should normaly set XKMAX=1, XKMIN=0.
C----------------------------------------------------------------------
C
C ELSE IF( MODE = 0 ) THEN
C ========================
C         
C Generation of the single Monte Carlo event
C The four momenta of the final state electron positron and photon
C are primarily encoded in 
C      COMMON / MOMBAB / QP(4),QM(4),QK(4),WTM,XWT(20)
C where QP and QM are four-momenta of positron and electron and QK
C is four-momentum of the photon (from time to time QK is zero).
C They are in GeV units and the z-axis points in the direction of
C the positron beam.
C WTM is the main model weight selected according to KEYSIN switch.
C All four possible model weights with various contributions on/off
C (see table 1) are encoded in XWT. For KEYWGT=0 we have WTM=1
C and XWT should not be used!
C Note that content of /MOMBAB/ is copied to standard multiphoton
C common blocks /MOMSET/ and /WGTALL/ of BHLUMI as well.
C
C
C ELSE IF( MODE = 1 ) THEN
C ========================
C       
C The total cross section corresponding to generated series of event,
C i.e. resulting from MC integrartion is calculated and stored in XPAR,
C together with a lot of auxiliary information, see table 3 below.
C This impressive set of output was nesessary for TH-5888 calculations.
C                         Table 3
C----------------------------------------------------------------------
C  Entry    Variable   Meaning
C----------------------------------------------------------------------
C  NPAR(10)  NEVGEN  Number of generated MC events
C  NPAR(20)  NEVGEN  Number of generated MC events
C  XPAR( 9)    SIG0  Born x-section  [nb]
C  XPAR(10)    XST2  Total x-section [nb]
C  XPAR(11)   RXST2  The relative error of XST2
C  XPAR(12)     XSS  Soft (k<k0)  x-section [nb]
C  XPAR(13)    RXSS  The relative error of XSS
C  XPAR(14)     XSH  Hard (k<k0)  x-section [nb]
C  XPAR(15)    RXSH  The relative error of XSH
C  XPAR(20)    SIGS+SIGH  Crude total MC x-section [nb] which 
C                    is necessary for rescaling histograms in 
C                    run with weighted events.
C  XPAR(21)          =0, error of XPAR(20) is zero
C  XPAR(22)    SIGS  Crude soft MC x-section [nb]
C  XPAR(23)          =0, error of XPAR(22) is zero
C  XPAR(24)    SIGH  Crude hard MC x-section [nb]
C  XPAR(25)          =0, error of XPAR(24) is zero
C----------------------------------------------------------------------
C For constant weight option KEYWGT=0 (convevience in rescaling histos)
C we put XPAR(20,21,22)=XPAR(10,11,12) !
C For MODE=1 program is called upon many times in the process of 
C rescaling histograms there is no output printed.
C
C ELSE IF( MODE = 2 ) THEN
C ========================                     
C 
C Only in this MODE=2 in addition to filling XPAR and NPAR as 
C for MODE=1 the values of various x-sections are printed on 
C standard output file.
C                         
C ENDIF
C ====
C
C
C
C======================================================================
C======================================================================
C History of corrections for the record (stj) 
C you may skip reading this.
C======================================================================
C===== Series of corrections from Jan. 91 to June 91 ==================
C======================================================================
C===== VERSION OF OLDBAB USED IN TH.5888/90 BENCHMARK =================
C===== DILOG and GMONIT moved to library ==============================
C===== REPI replaced with modern version ==============================
C===== Adjusted KEYRAD and KEYPOT/KEYWGT ==============================
C===== Outout XPAR(11) is relative error ==============================
C===== KINBIS translates MOMBAB into MOMSET ===========================
C======================================================================
C===== Series of corrections from Jan. 90 to Dec. 90 ==================
C======================================================================
C  (0) IMODE=-1,0,1,2   initialization, production, output
C  (1) new random number generator
C  (2) QP,QM,PK moved to /MOMBAB/
C  (3) QP,QM,PK in GeV units
C  (4) negative weights in soft part treated properly
C  (5) possibility of weighted events, including negative
C  (6) reorganized QED matrix element with up-down interf. isolated 
C======================================================================
C======================================================================
C
C        Here starts 'original' source code
C
C------------------------------------------- REMARKS ---------------
C
C
C SIMULATION OF RADIATIVE BHABHA SCATTERING IN 2ND & 3RD ORDER Q.E.D.
C WRITTEN BY R.KLEISS (LEIDEN & DESY) OCTOBER '82
C
C
C********************************************************************
C
C A DETAILED WRITEUP, TOGETHER WITH PHYSICS INTERPRETATION
C CAN BE FOUND IN  :
C         F.A. BERENDS AND R. KLEISS , NUCL. PHYS. B228(1983)537 .
C IF YOU USE THIS PROGRAM, PLEASE MAKE REFERENCE TO THE ABOVE PAPER !
C
C********************************************************************
C
C
C BEFORE CALL OF FIRST EVENT, COMMON 'PARM01' MUST BE FILLED WITH :
C EBEAM = BEAM ENERGY (IN GEV) ;
C THMIN = MINIMUM SCATTERING ANGLE OF ELECTRON/POSITRON (IN DEGREES)
C THMAX = MAXIMUM SCATTERING ANGLE OF ELECTRON/POSITRON (IN DEGREES)
C XKMIN = MINIMUM BREMSSTRAHLUNG ENERGY (IN UNITS OF EBEAM) ;
C XKMAX = MAXIMUM BREMSSTRAHLUNG ENERGY (IN UNITS OF EBEAM) .
C XKMIN=0,XKMAX=1 ARE ALLOWED. THMIN=0,THMAX=180 ARE NOT ALLOWED.
C IF XKMIN < 0.01 IT IS ASSUMED TO BE 0.
C QP(I) = FOUR MOMENTUM OF OUTGOING POSITRON (I=1,2,3,4) ;
C QM(I) = FOUR MOMENTUM OF OUTGOING ELECTRON (I=1,2,3,4) ;
C QK(I) = FOUR MOMENTUM OF PHOTON            (I=1,2,3,4) .
C FOUR MOMENTA ARE IN UNITS OF EBEAM. THE INCOMING MOMENTA ARE:
C POSITRON (0,0,1,1) , ELECTRON (0,0,-1,1) .
C--------------------------------------------------------------------
C LATEST UPDATE: APRIL 5,1983.
C--------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      PARAMETER(PI= 3.1415926535897932D0, ALFINV=137.03604D0)   
      PARAMETER(GNANOB =  389385D0) 
CCCCCCCCC                       (==( 002 )==)
      DIMENSION XPAR(*),NPAR(*)
CCCCCCCCC                       (==( 003 )==)
      COMMON / MOMBAB / QP(4),QM(4),QK(4),WTM,XWT(20)
      COMMON / PARM01 / EBEAM,THMIN,THMAX,XKMIN,XKMAX,CMIN,CMAX
      COMMON / PARM02 / SIG0,SIGS,SIGH,SIGT,WT
      COMMON / INOUT  / NINP,NOUT 
      COMMON / CMONIT/ AVERWT,ERRELA,NEVTOT,NEVACC,NEVNEG,NEVOVE,NEVZER
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G
C Communicates with VARRAN      (==( 007 )==)
      COMMON / RANPAR / KEYRND
C True transfer (==( 106 )==) 
      COMMON / TRANSR / TRAN,TRANP,TRANQ
      DOUBLE PRECISION DRVEC(100)
      DIMENSION WT(18)
CCCCCCCCC                       (==( 001 )==)
C     ===================
      IF(MODE.EQ.-1) THEN
C     ===================
C-------------------------------------------- INITIALIZATION --------
C ...BX-formats for nice and flexible outbuts
      BXOPE =  '(//1X,15(5H*****)    )'
      BXTXT =  '(1X,1H*,                  A48,25X,    1H*)'
      BXL1I =  '(1X,1H*,I17,                 16X, A20,A12,A7, 1X,1H*)'
      BXL1F =  '(1X,1H*,F17.8,               16X, A20,A12,A7, 1X,1H*)'
      BXL2F =  '(1X,1H*,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H*)'
      BXL1G =  '(1X,1H*,G17.8,               16X, A20,A12,A7, 1X,1H*)'
      BXL2G =  '(1X,1H*,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H*)'
      BXCLO =  '(1X,15(5H*****)/   )'
      CALL GMONIT(-1,70,0D0,1D0,1D0)
      CALL GMONIT(-1,71,0D0,1D0,1D0)
      CALL GMONIT(-1,72,0D0,1D0,1D0)
      CALL GMONIT(-1,41,0D0,1D0,1D0)
      CALL GMONIT(-1,42,0D0,1D0,1D0)
      CALL GMONIT(-1,51,0D0,1D0,1D0)
      CALL GMONIT(-1,52,0D0,1D0,1D0)
CCCCCCCCC                       (==( 002 )==)
      KEYOPT=NPAR(1)
      KEYRAD=NPAR(2) 
      CMSENE=XPAR(1)
      THMIN =XPAR(2)
      THMAX =XPAR(3)
      XK0   =XPAR(4)
      XKMAX =XPAR(5)
      XKMIN =XPAR(6)
      EBEAM =CMSENE/2
      KEYRND = MOD(KEYOPT,10)   
      KEYWGT = MOD(KEYOPT,100)/10   
      KEYSIN = MOD(KEYRAD,10)
      KEYMIR = MOD(KEYRAD,100)/10   
CCCCCCCCC                       (==( 002 )==)
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '     **********************      '
      WRITE(NOUT,BXTXT) '     **    O L D B I S   **      '
      WRITE(NOUT,BXTXT) '     **********************      '
      WRITE(NOUT,BXTXT) ' O(alpha) Monte Carlo Program    '
      WRITE(NOUT,BXTXT) ' for the small-angle Bhabha scat.'
      WRITE(NOUT,BXTXT) '         --------------          '
      WRITE(NOUT,BXTXT) '         Important Note          '
      WRITE(NOUT,BXTXT) '         --------------          '
      WRITE(NOUT,BXTXT) ' This is an improved version     '
      WRITE(NOUT,BXTXT) ' of the OLDBAB M.C. program of   '
      WRITE(NOUT,BXTXT) ' F.A. BERENDS AND R. KLEISS      '
      WRITE(NOUT,BXTXT) ' [1] NUCL. PHYS. B228 (1983) 537 '
      WRITE(NOUT,BXTXT) ' Changes were done by S. Jadach, '
      WRITE(NOUT,BXTXT) ' E. Richter-Was and other        '
      WRITE(NOUT,BXTXT) ' authors of the paper            '
      WRITE(NOUT,BXTXT) ' [2] Phys. Lett. B253 (1991) 469 '
      WRITE(NOUT,BXTXT) ' All modifications are desribed  '
      WRITE(NOUT,BXTXT) ' in detail in the source code.   '
      WRITE(NOUT,BXTXT) ' PLEASE CITE references [1,2] !  '
      WRITE(NOUT,BXCLO)
C
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '* This program is part of       *'
      WRITE(NOUT,BXTXT) '*   BHLUMI version 2.01         *'
      WRITE(NOUT,BXTXT) '*   September      1991         *'
      WRITE(NOUT,BXTXT) '*         AUTHORS               *'
      WRITE(NOUT,BXTXT) '* S. Jadach, E. Richter-Was     *'
      WRITE(NOUT,BXTXT) '*    B.F.L. Ward, Z. Was        *'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)
C
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '         OLDBIS input            '
      WRITE(NOUT,BXL1I) KEYOPT,     ' option    switch  ','KEYOPT','N1'
      WRITE(NOUT,BXL1I) KEYRND,     ' VARRAN    switch  ','KEYRND','  '
      WRITE(NOUT,BXL1I) KEYWGT,     ' weighting switch  ','KEYWGT','  '
      WRITE(NOUT,BXL1I) KEYRAD,     ' radiation switch  ','KEYRAD','N2'
      WRITE(NOUT,BXL1I) KEYSIN,     ' interf. s-chan.   ','KEYSIN','  '
      WRITE(NOUT,BXL1I) KEYMIR,     ' QM emiss.only/test','KEYMIR','  '
      WRITE(NOUT,BXL1F) CMSENE,     ' CMSENE            ','CMSENE','X1'
      WRITE(NOUT,BXL1F) THMIN ,     ' THMIN             ','THMIN ','X2'
      WRITE(NOUT,BXL1F) THMAX ,     ' THMAX             ','THMAX ','X3'
      WRITE(NOUT,BXL1F) XK0   ,     ' XK0               ','XK0   ','X4'
      WRITE(NOUT,BXL1F) XKMAX ,     ' XKMAX             ','XKMAX ','X5'
      WRITE(NOUT,BXL1F) XKMIN ,     ' XKMIN             ','XKMIN ','X6'
      WRITE(NOUT,BXCLO)
      NEVGEN=0
      WMAX=2.D0
      DO 1 I=1,18
    1 WT(I)=0.
      WT(13)=WMAX
      WT(17)=WMAX 
      AMEL  = 0.511D-3
      XM2=(AMEL/EBEAM)**2
      EM2=1.D0+.5D0*XM2
      XL2=DLOG(2.D0/XM2)
      TWOPI= 2D0*PI
      CMIN=DCOS(THMAX*TWOPI/360.D0)
      CMAX=DCOS(THMIN*TWOPI/360.D0)
      DMIN=1.D0/(1.D0-CMIN)
      DMAX=1.D0/(1.D0-CMAX)
      XLCMIN=XL2+DLOG(1.D0-CMIN)
      XSOFT=1.D0
      IF(XKMIN.LT.XK0) GOTO 2
      XSOFT=0.D0
      XK0=XKMIN
    2 XKRAT=XKMAX/XK0
      ALFA  =  1D0/ALFINV       
      SIGBS =  8D0*PI*ALFA**2/(2D0*EBEAM)**2 *GNANOB
      SIG0  =  SIGBS/16D0  *(
     $ 16*(DMAX-DMIN) + 16*DLOG(DMIN/DMAX) + 9*(CMAX-CMIN)
     $ + (CMAX**2-CMIN**2) + 1.D0/3.D0*(CMAX**3-CMIN**3)         )
      SIGH =  SIGBS * (ALFA/PI)  *4D0 *DLOG(XKRAT)*
     $ ((XL2+1.D0+DLOG(1.D0-CMAX))*DMAX - (1.D0+XLCMIN)*DMIN )  
CC[[[[[[[[[[[[   (==( 100 )==)
C     SIGS=0.D0
C     DO 3 I=1,100
C     R= FLOAT(I-1)/99.D0
C     C=1.-1./(R*DMAX+(1.-R)*DMIN)
C     IF(KEYSIN.EQ.0) THEN
C       CALL VIRSOF(EBEAM,XK0,C,BORN,CORR)
C     ELSE
C       CALL VIRSEL(KEYSIN,EBEAM,XK0,C,BORN,CORR)
C     ENDIF
C&&&3 SIGS=SIGS+BORN*(1.+CORR)*(1.-C)**2
C   3 SIGS=SIGS+BORN*(1.D0     )*(1.D0-C)**2
C the above assures full indep. of xsect. on k0, stj, sept. 89.
C     WRITE(NOUT,*) ' /////CORR IN SIGS EXCLUDED/////'
C     SIGS=SIGS*TWOPI*(DMAX-DMIN)/100.D0*XSOFT
CC]]]]]]]]]]]]
C SIGS is a dummy parameter
C convention: sigma_sof_crude = sigBS/(1-c)**2 (stj)
      SIGS = SIGBS *(DMAX-DMIN)
      SIGT=SIGS+SIGH
      YSOFT=SIGS/SIGT
      ZSOFT=0.D0
      IF(SIGS.NE.0.D0)  ZSOFT=(DMAX-DMIN)/SIGS*TWOPI
      DTOT=SIGT/SIG0-1.
      WRITE(NOUT, 4) XM2   ,EM2   ,XL2   ,TWOPI ,CMIN  ,CMAX
     .       ,DMIN  ,DMAX  ,XLCMIN,XK0   ,XSOFT ,XKRAT
     .       ,WMAX  ,SIG0  ,SIGH  ,SIGS  ,SIGT  ,YSOFT
     .       ,ZSOFT ,DTOT
    4 FORMAT('0',50('=')/,
     . ' INITIALIZATION FOR BHABHA SCATTERING'/,6(' ',4D15.6/))
      WRITE(NOUT,5)EBEAM,THMIN,THMAX,XK0,XKMAX,SIG0,SIGS,SIGH,SIGT,DTOT
    5 FORMAT(
     . '                        BEAM ENERGY =',F15.6,' GEV'/,
     . '           MINIMUM SCATTERING ANGLE =',F15.6,' DEGREES'/,
     . '           MAXIMUM SCATTERING ANGLE =',F15.6,' DEGREES'/,
     . ' MINIMUM HARD BREMSSTRAHLUNG ENERGY =',F15.6/,
     . ' MAXIMUM HARD BREMSSTRAHLUNG ENERGY =',F15.6/,
     . '         LOWEST ORDER CROSS SECTION =',D15.6,' NB'/,
     . ' APPROX. CROSS SECTION IN SOFT PART =',D15.6,' NB'/,
     . ' APPROX. CROSS SECTION IN HARD PART =',D15.6,' NB'/,
     . '        APPROX. CROSS SECTION TOTAL =',D15.6,' NB'/,
     . '           APPROX. TOTAL CORRECTION =',F15.6)
!-------------------------------------------------------------
! This is Generator Identificator
      IDGEN =  1        
! Important histo which remembers total x-section 
      CALL GMONIT(  -1, IDGEN,0D0,2*SIGT,1D0)          
C     ======================
      ELSEIF(MODE.EQ.0) THEN
C     ======================
      NEVGEN=NEVGEN+1
      CALL GMONIT(  0, IDGEN, SIGT, 2*SIGT,1D0)
    6 CONTINUE
      WTM=1D0
      WTK=1D0
C-------------------------------------------- CHOOSE HARD OR SOFT ---
CCCCCCCCC                       (==( 007 )==)
      CALL VARRAN(DRVEC,1)
      IF(DRVEC(1).LT.YSOFT) GOTO 11
C-------------------------------------------- HARD PHOTON PART ------
      WT(1)=WT(1)+1.D0
C-------------------------------------------- GENERATE K VALUE ------
      CALL VARRAN(DRVEC,1)
      XK=XK0*XKRAT**DRVEC(1)
C-------------------------------------------- GENERATE C VALUE ------
    7 CONTINUE
      CALL VARRAN(DRVEC,2)
      C=1.D0-1.D0/(DMIN+DRVEC(1)*(DMAX-DMIN))
      R=DRVEC(2)*XLCMIN/(XL2+DLOG(1.D0-C))
      WT(2)=WT(2)+1.D0
      IF(R.GT.1.D0) GOTO 7
      CM=2.D0*(1.D0-C)
      SC=DSQRT(1.D0-C*C)
C-------------------------------------------- GENERATE FI VALUE -----
      FI=R*TWOPI
C-------------------------------------------- GENERATE U VALUE ------
      D=XM2/CM
      CALL VARRAN(DRVEC,1)
      R=-1.D0+2.D0*DRVEC(1)
      V=(D/(1.D0+D))**DABS(R)
      U=((1.D0+D)*V-D)/(1.D0+V)
      E2=U*(1.D0-U)*CM
      EV=DSQRT(1.D0-E2)
      E2=XM2+E2
      IF(R.LT.0.D0) U=1.D0-U
C-------------------------------------------- GENERATE C1 VALUE -----
      CALL VARRAN(DRVEC,2)
      R=DRVEC(1)
      VC=2.D0*E2*(1.D0-R)/(E2+2.D0*EV*(EM2+EV)*R)
      C1=1.D0-VC
      SC1=DSQRT(VC*(2.D0-VC))
C-------------------------------------------- GENERATE F1 VALUE -----
      F1=TWOPI*DRVEC(2)
      CF1=DCOS(F1)
      SF1=DSIN(F1)
C-------------------------------------------- CONSTRUCT QK DIRECTION
      UC=-1.D0+U-U*C
      QK1=(UC*SC1*CF1-U*SC*C1)/EV
      QK2=SC1*SF1
      QK3=(U*SC*SC1*CF1+UC*C1)/EV
      CG=C*QK3+SC*QK1
C-------------------------------------------- REJECT CT VALUES ------
      XKM=1.D0-XK
      X=2.D0*XKM/(2.D0-XK+XK*CG)
      XT=2.D0-X-XK
      CT=(X*C+XK*QK3)/XT           
CCCCC               (==( 009 )=)
* This correction is important for calorimetric trigger!!!!
****      IF(CT.LT.CMIN.OR.CT.GT.CMAX) WTK=0D0
      WT(3)=WT(3)+ WTK
C-------------------------------------------- CALCULATE WEIGHT ------
      S =4.D0
      S1=4.D0*XKM
      T =-2.D0*X *(1.D0-C )
      T1=-2.D0*XT*(1.D0-CT)
      U =-2.D0*XT*(1.D0+CT)
      U1=-2.D0*X *(1.D0+C )
      X1=XK*(EM2-QK3)
      X2=XK*(EM2+QK3)
      DY=.5D0*XM2*XK/XKM
      Y1=2.D0*(1.D0-XT)+DY
      Y2=2.D0*(1.D0-X )+DY
C Define true transfer (==( 106 )==) 
      tran= ebeam**2*(-T)
CC[[[[[[[[[[[[         (==(102)==)
*.....KEYSIN=0   original OLDBAB 
*.....KEYSIN=1   T NONINTERFERENCE CHANNEL
*.....KEYSIN=2   T CHANNEL WITH INTERF.
*.....KEYSIN=3 S+T CHANNEL
      WTH30=(S*S1*(S*S+S1*S1)+T*T1*(T*T+T1*T1)+U*U1*(U*U+U1*U1))
     .  /(4.*S**3*S1)
     .  *(1.D0-(S*Y1*Y2+S1*X1*X2+U*X2*Y1+U1*X1*Y2)
     .         /(T*X2*Y2+T1*X1*Y1))
     .  *(1.D0-XM2*XK/(1.D0+XKM*XKM)
     .        *(XKM/X1+XKM/X2+1.D0/Y1+1.D0/Y2))
      XWT(10)=  WTK*WTH30
      XWT(13)=  WTK*WTH30
C Here t-channel only and no interferences
      XTCH= (S**2 +U**2 )/(-X1*Y1*T1)*(1+2*XM2*X1/Y1/T1)
     $     +(S1**2+U1**2)/(-X1*Y1*T1)*(1+2*XM2*Y1/X1/T1)
     $     +(S**2 +U1**2)/(-X2*Y2*T )*(1+2*XM2*X2/Y2/T )
     $     +(S1**2+U**2 )/(-X2*Y2*T )*(1+2*XM2*Y2/X2/T )
      XES2= 4*S**2*(-1/T1/X1/Y1 -1/T/X2/Y2)
      WTH1 = XTCH/XES2
      XWT(11)=  WTK*WTH1
C t-chanel only PLUS interferences
      XTCH= (S**2 +U**2 )/(-X1*Y1*T1)*(1+2*XM2*X1/Y1/T1)
     $     +(S1**2+U1**2)/(-X1*Y1*T1)*(1+2*XM2*Y1/X1/T1)
     $     +(S**2 +U1**2)/(-X2*Y2*T )*(1+2*XM2*X2/Y2/T )
     $     +(S1**2+U**2 )/(-X2*Y2*T )*(1+2*XM2*Y2/X2/T )
      XTCHI= XTCH
     $  +1/(T*T1)*(S**2+U**2+S1**2+U1**2)
     $     *( S/X1/X2 +S1/Y1/Y2 +U/X1/Y2 +U1/X2/Y1)
      XES2= 4*S**2*(-1/T1/X1/Y1 -1/T/X2/Y2)  
      WTH2 = XTCHI/XES2
      XWT(12)=  WTK*WTH2
      IF(KEYSIN.GT.3.OR.KEYSIN.LT.0) THEN
       WRITE(NOUT,*) ' ++++ WRONG KEYSIN ',KEYSIN
       STOP
      ENDIF
      WTM   =  XWT(10+KEYSIN)
      XWT(2)=  WTM
CC]]]]]]]]]]]]
C-------------------------------------------- CONSTRUCT MOMENTA -----
      CFI=DCOS(FI)
      SFI=DSIN(FI)
      QK(4)=XK
      QK(3)=XK*QK3
      QK(2)=XK*(QK2*CFI-QK1*SFI)
      QK(1)=XK*(QK2*SFI+QK1*CFI)
      QP(4)=X
      QP(3)=X*C
      QP(2)=-X*SC*SFI
      QP(1)=X*SC*CFI
      DO 8 I=1,4
    8 QM(I)=-QP(I)-QK(I)
      QM(4)=2.D0+QM(4)
C-------------------------------------------- REJECT W VALUES -------
      WT(4)=WT(4)+WTM
      WT(5)=WT(5)+WTM*WTM
      IF(WTM.LT.0.D0)   WT(11)=WT(11)+1.D0
      IF(WTM.GT.WMAX)   WT(12)=WT(12)+1.D0
      IF(WTM.LT.WT(13)) WT(13)=WTM
      IF(WTM.GT.WT(14)) WT(14)=WTM  
      CALL VARRAN(DRVEC,1)
      RN=DRVEC(1)
C principal weight
      CALL GMONIT(0,72,WTM,WMAX,RN)
      CALL GMONIT(0,70,WTM,WMAX,RN)
C auxiliary weights
      CALL GMONIT(0,51,WTH1*WTK,WMAX,RN)
      CALL GMONIT(0,52,WTH2*WTK,WMAX,RN)
C ...
      IF(KEYWGT.EQ.0.AND.RN*WMAX.GT.WTM) GOTO 6
      IF(KEYWGT.EQ.0) WTM=1D0
      WT(6)=WT(6)+1.D0
C-------------------------------------------- REFLECTION POINT ------
      CALL VARRAN(DRVEC,1)
      RN=DRVEC(1)
CCCCCCCCCC         (==( 106 )==)
C Reflection QP <=> QM can be optionally suspended
      IF(KEYMIR.EQ.1) RN=1
      IF(RN.GT.0.5D0) GOTO 10
      DO 9 I=1,3
      QK(I)=-QK(I)
      QPI=QP(I)
      QP(I)=-QM(I)
    9 QM(I)=-QPI
      QPI=QP(4)
      QP(4)=QM(4)
      QM(4)=QPI
   10 CONTINUE
CCCCCCCCCC      (==( 008 )==)
      DO 50 I=1,4
      QP(I)=EBEAM*QP(I)
      QM(I)=EBEAM*QM(I)
   50 QK(I)=EBEAM*QK(I)
      CALL KINBIS
      RETURN
C-------------------------------------------- END OF HARD PART ------
C-------------------------------------------- SOFT PHOTON PART ------
   11 CONTINUE
      WT(7)=WT(7)+1.D0
C-------------------------------------------- GENERATE C VALUE ------
      CALL VARRAN(DRVEC,1)
      C=1.D0-1.D0/(DMIN+DRVEC(1)*(DMAX-DMIN))
C-------------------------------------------- CALCULATE WEIGHT ------
C Define true transfer (==( 106 )==) 
      tran= 2*ebeam**2*(1-C)
CC[[[[[[[[[       (==(102)==)
      CALL VIRSOF(EBEAM,XK0,C,BORN,CORR) 
      XWT(10)=BORN*(1.D0+CORR)*(1.D0-C)**2*ZSOFT
      DO 52 K=1,3
      CALL VIRSEL(K,EBEAM,XK0,C,WBORN,CORR)
      XWT(10+K) =WBORN*(1.D0+CORR)*(1.D0-C)**2
   52 CONTINUE
      WTM   =  XWT(10+KEYSIN)
      XWT(2)=  WTM
C]]]]]]]]]]
      WT(8)=WT(8)+WTM
      WT(9)=WT(9)+WTM*WTM
      IF(WTM.LT.0.D0)   WT(15)=WT(15)+1.D0
      IF(WTM.GT.WMAX)   WT(16)=WT(16)+1.D0
      IF(WTM.LT.WT(17)) WT(17)=WTM
      IF(WTM.GT.WT(18)) WT(18)=WTM
C-------------------------------------------- REJECT W VALUES -------
      CALL VARRAN(DRVEC,1)
      RN=DRVEC(1)  
C principal weight
      CALL GMONIT(0,71,WTM,WMAX,RN)
      CALL GMONIT(0,70,WTM,WMAX,RN)
C auxiliary weights
      CALL VIRSEL(1,EBEAM,XK0,C,WBORN,CORR)
      WTM1=WBORN*(1.D0+CORR)*(1.D0-C)**2
      CALL VIRSEL(KEYSIN,EBEAM,XK0,C,WBORN,CORR)
      WTM2=WBORN*(1.D0+CORR)*(1.D0-C)**2
      CALL GMONIT(0,41,WTM1,WMAX,RN)
      CALL GMONIT(0,42,WTM2,WMAX,RN)
C  (==( 101 )==)
      IF(KEYWGT.EQ.0.AND.RN*WMAX.GT.WTM) GOTO 6
      IF(KEYWGT.EQ.0) WTM=1D0
      WT(10)=WT(10)+1.D0
C-------------------------------------------- GENERATE FI VALUE -----
      FI=R*TWOPI
C-------------------------------------------- CONSTRUCT MOMENTA -----
      CFI=DCOS(FI)
      SFI=DSIN(FI)
      SC=DSQRT(1.D0-C*C)
      QK(4)=0.D0
      QK(3)=0.D0
      QK(2)=0.D0
      QK(1)=0.D0
      QP(4)=1.D0
      QP(3)=C
      QP(2)=SC*SFI
      QP(1)=SC*CFI
      DO 12 I=1,3
   12 QM(I)=-QP(I)
      QM(4)=1.D0
      DO 15 I=1,4
      QP(I)=EBEAM*QP(I)
   15 QM(I)=EBEAM*QM(I)
      CALL KINBIS
C     ===================================
      ELSEIF(MODE.EQ.1.OR.MODE.EQ.2) THEN
C     ===================================
CCCCCCCCC                       (==( 004 )==)
      NPAR(10)= NEVGEN
      NPAR(20)= NEVGEN
      CALL GMONIT(1,70,AWT70,DWT70,DUMM3)
      CALL GMONIT(1,71,AWT71,DWT71,DUMM3)
      CALL GMONIT(1,72,AWT72,DWT72,DUMM3)
      XST    =  SIGT*AWT70
      DXST   =  XST *DWT70
      XSS    =  SIGS*AWT71
      DXSS   =  XSS *DWT71
      XSH    =  SIGH*AWT72
      DXSH   =  XSH *DWT72
      XST2   =  XSS+XSH
      DXST2  =  SQRT(DXSS**2 +DXSH**2)
      XPAR( 9)= SIG0
C for unweighted events, WTM=1D0
      XPAR(10)=  XST2
      XPAR(11)= DXST2/XST2
      XPAR(12)=  XSS
      XPAR(13)= DXSS    
      XPAR(14)=  XSH
      XPAR(15)= DXSH    
C for WEIGHTED events
      XPAR(20)= SIGS+SIGH
      XPAR(21)= 0D0
      XPAR(22)= SIGS
      XPAR(23)= 0D0
      XPAR(24)= SIGH
      XPAR(25)= 0D0
C for unweighted events, WTM=1D0
      IF(KEYWGT.EQ.0) THEN
        DO 313 I=20,25
  313   XPAR(I)=XPAR(I-10)
      ENDIF      
C ...
      CALL GMONIT(1,41,AWT41,DWT41,DUMM3)
      CALL GMONIT(1,51,AWT51,DWT51,DUMM3)
      XSA    =  SIGS*AWT41
      DXSA   =  XSA *DWT41
      XAH    =  SIGH*AWT51
      DXAH   =  XSH *DWT51
      XSA1   =  XSA+XAH
      DXSA1  =  SQRT(DXSA**2 +DXAH**2) 
      XPAR(41)=  XSA1
      XPAR(51)= DXSA1    
      XPAR(60)=  XSA
      XPAR(70)= DXSA
      XPAR(61)=  XAH
      XPAR(71)= DXAH
      CALL GMONIT(1,42,AWT42,DWT42,DUMM3)
      CALL GMONIT(1,52,AWT52,DWT52,DUMM3)
      XSA    =  SIGS*AWT42
      DXSA   =  XSA *DWT42
      XAH    =  SIGH*AWT52
      DXAH   =  XAH *DWT52
      XSA2   =  XSA+XAH
      DXSA2  =  SQRT(DXSA**2 +DXAH**2)  
      XPAR(42)=  XSA2
      XPAR(52)= DXSA2
      IF(MODE.EQ.1) RETURN      
C     ====================
CCCCCCC        (==( 005 )==)
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '  OLDBIS:        WINDOW A        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '   X.sect. in [nb] units         '
      WRITE(NOUT,BXL1I) NEVGEN,     ' generated events  ','NEVGEN','A0'
      WRITE(NOUT,BXL2F) XSS ,DXSS  ,'Xsec.      soft    ','      ','A1'
      WRITE(NOUT,BXL2F) XSH ,DXSH  ,'Xsec.      hard    ','      ','A2'
      WRITE(NOUT,BXL2F) XST ,DXST  ,'Xsec. straight     ',' total','A3'
      WRITE(NOUT,BXL2F) XST2,DXST2 ,'Xsec. clever       ',' total','A4'
      WRITE(NOUT,BXTXT) '   More on weights etc...        '
      WRITE(NOUT,BXL1F) SIGS ,      'crude Xs. soft     ','      ','A5'
      WRITE(NOUT,BXL1F) SIGH ,      'crude Xs. hard     ','      ','A6'
      DWT71=DWT71*AWT71
      DWT72=DWT72*AWT72
      WRITE(NOUT,BXL2F) AWT71,DWT71,'aver. wt. soft     ','      ','A7'
      WRITE(NOUT,BXL2F) AWT72,DWT72,'aver. wt. hard     ','      ','A8'
      WRITE(NOUT,BXL1F) SIG0 ,      'Born  Xs.          ','      ','A9'
      WRITE(NOUT,BXCLO)  
C ...        
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '  OLDBIS:        WINDOW B        '
      WRITE(NOUT,BXTXT) '     auxiliary information       '
      WRITE(NOUT,BXL2F) XSA1,DXSA1 ,'KSIN=1, no interf. ',' total','B1'
      WRITE(NOUT,BXL2F) XSA2,DXSA2 ,'KSIN=2, with interf',' total','B2'
      WRITE(NOUT,BXCLO)

C old output routine---------
      CALL BABINF
      ENDIF
      END
      SUBROUTINE BABINF
C-------------------------------------------- REMARKS ---------------
C INFO ROUTINE TO BE CALLED AFTER 'BHABHA' HAS RUN EVENTS. IT CALCULA
C THE EXACT CROSS SECTION CORRESPONDING TO THE GENERATED EVENT SAMPLE
C AND SOME STATISTICS ON THE GENERATION.
C WT( 1) = NO.OF STARTS IN HARD-PHOTON PART;
C WT( 2) = NO. OF TRIALS INSIDE C-GENERATION W.R.P. LOOP;
C WT( 3) = NO. OF HARD PHOTON TRIALS SURVIVING C & CT CUTS;
C WT( 4) = SUM OF HARD PHOTON WEIGHTS;
C WT( 5) = SUM OF SQUARED HARD PHOTON WEIGHTS;
C WT( 6) = NO. OF ACCEPTED HARD PHOTON EVENTS;
C WT( 7) = NO. OF STARTS OF SOFT-PHOTON PART;
C WT( 8) = SUM OF SOFT PHOTON WEIGHTS;
C WT( 9) = SUM OF SQUARED SOFT PHOTON WEIGHTS;
C WT(10) = NO. OF ACCEPTED SOFT PHOTON EVENTS;
C WT(11) = NO. OF HARD EVENTS WITH W < 0;
C WT(12) = NO. OF HARD EVENTS WITH W > WMAX;
C WT(13) = MINIMUM GENERATED WEIGHT IN HARD PART;
C WT(14) = MAXIMUM GENERATED WEIGHT IN HARD PART;
C WT(15) = NO. OF SOFT EVENTS WITH W < 0;
C WT(16) = NO. OF SOFT EVENTS WITH W > WMAX;
C WT(17) = MINIMUM GENERATED WEIGHT IN SOFT PART;
C WT(18) = MAXIMUM GENERATED WEIGHT IN HARD PART.
C--------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      COMMON / PARM02 / SIG0,SIGS,SIGH,SIGT,WT
      DIMENSION WT(18)
      COMMON / INOUT  / NINP,NOUT 
      DATA DH/0.D0/,DHE/0.D0/,DDHE/0.D0/,EC/0.D0/
      DATA ECT/0.D0/,EWH/0.D0/,DS/0.D0/,DSE/0.D0/
      DATA DDSE/0.D0/,EWS/0.D0/,DT/0.D0/,DTE/0.D0/,DDTE/0.D0/
      IF(WT(1).EQ.0.D0) GOTO 1
      SIGHE=SIGH*WT(4)/WT(1)
      DSIGHE=SIGH/WT(1)*DSQRT(WT(5)-WT(4)**2/WT(1))
      DH  =SIGH/SIG0
      DHE =SIGHE/SIG0
      DDHE=DSIGHE/SIG0
      EC=WT(1)/WT(2)
      ECT=WT(3)/WT(1)
      EWH=WT(6)/WT(3)
    1 CONTINUE
      IF(WT(7).EQ.0.D0) GOTO 2
      SIGSE=SIGS*WT(8)/WT(7)
      DSIGSE=SIGS/WT(7)*DSQRT(WT(9)-WT(8)**2/WT(7))
      DS  =SIGS/SIG0
      DSE =SIGSE/SIG0
      DDSE=DSIGSE/SIG0
      EWS=WT(10)/WT(7)
    2 CONTINUE
      DT  =SIGT/SIG0
      DTE =DHE+DSE
C&&&& DDTE=DDHE+DDSE
      DDTE=SQRT(DDHE**2+DDSE**2)
      SIGTOT=DTE*SIG0
      SIGSFT=DSE*SIG0
      SIGHRD=DHE*SIG0
C-----
      WRITE(NOUT,3) SIG0,SIGH,SIGS,SIGT,(WT(I),I=1,18)
    3 FORMAT(1H0,90(1H-)/,' BHABHA SAMPLE STATISTICS'/,(4D15.6))
      WRITE(NOUT,1004) SIG0,DH,DHE,DDHE,EC,ECT,EWH,DS,DSE,DDSE
      WRITE(NOUT,1005) EWS,DT,DTE,DDTE,SIGTOT,SIGSFT,SIGHRD
 1004 FORMAT(
     . '          LOWEST ORDER CROSS SECTION =',D15.6,' NB = UNIT'/,
     . '   APPROXIMATED HARD PHOTON XSECTION =',F15.6/,
     . '          EXACT HARD PHOTON XSECTION =',F15.6/,
     . '                         UNCERTAINTY =',F15.6/,
     . ' W.R.P EFFICIENCY IN INTERNAL C LOOP =',F15.6/,
     . '   '      '      OF C/CT RESTRICTION =',F15.6/,
     . '   '      '         FOR HARD WEIGHTS =',F15.6/,
     . '   APPROXIMATED SOFT PHOTON XSECTION =',F15.6/,
     . '          EXACT SOFT PHOTON XSECTION =',F15.6/,
     . '                         UNCERTAINTY =',F15.6)
 1005 FORMAT(
     . '  W.R.P. EFFICIENCY FOR SOFT WEIGHTS =',F15.6/,
     . '    APPROXIMATED TOTAL CROSS SECTION =',F15.6/,
     . '           EXACT TOTAL CROSS SECTION =',F15.6/,
     . '                         UNCERTAINTY =',F15.6/,
     . ' ----------> TOTAL CROSS SECTION =====',D15.6,' NB'/,
     . ' ---------->  SOFT CROSS SECTION =====',D15.6,' NB'/,
     . ' ---------->  HARD CROSS SECTION =====',D15.6,' NB')
      WRITE(NOUT,5) (WT(I),I=11,18)
    5 FORMAT('0 GENERATED WEIGHTS:'/,
     . '          < 0 IN HARD PART =',F15.6/,
     . '       > WMAX IN HARD PART =',F15.6/,
     . '      MINIMUM IN HARD PART =',F15.6/,
     . '      MAXIMUM IN HARD PART =',F15.6/,
     . '          < 0 IN SOFT PART =',F15.6/,
     . '       > WMAX IN SOFT PART =',F15.6/,
     . '      MINIMUM IN SOFT PART =',F15.6/,
     . '      MAXIMUM IN SOFT PART =',F15.6)
      END
      SUBROUTINE VIRSEL (KEY,EBEAM,XK0,COST,WBORN,CORR)
C     ************************************************
C new routine, not present in OLDBAB     (==(102)==)
C......COMPACT FORMULA FOR CROSS SECTION IN BHABHA PROCESS
C......BORN + SOFT +VIRTUAL CORRECTIONS... ONLY QED VERSION
C......MOST OF THE FORMULAS FROM
C......BOHM,DENNIER,HOLLIK,NUCL.PHYS.B304(1988),687
C......COMPLETED E. RICHTER-WAS. APRIL 1990.........
*.....AMGAM IS DUMMY PARAMETR
*.....KEY=1   T NONINTERFERENCE CHANNEL
*.....KEY=2   T CHANNEL
*.....KEY=3 S+T CHANNEL
*     ********************
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
      DATA AMGAM,AMEL /1D0, 0.511D-3/
 
      ALFA=1D0/ALFINV
      ALFPI=ALFA/PI
      AKMAX=XK0*EBEAM
      S=4D0*EBEAM**2
      T=-1/2D0*S*(1D0-COST)
      U=-1/2D0*S*(1D0+COST)
      DLT=DLOG(-T/AMEL**2)
      DLS=DLOG( S/AMEL**2)
      DLU=DLOG(-U/AMEL**2)
      DLTG=DLOG(-T/AMGAM**2)
      DLSG=DLOG( S/AMGAM**2)
      DLKG=DLOG(2D0*AKMAX/AMGAM)
      DLKGM=DLOG(AMEL**2*AKMAX**2/AMGAM**2/EBEAM**2)
      DLUS=DLOG(-U/S)
      DLTS=DLOG(-T/S)
      DLTU=DLOG(T/U)
C........S + T CHANNEL...............
      IF(KEY.EQ.3) THEN
C....BORN LEVEL
      XT1=U/T
      XT3=S/T
      XS1=U/S
      XS2=T/S
      XT=XT1**2+XT3**2
      XS=XS1**2+XS2**2
      XST=2D0*XS1*XT1
      BORN=XT+XS+XST
C....VIRTUAL CORRECTIONS
      VIRT =ALFPI*(
     $      -1D0+PI*PI/12D0+1/4D0*DLT+1/4D0*DLT**2
     $      +1/2D0*DLTG*(1D0-DLT))
      DVIRT=4D0*XT*VIRT
      VIRS =ALFPI*(
     $      -1D0+PI*PI/3D0+1/4D0*DLS+1/4D0*DLS**2
     $      +1/2D0*DLSG*(1D0-DLS))
      DVIRS=4D0*XS*VIRS
      DVIRST  =2D0*(VIRT+VIRS)*XST
C.....BOX CORRECTIONS
      BOXS1= 2D0*ALFPI*(-DLTU*DLSG
     $      +S/2D0/(S+T)*DLTS-S*(S+2D0*T)/4D0/(S+T)**2*DLTS**2)
      BOXS2= 2D0*ALFPI*(-DLTU*DLSG
     $      -S/2D0/(S+U)*DLUS+S*(S+2D0*U)/4D0/(S+U)**2*DLUS**2)
      BOXT1= 2D0*ALFPI*(+DLUS*DLTG
     $      -T/2D0/(T+S)*DLTS-T*(T+2D0*S)/4D0/(T+S)**2*DLTS**2)
      BOXT3= 2D0*ALFPI*(+DLUS*DLTG
     $ +T/2D0/(T+U)*DLTU+T*(T+2D0*U)/4D0/(T+U)**2*(DLTU**2+PI*PI))
      DBOXT=BOXT1*XT1**2+BOXT3*XT3**2
      DBOXS=BOXS1*XS1**2+BOXS2*XS2**2
      DBOXST=1/2D0*(BOXS1+BOXT1)*XST
C.....SOFT CORRECTIONS
      XS=2D0*ALFPI*(
     $      (DLS-1D0)*DLKGM+1/2D0*DLS**2-PI*PI/3D0       )
      XT=2D0*ALFPI*(
     $      (DLT-1D0)*DLKGM+1/2D0*DLT**2
     $      -1/2D0*DLTS**2+PI*PI/6D0
     $      -2D0*DILOGY(1D0+2D0*EBEAM/DSQRT(-T))
     $      -2D0*DILOGY(1D0-2D0*EBEAM/DSQRT(-T))   )
      XU=2D0*ALFPI*(
     $      (DLU-1D0)*DLKGM+1/2D0*DLU**2
     $      -1/2D0*DLUS**2+PI*PI/6D0
     $      -2D0*DILOGY(1D0+2D0*EBEAM/DSQRT(-U))
     $      -2D0*DILOGY(1D0-2D0*EBEAM/DSQRT(-U))   )
      SOFT=XS+XT-XU
C.....ADDED TOGETHER
      CORS=DVIRS+DBOXS
      CORT=DVIRT+DBOXT
      CORST=DVIRST+DBOXST
      CORSOFT=BORN*SOFT
      BORNCOR=BORN+CORSOFT+CORS+CORT+CORST
      CORR=(BORNCOR-BORN)/BORN
      ENDIF
C.....end of...S + T CHANNEL...end of............
C........ T CHANNEL...............
      IF(KEY.EQ.2) THEN
C....BORN LEVEL
      XT1=U/T
      XT3=S/T
      XT=XT1**2+XT3**2
      BORN=XT
C....VIRTUAL CORRECTIONS
      VIRT =ALFPI*(
     $      -1D0+PI*PI/12D0+1/4D0*DLT+1/4D0*DLT**2
     $      +1/2D0*DLTG*(1D0-DLT))
      DVIRT=4D0*XT*VIRT
C.....BOX CORRECTIONS
      BOXT1= 2D0*ALFPI*(+DLUS*DLTG
     $      -T/2D0/(T+S)*DLTS-T*(T+2D0*S)/4D0/(T+S)**2*DLTS**2)
      BOXT3= 2D0*ALFPI*(+DLUS*DLTG
     $ +T/2D0/(T+U)*DLTU+T*(T+2D0*U)/4D0/(T+U)**2*(DLTU**2+PI*PI))
      DBOXT=BOXT1*XT1**2+BOXT3*XT3**2
C.....SOFT CORRECTIONS
      XS=2D0*ALFPI*(
     $      (DLS-1D0)*DLKGM+1/2D0*DLS**2-PI*PI/3D0       )
      XT=2D0*ALFPI*(
     $      (DLT-1D0)*DLKGM+1/2D0*DLT**2
     $      -1/2D0*DLTS**2+PI*PI/6D0
     $      -2D0*DILOGY(1D0+2D0*EBEAM/DSQRT(-T))
     $      -2D0*DILOGY(1D0-2D0*EBEAM/DSQRT(-T))   )
      XU=2D0*ALFPI*(
     $      (DLU-1D0)*DLKGM+1/2D0*DLU**2
     $      -1/2D0*DLUS**2+PI*PI/6D0
     $      -2D0*DILOGY(1D0+2D0*EBEAM/DSQRT(-U))
     $      -2D0*DILOGY(1D0-2D0*EBEAM/DSQRT(-U))   )
      SOFT=XS+XT-XU
C.....ADDED TOGETHER
      CORT=DVIRT+DBOXT
      CORSOFT=BORN*SOFT
      BORNCOR=BORN+CORSOFT+CORT
      CORR=(BORNCOR-BORN)/BORN
      ENDIF
C....end of....T  CHANNEL........end of........
C........ T non interference CHANNEL...............
      IF(KEY.EQ.1) THEN
C....BORN LEVEL
      XT1=U/T
      XT3=S/T
      XT=XT1**2+XT3**2
      BORN=XT
C....VIRTUAL CORRECTIONS
      VIRT =ALFPI*(
     $      -1D0+PI*PI/12D0+1/4D0*DLT+1/4D0*DLT**2
     $      +1/2D0*DLTG*(1D0-DLT))
C.....SOFT CORRECTIONS
      SOFT=2D0*ALFPI*(
     $      (DLT-1D0)*DLKGM+1/2D0*DLT**2
     $      -1/2D0*DLTS**2+PI*PI/6D0
     $      -2D0*DILOGY(1D0+2D0*EBEAM/DSQRT(-T))
     $      -2D0*DILOGY(1D0-2D0*EBEAM/DSQRT(-T))   )
C.....ADDED TOGETHER
      CORR=SOFT+4D0*VIRT
      ENDIF
C....end of....T non interference  CHANNEL........end of........   
      WBORN = BORN/8D0
      END
 
      SUBROUTINE VIRSOF(EB,XK0,X,BORN,CORR)
C-------------------------------------------- REMARKS ---------------
C BHABHA SCATTERING DIFFERENTIAL CROSS SECTION WITH CORRECTIONS:
C 1) SELF-ENERGY DIAGRAMS FOR IN- AND OUTGOING LEPTONS;
C 2) VERTEX CORRECTION DIAGRAMS;
C 3) BOX DIAGRAMS (TWO-PHOTON EXCHANGE);
C 4) SOFT BREMSSTRAHLUNG (PHOTON ENERGY < XK0*EBEAM );
C 5) VACUUM POLARIZATION (PHOTON SELF-ENERGY DIAGRAMS);
C 6) INTERFERENCE BETWEEN PHOTON AND Z0 EXCHANGE GRAPHS.
C FORMULA AND CONVENTIONS TAKEN FROM:
C F.A.BERENDS ET AL, NUCL.PHYS.B68(1974)541.
C EB   = BEAM ENERGY IN GEV;
C XK0  = CUTOFF ON SOFT BREMSSTRAHLUNG ENERGY;
C X    = COSINE OF POLAR SCATTERING ANGLE OF POSITRON;
C BORN = LOWEST-ORDER DIFFERENTIAL CROSS SECTION IN NANOBARN;
C CORR = TOTAL OF CORRECTIONS GIVEN ABOVE.
C ASSUMED VALUES: 90 GEV FOR THE Z0 MASS, .23 FOR SIN**2(TH).
C--------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      COMMON / INOUT  / NINP,NOUT 
C-------------------------------------------- LOWEST ORDER ----------
      X2=X*X
      X3=X2*X
      X4=X3*X
      XNUM=9.D0+6.D0*X2+X4
      BORN=1.295958D0/EB**2*XNUM/(1.-X)**2
C-------------------------------------------- CORRECTIONS (NO VAC.POL
      U=2.D0*DLOG(3.9139D03*EB)
      V=U+DLOG(1.+X)-0.6931472D0
      W=U+DLOG(1.-X)-0.6931472D0
      CORR=2.3228D-03*(
     .-4.D0*(1.D0-U+V-W)*DLOG(XK0) - 6.5797D0 + U*U - V*V + W*W
     .+ 2.*DILOGY((1.D0+X)/2.D0) - 2.*DILOGY((1.D0-X)/2.D0)
     .+ ( ( 1.-12.*X+12.*X2 -4.*X3 +3.*X4)*U
     .   -( 5. -7.*X +3.*X2    -X3       )*V
     .   +(31. +5.*X +9.*X2 +3.*X3       )*W
     .   +( 3. +7.*X -5.*X2 -3.*X3 -2.*X4)*U*U*.5
     .   +( 3. -3.*X    +X2    -X3       )*V*V
     .   -( 9. +7.*X+11.*X2 +5.*X3       )*W*W*.5
     .   -( 2.    -X           -X3       )*U*V*X*2.
     .   -(21. +3.*X +9.*X2 -3.*X3 +2.*X4)*U*W
     .   +( 6. +5.*X +4.*X2    +X3       )*V*W*2.
     .   -(36.      +24.*X2        +4.*X4)
     .   +(18.-15.*X+12.*X2 -3.*X3 +4.*X4)*3.2899)/XNUM)
C<<<<<--------------------------------------- VACUUM POLARIZATION ---
C[[[[[       (==( 103 )=)
      CORR=CORR+
     .(6.*X -6.*X2 +2.*X3 -2.*X4)*REPI(4.*EB*EB)/XNUM+
     .(-18.-6.*X-6.*X2-2.*X3)*REPI(-2.*EB*EB*(1.-X))/XNUM
C<<<<---------------------------------------- INTERFERENCE WITH Z0 --
C&&&  CORR=CORR + XSWEAK(EB,9.D01,.23D0,X)/BORN
CCCC  CORR=CORR + XSWEAK(EB,92D0 ,0.2288D0,X)/BORN
CCCCCCCCCCC (==( 104 )==)
CCCC  Note outdated Z mass, width and SINW2 here !!!!
      CORR=CORR + YYWEAK(EB,92D0,2.4533D0,0.2288D0,X)
      Y=BORN*(1.+CORR)*1000.D0
CCC   IF(CORR.LT.-1.D0) WRITE(NOUT,1) CORR,EB,XK0,X
CCC 1 FORMAT(' ***VIRSOF WARNING : ',4D15.6)
      END    
CCCCCCCCCCCCCCCC (==( 104 )==)
C slightly improved version of XSWEAK with finite width of Z resonance
CCCCCCCCCCCCCCCC
C BHABHA CROSS SECTION DUE TO INTERFERENCE OF THE PHOTON AND Z0 GRAPH
C EB  = BEAM ENERGY (GEV)
C XMZ0= Z0 MASS (GEV)
C GAMZ0= Z0 WIDTH (GEV)
C SIN2= SIN(WEAK MIXING ANGLE)**2
C X   = COSINE OF POLAR SCATTERING ANGLE OF POSITRON
      FUNCTION YYWEAK(EB,XMZ0,GAMZ0,SIN2,X)
*     *************************************
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      XM=1.-X
      A2=1./(16.*SIN2*(1.-SIN2))
      VP=((1.-4.*SIN2)**2+1.)*A2
      VM=VP-2.*A2
      SS=4*EB**2
      TT= SS*(1-X)/2D0
C (==( 104 )==)
CC    XS=SS/(SS-XMZ0**2 )
      XS=SS*(SS-XMZ0**2)/((SS-XMZ0**2)**2+(XMZ0*GAMZ0)**2)
      XT=TT/(TT+XMZ0**2 )
      XP=(1.+X)**2
CC    FACT=5.183833D0/EB**2
CC    BORN=1.295958D0/EB**2*XNUM/(1.-X)**2
      XNUM=9.D0+6.D0*X**2+X**4
      FACT= 4D0/XNUM*XM**2
      YYWEAK=FACT*( XS*( VP*XP + VM*XM**2 )/4.
     $            + XT*( VP*XP + VM*4.    )/(XM*XM)
     $        -(XS+XT)*( VP*XP            )/(2.*XM) )
      END
C==================== not used anymore...
      FUNCTION XSWEAK(EB,XMZ0,SIN2,X)
C-------------------------------------------- REMARKS ---------------
C BHABHA CROSS SECTION DUE TO INTERFERENCE OF THE PHOTON AND Z0 GRAPH
C EB  = BEAM ENERGY (GEV)
C XMZ0= Z0 MASS (GEV)
C SIN2= SIN(WEAK MIXING ANGLE)**2
C X   = COSINE OF POLAR SCATTERING ANGLE OF POSITRON
C THE COUPLING CONSTANTS OF THE ELECTRONS TO THE Z0 ARE CALCULATED
C ACCORDING TO THE STANDARD SU(2)*U(1) MODEL, USING SIN2. THE MASS
C XMZ0 IS TREATED AS AN ADDITIONAL, INDEPENDENT PARAMETER.
C NEITHER THE PURE Z0 CHANNEL, NOR THE EFFECTS OF A NONZERO Z0 WIDTH
C ARE TAKEN INTO ACCOUNT ---> THIS ROUTINE IS NOT GOOD FOR LEP/SLC.
C--------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      DATA INIT /0/
      IF(INIT.NE.0) GOTO 1
      INIT=1
C-------------------------------------------- INITIALIZATION --------
      A2=1./(16.*SIN2*(1.-SIN2))
      VP=((1.-4.*SIN2)**2+1.)*A2
      VM=VP-2.*A2
      CHIQ=(XMZ0/EB)**2/2.
      XS=1./(1.-CHIQ/2.)
      FACT=5.183833D0/EB**2
C--------------------------------------------------------------------
    1 XM=1.-X
      XP=(1.+X)**2
      XT=1./(1.+CHIQ/XM)
      XSWEAK=FACT*( XS*( VP*XP + VM*XM**2 )/4.
     .            + XT*( VP*XP + VM*4.    )/(XM*XM)
     .        -(XS+XT)*( VP*XP            )/(2.*XM) )
      RETURN
      END

      SUBROUTINE KINBIS       
C     *****************
C TRANSTATES /MOMBAB/ INTO /MOMSET/     
C     ************************
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
C from OLDBIS
      COMMON / PARM01 / EBEAM,THMIN,THMAX,XKMIN,XKMAX,CMIN,CMAX
      COMMON / MOMBAB / QP(4),QM(4),QK(4),WTM,XWT(20)     
C to BHLUMI
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300) 
      COMMON / TRANSR / TRAN,TRANP,TRANQ
      AMEL  =  0.5111D-3
C beams
      P1(4)= EBEAM
      P1(3)= DSQRT(EBEAM**2-AMEL**2)
      P1(2)= 0D0
      P1(1)= 0D0
C...
      Q1(4)= EBEAM
      Q1(3)=-DSQRT(EBEAM**2-AMEL**2)
      Q1(2)= 0D0
      Q1(1)= 0D0
C final electrons
      DO 10 K=1,4
      P2(K)=QP(K)
      Q2(K)=QM(K)
   10 CONTINUE      
C photon, if present  
      NPHOT=0
      IF(QK(4).NE.0D0) THEN
       NPHOT=1
       DO 20 K=1,4
   20  PHOT(1,K)=QK(K)
      ENDIF
C Transfers    
      tranp= 2*(p1(4)*p2(4)-p1(3)*p2(3)-p1(2)*p2(2)-p1(1)*p2(1))
      tranq= 2*(q1(4)*q2(4)-q1(3)*q2(3)-q1(2)*q2(2)-q1(1)*q2(1))
C weights
      WTCRU1=1D0
      WTCRU2=1D0
      WTMOD    = WTM
      DO 70 K=1,20
   70 WTSET(K) = XWT(K)                   
C----------------------------------------------------------------C
C                      The end of OLDBIS                         C
C----------------------------------------------------------------C
      END      
         

