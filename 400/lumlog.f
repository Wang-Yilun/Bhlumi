      SUBROUTINE LUMLOG(MODE,XPAR,NPAR)
!     *********************************  
! =================================================================== 
! ===================================================================
! ===BBB=======BBB==BBB===BBB===BBB===BBB=======BBBBBB=====BBBBBBB===
! ===BBB=======BBB==BBB===BBBB=BBBB===BBB======BBBBBBBB===BBBBBBBBB==
! ===BBB=======BBB==BBB===BBBBBBBBB===BBB======BBB==BBB===BBB========
! ===BBB=======BBB==BBB===BBBBBBBBB===BBB======BBB==BBB===BBB==BBBB==
! ===BBBBBBBB==BBBBBBBB===BBB=B=BBB===BBBBBBB==BBB==BBB===BBB===BBB==
! ===BBBBBBBB===BBBBBB====BBB===BBB===BBBBBBB===BBBBBB=====BBBBBB====
! ===================================================================
!
!           ************************************************** 
!           *       **********************************       *
!           *       *      *******************       *       *
!           *       *      *                 *       *       *
!           *       *      *   L U M L O G   *       *       *
!           *       *      *                 *       *       *
!           *       *      *******************       *       *
!           *       **********************************       *
!           **************************************************    
!
! ----------------------------------------------------------------C
!                                                                C
!                       LUMLOG                                   C
!                                                                C
!          COLLLINEAR LEADING-LOG MONTE CARLO FOR                C
!               LOW-ANGLE BHABHA SCATTERING                      C
!                     NOVEMBER 1990                              C
!                 last update  5 feb. 91 (pairs)                 C
!                 last update 14 feb. 91 born in robol6          C
!                                     bug in modelu              C
!                 last update  8 apr. 91 cosmetics (sj)          C
!                 last update 26 aug. 91 cosmetics (sj)          C
!                 last update    May. 94 Final state rad. (sj)   C
!                 last update   June. 95 LL emulation     (sj)   C
!         AUTHORS:                                               C
!         S. Jadach, E. Richter-Was, Z. Was, B.F.L. Ward         C
!                                                                C
! The user is kindly requested to cite preprint TH.5995 (1991)   C
! Phys. Lett. B260 (1991) 438 of the same authors.               C
! (Note that TH.5888 is now Phys. Lett. B253 (1991) 469).        C
! ----------------------------------------------------------------C
!
! Note that LUMLOG contains originally two MC program: BHALOG which
! is present here and MULTILOG which is not inluded here.
! The series of comparisons between these two programs (TH-5995)
! has leaded the 0.02% technical precision estimate for LUMLOG.
!
! ----------------------------------------------------------------------
!                 INPUT and OUTPUT of LUMLOG
! ----------------------------------------------------------------------
! All input and output goes through parameters in 
!                 CALL LUMLOG(MODE,XPAR,NPAR)
! and through /MOMSET/ and /WGTALL/ common blocks.
! In the following we shall  briefly indicate the meaning of the
! above parameters/variables.
!
! IF( MODE =-1 ) THEN
! ===================
! Initialisation is performed, all input parameters are transfered
! through XPAR and NPAR.
! In the following table we indicate the meaning of NPAR, XPAR 
! entries for LUMLOG subgenerator in the initialization mode.
!         Table        Input parameters of LUMLOG
! ----------------------------------------------------------------------
!  Entry    Variable   Meaning
! ----------------------------------------------------------------------
!  NPAR( 1)  KEYOPT =1000*KEYGEN +10*KEYWGT +KEYRND 
!                   General option switch 
!            KEYGEN =2 for this sub-generator
!            KEYRND =1,2 type of random number generator RANMAR,RANECU
!            KEYWGT =1 for variable weight WTM
!            KEYWGT =2 generation down to zero angle
!  NPAR( 2)  KEYRAD =100*KEYFIN +10*KEYTES+ KEYBLO , QED option switch 
!            KEYFIN =  0    No explicit final state emission
!            KEYFIN =  1    Final state emission included
!            KEYTES =  0    normal str. function DEFAUT!!!
!            KEYTES =  1    test str. functions (1-z)**(-1/2)
!            KEYBLO =  3    LLog  ln(s'*xi_dot/me**2)  -1)  DEFAUT!!!
!            KEYBLO =  4    LLog  ln(s*xiA/me**2)   
!  XPAR( 1)  CMSENE Total center mass energy [GeV]
!  XPAR( 2)   TMINL Theta_minimum [degr]
!  XPAR( 3)   TMAXL Theta_maximum [degr]
!  XPAR( 4)     XK0 Dimensionless infrared cut-off (on real soft 
!                   photons) relevant for un-exponentiated case.
!                   Range 0.000001<XKO<0.0001 recommeded.
!  XPAR( 5)   XKMAX Determines minimum effective mass s' of the
!                   final state electron-positron, s'>s*(1-XKMAX).
!                   XKMAX=1 is alowed and recommended.
! ----------------------------------------------------------------------
!
! ELSE IF( MODE = 0 ) THEN
! ========================
! Generation of the single Monte Carlo event. 
! Only variable weight events are produced (!)
! (the user may turn them, if he wishes, into WT=1 events by himself).
! The four momenta of the final state electron, positron and photon
! are encoded in 
!      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
! where P1 and Q1 are four-momenta of positron and elecron beams.
! P2 and Q2 are four-momenta of outgoing positron and electron.
!                       IMPORTANT 
!                       ---------
! Remember that P2 and Q2 from LUMLOG do not represent normal 'bare'
! electron and positron but 'dressed' ones, i.e. they sum four-momenta
! of electron and ALL photons collinear with them 
! (as calorimeter does). 
! (N.B. Collinearity relation extends to photons with transverse 
! momentum up to characteristic LL scale pT_max = sqrt(Q**2)).
!                       ---------
! As a result, the list PHOT(100,4) of photon four-momenta is empty 
! and NPHOT=0 (the number of real photons in PHOT).  
! The principal weight WTM of the event is placed in
!      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)   
! It is usually of interest to use 'paralel weights' from WTSET.
! The event weight is constructed then as WT= WTCRU1*WTCRU2*WTSET(J).
! Which J is alowed and what version of the QED matrix element 
! it represents is summarized in the table below.
! N.B. principal weight WTM= WTCRU1*WTCRU2*WTSET(4).
! All calculation in LUMLOG are Leading-Log initial state bremss. type.
!       Table of WTSET entries for LUMLOG
! ----------------------------------------------------------------------
!  Entry      Type of the used electron (non-singl) structure function
! ----------------------------------------------------------------------
!             QED order   Exponentiation       pairs 
!             ---------------------------------------------------------
!  WTSET( 1)  Zero-th     exponentiated        No
!  WTSET( 2)  First       exponentiated        No
!  WTSET( 3)  Second      exponentiated        No
!  WTSET( 4)  Third       exponentiated        No <<== PRINCIPAL WEIGHT
!  WTSET( 5)  Third       exponentiated       Yes
!  WTSET(11)  Zero-th     not exponentiated    No
!  WTSET(12)  First       not exponentiated    No
!  WTSET(13)  Second      not exponentiated    No
!  --------------------------------------------------------------------
!     LL emulation of the full matrix element (June 95)
!  --------------------------------------------------------------------
! All beta's:  O(alf0),O(alf1),O(alf2)
!      WTSET( 40) =   Total O(alf0)
!      WTSET( 41) =   Total O(alf1)
!      WTSET( 42) =   Total O(alf2) <<== PRINCIPAL emulation of multiph.
! Individual beta's in various orders.
! O(alf1)
!      WTSET( 43) =     beta0
!      WTSET( 44) =     beta1
!      WTSET( 45) =     beta1 upper  line component
!      WTSET( 46) =     beta1 lower  line component
! O(alf2)
!      WTSET( 47) =     beta0
!      WTSET( 48) =     beta1
!      WTSET( 49) =     beta2
!      WTSET( 50) =     bt11u =  beta1 upper line component
!      WTSET( 51) =     bt11l =  beta1 lower line component
!      WTSET( 52) =     bt2ul =  beta2 upper*lower component
!      WTSET( 53) =     bt20u =  beta2 upper component
!      WTSET( 54) =     bt20l =  beta2 lower component
!  --------------------------------------------------------------------
! 
! ELSE IF( MODE = 1 ) THEN
! ========================
! The total cross section corresponding to generated series of event,
! i.e. resulting from MC integrartion is calculated and stored in XPAR
! and NPAR, see table below.
!  --------------------------------------------------------------------
!  Entry    Variable   Meaning
!  --------------------------------------------------------------------
!  NPAR(20)  NEVGEN  Number of generated MC events
!  XPAR(20)    XCRU  Crude total MC x-section [nb] which is necessary
!                    for rescaling histograms 
!                    in run with weighted events.
!  XPAR(21)          =0, error of XPAR(20), it is identicaly zero
!  XPAR(22)    BORN  Born x-cextion [nb]
!  XPAR(25)    SIG0  Miscelanous
!  --------------------------------------------------------------------
! For MODE=1 program is called upon many times in the process of 
! rescaling histograms, therefore, there is no output printed in 
! this mode.
!
! ELSE IF( MODE = 2 ) THEN
! ========================                     
! 
! Only in this MODE=2 in addition to filling XPAR and NPAR as 
! for MODE=1 the values of various x-sections are printed on 
! the standard output file.
!                
! ENDIF
! ====
!
!     ***********************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      DIMENSION XPAR(*),NPAR(*)
! 
      CALL BHALOG(MODE,XPAR,NPAR)
      END

      SUBROUTINE BHALOG(MODE,XPAR,NPAR)
!     ***********************************  
! ----------------------------------------------------------------C
!            BHALOG is part of LUMLOG library                    C
! ----------------------------------------------------------------C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(ALFA   =  1D0/137.03604D0,   ALF1=  ALFA/PI)
      PARAMETER(GNANOB =  389385D0) 
      PARAMETER(AMEL   =  0.511D-3)       
      PARAMETER(X0     =  1D-5)       
      DIMENSION XPAR(*),NPAR(*)
      COMMON / INOUT  / NINP,NOUT 
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      COMMON / BHPARZ / AMAZ,GAMMZ,SINW2,GV,GA
      COMMON / RANPAR / KEYRND
      LOGICAL LCONE1,LCONE2,LCONE,LENERG     
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G
      DOUBLE PRECISION DRVEC(100)
      EXTERNAL FUNSKI

      IF( MODE.EQ.-1) THEN
!     ********************    
! ...BX-formats for nice and flexible outputs
      BXOPE =  '(//1X,15(5H*****)    )'
      BXTXT =  '(1X,1H*,                  A48,25X,    1H*)'
      BXL1I =  '(1X,1H*,I17,                 16X, A20,A12,A7, 1X,1H*)'
      BXL1F =  '(1X,1H*,F17.8,               16X, A20,A12,A7, 1X,1H*)'
      BXL2F =  '(1X,1H*,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H*)'
      BXL1G =  '(1X,1H*,G17.8,               16X, A20,A12,A7, 1X,1H*)'
      BXL2G =  '(1X,1H*,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H*)'
      BXCLO =  '(1X,15(5H*****)/   )'    
! .........
      CMSENE = XPAR(1)
      TMINL  = XPAR(2)
      TMAXL  = XPAR(3)  
      XK0    = XPAR(4)
      XKMAX  = XPAR(5)
      KEYOPT = NPAR(1)
      KEYRAD = NPAR(2)
!
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '*     *****************         *'
      WRITE(NOUT,BXTXT) '*     ***   LUMLOG  ***         *'
      WRITE(NOUT,BXTXT) '*     *****************         *'
      WRITE(NOUT,BXTXT) '*(P.L. B260 (1991) 438, TH-5995)*'
      WRITE(NOUT,BXTXT) '* This program is now part of   *'
      WRITE(NOUT,BXTXT) '*   BHLUMI version 4.04         *'
      WRITE(NOUT,BXTXT) '*   September      1996         *'
      WRITE(NOUT,BXTXT) '*         AUTHORS               *'
      WRITE(NOUT,BXTXT) '* S. Jadach, E. Richter-Was     *'
      WRITE(NOUT,BXTXT) '*    B.F.L. Ward, Z. Was        *'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)
!
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '   ===== LUMLOG/BHALOG  ======   '
      WRITE(NOUT,BXTXT) '   initialization starts....     '
! .........        
      KEYBLO = MOD(KEYRAD,10) 
      KEYTES = MOD(KEYRAD,100)/10    
      KEYFIN = MOD(KEYRAD,1000)/100    
      KEYRND = MOD(KEYOPT,10)   
      KEYWGT = MOD(KEYOPT,100)/10   
      IF(TMAXL.GE. 180D0) TMAXL=179.9999999D0
      XIA    = (1-COS(TMINL*PI/180))/2  
      XIB    = (1-COS(TMAXL*PI/180))/2  
      RAXI   = XIB/XIA
      DSIG0 = 4D0*ALFA**2*PI/CMSENE**2*GNANOB    
!  ETA is a dummy parameter present only in crude str. funct.
      ETA  = ALF1*DLOG((CMSENE/AMEL)**2*XIA)
! for callibration test
      IF(KEYTES.EQ.1)  ETA  = 0.5D0    
      ZMIN = AMEL/CMSENE   
      IF(XKMAX.LT.1D0) ZMIN = 1D0-XKMAX
!-----------------------------------------------------------------
      IF(KeyWgt .EQ. 1) THEN
!-- OLD method: no generation of transfer below Xia
        bornc = dsig0*(1/xia-1/xib)      
      ELSEIF( KeyWgt .EQ. 2 .AND. abs(Xib/Xia-1).GT.0.1d0 ) THEN
!-- NEW feature: Generation of transfer below Xia
        bornc = dsig0*(2/xia-1/xib)      
      ELSE
        write(6,*) '+++++LUMLOG: wrong KeyWgt= ',KeyWgt
      ENDIF
!-----------------------------------------------------------------
      BORN  = DSIG0*FOBOR(XIA,XIB)
!... Z correction   
      AMAZ  = 91.187d0
      GAMMZ =  2.490d0
      SINW2 = 0.2319d0
      GA = -1/(4*DSQRT(SINW2*(1-SINW2)))
      GV = GA*(1-4*SINW2)

      WRITE(NOUT,BXL1I) KEYOPT,     ' options   switch  ','KEYOPT','N1'
      WRITE(NOUT,BXL1I) KEYRND,     ' VARRAN    switch  ','KEYRND','  '
      WRITE(NOUT,BXL1I) KEYWGT,     ' weight    switch  ','KEYWGT','  '
      WRITE(NOUT,BXL1I) KEYRAD,     ' radiation switch  ','KEYRAD','N2'
      WRITE(NOUT,BXL1I) KEYFIN,     ' KEYFIN            ','KEYFIN','  '
      WRITE(NOUT,BXL1I) KEYTES,     ' KEYTES            ','KEYTES','  '
      WRITE(NOUT,BXL1I) KEYBLO,     ' KEYBLO            ','KEYBLO','  '
      WRITE(NOUT,BXL1I) KEYRND,     ' KEYRND            ','KEYRND','  '
      WRITE(NOUT,BXL1F) CMSENE,     ' CMSENE [GeV]      ','CMSENE','X1'
      WRITE(NOUT,BXL1F) TMINL ,     ' TMIN  [degr.]     ','TMIN  ','X2'
      WRITE(NOUT,BXL1F) TMAXL ,     ' TMAX  [degr.]     ','TMAX  ','X3'
      WRITE(NOUT,BXL1F) XK0   ,     ' xk0 cut (no exp.) ','XK0   ','X4'
      WRITE(NOUT,BXL1F) XKMAX ,     ' minimum sprim/s   ','XKMAX ','X5'
      WRITE(NOUT,BXL1F) RAXI  ,     ' RAXI              ','RAXI  ','  '
      WRITE(NOUT,BXL1F) XIA   ,     ' XIA=(1-cosTMINL)/2','XIA   ','  '
      WRITE(NOUT,BXL1F) XIB   ,     ' XIB=(1-cosTMAXL)/2','XIB   ','  '
      WRITE(NOUT,BXL1F) ETA   ,     ' ETA               ','      ','  '
      WRITE(NOUT,BXL1F) ZMIN  ,     ' ZMIN              ','      ','  '
      WRITE(NOUT,BXL1F) BORNC ,     ' Born crude        ','BORNC ','  '
      WRITE(NOUT,BXL1F) BORN  ,     ' Born exact        ','BORNX ','  '
      WRITE(NOUT,BXL1F) AMAZ  ,     ' Z-mass GeV        ','AMAZ  ','  '
      WRITE(NOUT,BXL1F) GAMMZ ,     ' Z-width GeV       ','GAMMZ ','  '
      WRITE(NOUT,BXL1F) SINW2 ,     ' weak mixing angle ','SINW2 ','  '
! .......... 
      CALL VESK2W(-1,FUNSKI,DUM1,DUM2,WT) 
      IDA=100
      DO 12 K=1,5   
      CALL GMONIT(-1,IDA+   K,0D0,1D0,1D0)  
      CALL GMONIT(-1,IDA+10+K,0D0,1D0,1D0)  
   12 CONTINUE

      WRITE(NOUT,BXTXT) '  end of initialization          '
      WRITE(NOUT,BXCLO)  
!-------------------------------------------------------------
! This is Generator Identificator
      IDGEN =  2        
! Important histo which remembers total x-section
      CALL VESK2W(0,FUNSKI,DUM1,DUM2,WTVES) ! to prevent 1/zero in next call
      CALL VESK2W(1,FUNSKI,AWT,EREL,ZCRUD)
      SIGTZ  = ZCRUD *BORNC
      CALL GMONIT(  -1, IDGEN,0D0,2*SIGTZ,1D0)
!-------------------------------------------------------------
      nevgen=0            
          
!     *********************** 
      ELSEIF( MODE.EQ.0) THEN
!     ***********************     
      nevgen=nevgen+1   
      CALL GMONIT(  0, IDGEN, SIGTZ, 2*SIGTZ,1D0)
!==============================================================
!   Crude distributions, generation of internal MC variables
!==============================================================
  200 CONTINUE
! generation of (x1,x2)      
      CALL VESK2W(0,FUNSKI,DUM1,DUM2,WTVES)       
      XMAX=  1-ZMIN   
      XX1=0D0
      XX2=0D0
      IF(TT1.GT.X0**ETA) XX1 = XMAX*DEXP(1/ETA*DLOG(TT1))
      IF(TT2.GT.X0**ETA) XX2 = XMAX*DEXP(1/ETA*DLOG(TT2))  
      X1 = XX1
      X2 = XX2
      Z1 = 1D0-XX1
      Z2 = 1D0-XX2 
!--------------------------------------------------------------------
! Generate t-channel transfer (the true one)
!--------------------------------------------------------------------
      CALL VARRAN(DRVEC,1)
      rn1=DRVEC(1)
      IF(KeyWgt .EQ. 1) THEN
!--------------------------------------------------------------------
!-- OLD method: no generation of transfer below Xia
        Xicom =  xia*xib/(xib*(1d0-rn1)+XIA*rn1)
        WTxi=1d0
      ELSEIF( KeyWgt .EQ. 2 .AND. abs(Xib/Xia-1).GT.0.1d0 ) THEN
!--------------------------------------------------------------------
!-- NEW feature: Generation of angle below theta minimum
        ynorm= 2/Xia -1/Xib
        prob1= (1/Xia) /ynorm
        prob2= (1/Xia -1/Xib)/ynorm
        IF(rn1.LT.prob1) THEN
!-- flat distribution below Xia
          Xicom= Xia*rn1/prob1
          WTxi = Xia**2/Xicom**2
!-- to prevent crash of kinematics  (1d-8 effect)
          IF(Xicom*CMSene**2.lt.(10*Amel)**2) GOTO 200
        ELSE
!-- normal 1/Xi**2 distribution above Xia
          rnx  = (rn1-prob1)/prob2
          Xicom = 1d0 / (  rnx/Xia   +(1d0-rnx)/Xib )
          WTxi  = 1d0
        ENDIF
      ELSE
        write(6,*) '+++++LUMLOG: wrong KeyWgt= ',KeyWgt
      ENDIF
!--------------------------------------------------------------------
! Translate Xicom into one of Xi1 or Xi2 
      IF(Z1.LE.Z2) THEN 
! generation of xi2 according to 1/xi**2 distribution  
         XI2 =  Xicom
! calculate the other xi1
         XI  =  XI2*Z2/(Z1+XI2*(Z2-Z1))
         XI1 =  XI *Z2/(Z1*(1-XI) +Z2*XI) 
! for callibration test
         IF(KEYTES.EQ.1) XI1 = XI2*(Z2/Z1)**2
      ELSE            
! generation of xi1 according to 1/xi**2 distribution   
         XI1 =  Xicom
! calculate the other xi2
         XI  =  XI1*Z1/(Z2+XI1*(Z1-Z2))
         XI2 =  XI *Z1/(Z2*(1-XI) +Z1*XI)   
! for callibration test
         IF(KEYTES.EQ.1) XI2 = XI1*(Z1/Z2)**2
      ENDIF         
!------------------------------------------------------------------
! Final state bremsstrahlung, crude distributions
!------------------------------------------------------------------
      x3=0
      x4=0
      wtfin3=1
      wtfin4=1
      IF(keyfin.eq.1) THEN
        CALL Brelos(eta,tt3,x3,wtfin3)
        CALL Brelos(eta,tt4,x4,wtfin4)
      ENDIF
!------------------------------------------------------------------
!  Angular Trigger
!------------------------------------------------------------------
! symmetric angular trigger   
      LCONE1 = XI1.GT.XIA.AND.XI1.LT.XIB
      LCONE2 = XI2.GT.XIA.AND.XI2.LT.XIB
      LCONE  = LCONE1.AND.LCONE2       
      LENERG = Z1*Z2 .GT. ZMIN  
!CC No cut on z1*z2 = 1-v any more (stj sept. 91)
!      WTRIG=0D0
!      IF(LCONE) WTRIG=1D0
!CC No angular cut any more (stj may 94)
!      WTRIG=1d0
!
      WTCRU1 = WTVES
      WTCRU2 = WTxi
      WTCRUD = WTCRU1*WTCRU2
!------------------------------------------------------------------
! calculate four-momenta for accepted events
!------------------------------------------------------------------
!  No angular cut any more (stj may 94)
!*******IF( LCONE ) CALL KINOLT
        CALL KINOLT
!==================================================================
!==================    MODEL Weights    ===========================
!==================================================================

      IF(keyfin.eq.0) THEN
! Initial state only
        CALL mdline
        CALL mdlinu
      ELSE
! Initial+FINAL state 
        CALL mdlife(wtfin3,wtfin4)
        CALL mdlifu(wtfin3,wtfin4)
        CALL mdlbhe(wtfin3,wtfin4) ! LL emulation of O(alf2)exp
        CALL mdlbhU(wtfin3,wtfin4) ! Upper line emission, for DEBUG
      ENDIF

      CALL GMONIT( 0,IDA+ 1,WTCRUD*WTSET( 1),1D0,0D0)
      CALL GMONIT( 0,IDA+ 2,WTCRUD*WTSET( 2),1D0,0D0)
      CALL GMONIT( 0,IDA+12,WTCRUD*(WTSET(2)-WTSET(1)),1D0,0D0)
      CALL GMONIT( 0,IDA+ 3,WTCRUD*WTSET( 3),1D0,0D0)
      CALL GMONIT( 0,IDA+13,WTCRUD*(WTSET(3)-WTSET(2)),1D0,0D0)
      CALL GMONIT( 0,IDA+ 4,WTCRUD*WTSET( 4),1D0,0D0)
      CALL GMONIT( 0,IDA+14,WTCRUD*(WTSET(4)-WTSET(3)),1D0,0D0)
      CALL GMONIT( 0,IDA+ 5,WTCRUD*WTSET( 5),1D0,0D0)
      CALL GMONIT( 0,IDA+15,WTCRUD*(WTSET(5)-WTSET(4)),1D0,0D0)

! This is principal weight (third order no pairs)
      WTMOD     =      WTCRUD*WTSET( 4)

!============================================================= 
! for calibration test
      IF(KEYTES.EQ.1)  WTMOD=WTCRUD 
!     ********
      ELSE
!     ********  
! final printout  
      NPAR(10)= NEVGEN
      NPAR(20)= NEVGEN
      CALL VESK2W(1,FUNSKI,AWT,EREL,ZCRUD)
      DWT  = AWT*EREL       
      XCRU = ZCRUD *BORNC
! ....... xsections in concecutive orders......
      CALL GMONIT( 1,IDA+ 1,AWT01,DWT01,DUMM)    
      XS01 = XCRU*AWT01
      DS01 = XCRU*AWT01*DWT01 
      CALL GMONIT( 1,IDA+ 2,AWT02,DWT02,DUMM)    
      XS02 = XCRU*AWT02
      DS02 = XCRU*AWT02*DWT02
      CALL GMONIT( 1,IDA+ 3,AWT03,DWT03,DUMM)    
      XS03 = XCRU*AWT03
      DS03 = XCRU*AWT03*DWT03
      CALL GMONIT( 1,IDA+ 4,AWT04,DWT04,DUMM)    
      XS04 = XCRU*AWT04
      DS04 = XCRU*AWT04*DWT04 
! ....... the differences between orders......
      CALL GMONIT( 1,IDA+12,AWT12,DWT12,DUMM)    
      RXS12 = XCRU*AWT12        /BORN
      RDS12 = XCRU*AWT12*DWT12  /BORN
      CALL GMONIT( 1,IDA+13,AWT13,DWT13,DUMM)    
      RXS13 = XCRU*AWT13        /BORN
      RDS13 = XCRU*AWT13*DWT13  /BORN
      CALL GMONIT( 1,IDA+14,AWT14,DWT14,DUMM)    
      RXS14 = XCRU*AWT14        /BORN
      RDS14 = XCRU*AWT14*DWT14  /BORN
! ... and pairs: third order + pairs, pairs only
      CALL GMONIT( 1,IDA+ 5,AWT05,DWT05,DUMM)    
      XS05 = XCRU*AWT05
      DS05 = XCRU*AWT05*DWT05 
      CALL GMONIT( 1,IDA+15,AWT15,DWT15,DUMM)    
      RXS15 = XCRU*AWT15        /BORN
      RDS15 = XCRU*AWT15*DWT15  /BORN
! .......
      XPAR(10)= XS04
      XPAR(11)= DS04/XS04    
! for WEIGHTED events
      XPAR(20)= XCRU
      XPAR(21)= 0D0     
      XPAR(22)= BORN 
! auxiliary information
      XPAR(25)= DSIG0      
      IF(MODE.EQ.1) RETURN      
!     ====================
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '        OUTPUT FROM              '
      WRITE(NOUT,BXTXT) '  LUMLOG/BHALOG: WINDOW A        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '   X.sect. in [nb] units         '
      WRITE(NOUT,BXTXT) '   for total generated sample    '
      WRITE(NOUT,BXL1I) NEVGEN,     'generated events   ','NEVGEN','A1'
      WRITE(NOUT,BXL2F) AWT ,DWT   ,'vesko2 int. estim. ','AWT   ','A2'
      WRITE(NOUT,BXL1F) XCRU  ,     'crude xsec (vesko2)','XCRU  ','A3'
      WRITE(NOUT,BXL1F) BORN  ,     'Born xsection      ','BORN  ','A4'
      WRITE(NOUT,BXTXT) '        ---  O(alf0)exp ---      '
      WRITE(NOUT,BXL2F) XS01,DS01,  'xsec. total        ','XS01  ','A5'
      WRITE(NOUT,BXL1F) DS01/XS01,  'rel. error         ','      ','  '
      WRITE(NOUT,BXL1F) XS01/BORN-1,'O(alf0)/Born-1     ','      ','  '
      WRITE(NOUT,BXTXT) '        ---  O(alf1)exp ---      '
      WRITE(NOUT,BXL2F) XS02,DS02,  'xsec. total        ','XS02  ','A6'
      WRITE(NOUT,BXL1F) DS02/XS02,  'rel. error         ','      ','  '
      WRITE(NOUT,BXL1F) XS02/BORN-1,'O(alf1)/Born-1     ','      ','  '
      WRITE(NOUT,BXL2F) RXS12,RDS12,'O(alf1-alf0)/Born  ','RXS12 ','  '
      WRITE(NOUT,BXTXT) '        ---  O(alf2)exp ---      '
      WRITE(NOUT,BXL2F) XS03,DS03,  'xsec. total        ','XS03  ','A7'
      WRITE(NOUT,BXL1F) DS03/XS03,  'rel. error         ','      ','  '
      WRITE(NOUT,BXL1F) XS03/BORN-1,'O(alf2)/Born-1     ','      ','  '
      WRITE(NOUT,BXL2F) RXS13,RDS13,'O(alf2-alf1)/Born  ','RXS13 ','  '
      WRITE(NOUT,BXTXT) '        ---  O(alf3)exp ---      '
      WRITE(NOUT,BXL2F) XS04,DS04,  'xsec. total        ','XS04  ','A8'
      WRITE(NOUT,BXL1F) DS04/XS04,  'rel. error         ','      ','  '
      WRITE(NOUT,BXL1F) XS04/BORN-1,'O(alf3)/Born-1     ','      ','  '
      WRITE(NOUT,BXL2F) RXS14,RDS14,'O(alf3-alf2)/Born  ','RXS14 ','  '
      WRITE(NOUT,BXTXT) '        ---  O(alf3)exp ---      '
      WRITE(NOUT,BXTXT) '        ---  plus pairs ---      '
      WRITE(NOUT,BXL2F) XS05,DS05,  'xsec. total        ','XS05  ','A9'
      WRITE(NOUT,BXL1F) DS05/XS05,  'rel. error         ','      ','  '
      WRITE(NOUT,BXL1F) XS05/BORN-1,'O(alf3+prs)/Born-1 ','      ','  '
      WRITE(NOUT,BXL2F) RXS15,RDS15,'pairs/Born         ','RXS15 ','  '
      WRITE(NOUT,BXCLO)  
      ENDIF
!     ********
      END       


      SUBROUTINE KINOLT        
!     *****************        
! construction of four-momenta 
! They are stored in /MOMSET/ which replaces /UTILUS/
!     ***********************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      COMMON / INOUT  / NINP,NOUT 
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI  
! final state fermions in BHLUMI output format
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      COMMON / TRANSR / TRAN,TRANX,TRANY
!
      DOUBLE PRECISION DRVEC(100)
      DATA ICONT /0/
!
      ICONT = ICONT+1 
      ENE  = CMSENE/2D0 
!...beams
      P1(4)=  ENE
      P1(3)=  ENE
      P1(2)=  0D0
      P1(1)=  0D0
      Q1(4)=  ENE
      Q1(3)= -ENE
      Q1(2)=  0D0
      Q1(1)=  0D0 
!...Explicit collinear photons 
!   in the initial state (parallel to beams) and in the final state
      NPHOT= 2
!...outgoing dressed (calorimetric) electrons  
      Z1   = 1D0-X1
      Z2   = 1D0-X2
      z3   = 1 -x3
      z4   = 1 -x4
      TRAN = CMSENE**2*Z1*Z2*XI
      CDOT = 1D0-2D0*XI 
      SDOT = DSQRT(DABS(4D0*XI*(1D0-XI))) 
      XPT  = 2D0*DSQRT(Z1*Z2)*SDOT
      CALL VARRAN(DRVEC,1)
      PHI  = 2D0*PI*DRVEC(1)
!...first electron
      P2(4) = 0.5D0*ENE*Z3 *( Z1+Z2   +CDOT*(Z1-Z2))
      P2(3) = 0.5D0*ENE*Z3 *( CDOT*(Z1+Z2) +Z1-Z2) 
      P2(2) = 0.5D0*ENE*Z3 * XPT*DSIN(PHI) 
      P2(1) = 0.5D0*ENE*Z3 * XPT*DCOS(PHI)
!...second electron
      Q2(4) = 0.5D0*ENE*Z4 *( Z1+Z2   -CDOT*(Z1-Z2))
      Q2(3) = 0.5D0*ENE*Z4 *( -CDOT*(Z1+Z2) +Z1-Z2) 
      Q2(2) =-0.5D0*ENE*Z4 * XPT*DSIN(PHI) 
      Q2(1) =-0.5D0*ENE*Z4 * XPT*DCOS(PHI)
!...Photon(s) parallel to first beam
      PHOT(1,4)=  ENE*(1-Z1)
      PHOT(1,3)=  ENE*(1-Z1)
      PHOT(1,2)=  0
      PHOT(1,1)=  0
!...Photon(s) parallel to second beam
      PHOT(2,4)=  ENE*(1-Z2)
      PHOT(2,3)= -ENE*(1-Z2)
      PHOT(2,2)=  0
      PHOT(2,1)=  0
!...Photons parallel to outgoing electron/positron
      keyfin = mod(keyrad,1000)/100    
      IF(keyfin.eq.1) THEN
        nphot= 4
        DO k=1,4
           phot(3,k)= (1-z3)/z3 *p2(k)
           phot(4,k)= (1-z4)/z4 *q2(k)
        ENDDO
      ENDIF
      END

      FUNCTION FUNSKI(T1,T2)
!     ************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
! FLEPS is the smallest floating on a given instalation  
      PARAMETER(X0=1D-10)
!c===      PARAMETER(FLEPS= 1D-300)
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES

      TT1=T1
      TT2=T2
      XMAX = 1D0-ZMIN  
!c===      SMALL = FLEPS**ETA
!c===      IF(SMALL.GT.1D-4) 
!c===     $   WRITE(6,*) ' ++++ FUNSKO: warning FLEPS**ETA=',SMALL
!c===      IF(T1.LT.SMALL.OR.T2.LT.SMALL) GOTO 900 
      XX1=0D0
      XX2=0D0
      IF(T1.GT.X0**ETA) XX1 = XMAX*DEXP(1/ETA*DLOG(T1))
      IF(T2.GT.X0**ETA) XX2 = XMAX*DEXP(1/ETA*DLOG(T2))  
! Jacobian factor due to change of variables t=>x
! eta*x**(eta-1) multiplied  (numerical stability)  see DD1*DD2
      RJAC = XMAX**(2D0*ETA)
!=== $     /(ETA*XX1**(ETA-1D0))/(ETA*XX2**(ETA-1D0))
! anticipated angular trigger
      Z1 = 1-XX1
      Z2 = 1-XX2   
      SLOPE = Z1/Z2 
      SRAXI = SQRT(XIB/XIA)
      IF(SLOPE.GT.  SRAXI) GOTO 900
      IF(SLOPE.LT.1/SRAXI) GOTO 900  
! Ordinary flux-factor      
      FLUX = 1D0/(Z1*Z2)
! Jacobian due to xi=>xi1 or xi=>xi2 change
      IF( Z1.LT.Z2) THEN
        DJAC= Z1/Z2
      ELSE
        DJAC= Z2/Z1
      ENDIF  
! Crude structure functions
! eta*x**(eta-1) divided out (numerical stability)  see RJAC
!===  DD1 = ETA*XX1**(ETA-1D0)*(1D0+(1D0-XX1)**2)/2D0
!===  DD2 = ETA*XX2**(ETA-1D0)*(1D0+(1D0-XX2)**2)/2D0
      DD1 =                    (1D0+(1D0-XX1)**2)/2D0
      DD2 =                    (1D0+(1D0-XX2)**2)/2D0
! Test distributions
      IF(KEYTES.EQ.1) THEN
!===    DD1 = ETA*XX1**(ETA-1D0)
!===    DD2 = ETA*XX2**(ETA-1D0)
        DD1 = 1D0       
        DD2 = 1D0   
      ENDIF
! below no trace of eta*x**(eta-1) any more!
      FUNSKI = FLUX*RJAC*DJAC*DD1*DD2
      RETURN
  900 FUNSKI = 0D0
      END
  

      FUNCTION BORNB(CMSENE,THA,THB)
!     *******************************
! BORN XSECTION pure t-channel     
! THA,THB are in radians, CMSENE in GEV units
! result in nanobarns
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      PARAMETER( PI = 3.1415926535897932D0, ALFINV = 137.03604D0)
      PARAMETER( ALFPI=  1D0/PI/ALFINV ,ALFA=1D0/ALFINV)
      PARAMETER( GNANOB=389.385D-30*1.D33 )
      EXTERNAL BORXI
      XIA= (1D0-DCOS(THA))/2D0
      XIB= (1D0-DCOS(THB))/2D0
      DSIG0 = 4D0*ALFA**2*PI/CMSENE**2*GNANOB    
      CALL GAUSJD(BORXI,XIA,XIB,-1D-6,RESULT)
      BORNB= RESULT *DSIG0
      END
      FUNCTION BORXI(XI)
!     ******************
! Integrand of BORNB
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BORXI=1D0/XI**2*(1D0+(1D0-XI)**2)/2D0
      END

      FUNCTION FOBOR(XIA,XIB)
!     ***********************
! BORN XSECTION for pure t-channel     
! XI=(1-cos(theta))/2
! result in DSIG0 units where
! DSIG0 = 4D0*ALFA**2*PI/CMSENE**2
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      EXTERNAL BORNXF
      CALL GAUSJD(BORNXF,XIA,XIB,-1D-6,RESULT)
      FOBOR= RESULT
      END
      FUNCTION BORNXF(XI)
!     ******************
! Integrand for FOBOR
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BORNXF=1D0/XI**2*(1D0+(1D0-XI)**2)/2D0
!----------------------------------------------------------------C
!                  The end of LUMLOG                             C
!----------------------------------------------------------------C
      END

      SUBROUTINE Brelos(eta,tau,x,wt)
!     *******************************
! Generates single x=1-z according to crude distribution
! D(gam,x) =gam*x**(eta-1)*chi(x)
! Note: wt = <wtint> = \int_0^1 D(gam,x) dx
!       Up to first order wt=exp(-3/4*eta)=1-3/4*eta
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER(x0=1d-5)
      DOUBLE PRECISION drvec(10)
! Inline function
      chi(x)= (1+(1-x)**2)/2
!
 100  CONTINUE
      CALL varran(drvec,2)
      tau = drvec(1)
      x = 0
      IF(tau.GT.x0**eta) x = exp( 1/eta*log(tau) )
      rn =  drvec(2)
      wtint =chi(x)
      IF(rn.gt.wtint) goto 100
      wt = 1 -eta/(1+eta) +eta/2/(2+eta)
      END



      SUBROUTINE mdline
!     *****************
! This is set of weights for the exponentiated case up to O(gam**3)exp.
! Explicit emission INITIAL STATE ONLY (final electron dressed)
!---------------------------------------------------------------------
!===  In order to assure numerical stability
!===  factors x**(beta-1) =[xmax*t**(1/eta)]**(beta-1)
!===  are divided off by hand, here and in in STRUFU.
!===  (Similar division of x**(eta-1) in FUNSKI is independent.)
!---------------------------------------------------------------------
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(ALFA   =  1D0/137.03604D0,   ALF1=  ALFA/PI)
      PARAMETER(AMEL   =  0.511D-3)       
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      COMMON / INOUT  / NINP,NOUT 
!
      KEYBLO = MOD(KEYRAD,10) 
      XMAX=  1-ZMIN
      xx1 =  x1
      xx2 =  x2
      z1  =1-x1
      z2  =1-x2
!---------------------------------------------------------------------
! Two types of the big-log definition
      IF(KEYBLO.EQ.3) THEN 
        SPRIM = CMSENE**2*Z1*Z2
        BETA  = ALF1*(DLOG(SPRIM*XI/AMEL**2) -1)
      ELSEIF(KEYBLO.EQ.4) THEN
        BETA  = ALF1* DLOG((CMSENE/AMEL)**2*XIA) 
![[[[[        BETA  = ALF1* DLOG((CMSENE/AMEL)**2*XIB) 
      ELSE
        WRITE(6,*) ' WRONG KEYBLO=', KEYBLO
      ENDIF
!---------------------------------------------------------------------
! Born Factor
      WTBOR = (1+(1-XI)**2)/2D0      
!---------------------------------------------------------------------
!===  DDCRU =  ETA*XX1**(ETA-1D0)*(1D0+(1D0-XX1)**2)/2D0
!=== $        *ETA*XX2**(ETA-1D0)*(1D0+(1D0-XX2)**2)/2D0 
      DDCRU =  ETA*               (1D0+(1D0-XX1)**2)/2D0
     $        *ETA*               (1D0+(1D0-XX2)**2)/2D0 
     $        *XMAX**( ETA-BETA)*TT1**(1-BETA/ETA)
     $        *XMAX**( ETA-BETA)*TT2**(1-BETA/ETA)
! zero order
      DDMZR     =  STRUFU(310,BETA,XX1)*STRUFU(310,BETA,XX2)
      WTSET( 1) =  DDMZR/DDCRU*WTBOR 
! first order
      DDMO1     =  STRUFU(311,BETA,XX1)*STRUFU(311,BETA,XX2)
      WTSET( 2) =  DDMO1/DDCRU*WTBOR 
! second order
      DDMO2     =  STRUFU(312,BETA,XX1)*STRUFU(312,BETA,XX2)
      WTSET( 3) =  DDMO2/DDCRU*WTBOR 
! third order
      DDMO3     =  STRUFU(313,BETA,XX1)*STRUFU(313,BETA,XX2)
      WTSET( 4) =  DDMO3/DDCRU*WTBOR 
! third order + pairs       
      BETR = -3D0*DLOG(1-BETA/3D0)
      DDCRR =  ETA*               (1D0+(1D0-XX1)**2)/2D0
     $        *ETA*               (1D0+(1D0-XX2)**2)/2D0 
     $        *XMAX**( ETA-BETR)*TT1**(1-BETR/ETA)
     $        *XMAX**( ETA-BETR)*TT2**(1-BETR/ETA)
      DDMO4     =  STRUFU(313,BETR,XX1)*STRUFU(313,BETR,XX2)
      WTSET( 5) =  DDMO4/DDCRR*WTBOR 
      END

      SUBROUTINE mdlife(wtfin3,wtfin4)
!     ********************************
! This is set of weights for the exponentiated case up to O(gam**3)exp.
! Explicit emission INITIAL+FINAL STATE
!---------------------------------------------------------------------
!===  In order to assure numerical stability
!===  factors x**(beta-1) =[xmax*t**(1/eta)]**(beta-1)
!===  are divided off by hand, here and in in STRUFU.
!===  (Similar division of x**(eta-1) in FUNSKI is independent.)
!---------------------------------------------------------------------
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(ALFA   =  1D0/137.03604D0,   ALF1=  ALFA/PI)
      PARAMETER(AMEL   =  0.511D-3)       
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      COMMON / INOUT  / NINP,NOUT 
! Inline function
      chi(x)= (1+(1-x)**2)/2
!
      KEYBLO = MOD(KEYRAD,10) 
      XMAX=  1-ZMIN
      z1  =1-x1
      z2  =1-x2
!---------------------------------------------------------------------
! Two types of the big-log definition
      IF(KEYBLO.EQ.3) THEN 
        sprim = cmsene**2*z1*z2
        beta  = alf1*(dlog(sprim*xi/amel**2) -1)
      ELSEIF(KEYBLO.EQ.4) THEN
        beta  = alf1* dlog((cmsene/amel)**2*xia) 
![[[[[        beta  = alf1* dlog((cmsene/amel)**2*xib) 
      ELSE
        WRITE(6,*) ' WRONG KEYBLO=', KEYBLO
      ENDIF
!---------------------------------------------------------------------
! Born Factor
      WTBOR = chi(xi)
      WTfin = wtfin3*wtfin4
!---------------------------------------------------------------------
!===  DDCRU =  eta*x1**(eta-1)*chi(x1)*eta*x2**(eta-1)*chi(x2)
!=== $        *eta*x3**(eta-1)*chi(x3)*eta*x4**(eta-1)*chi(x4)
      DDCRU =  eta*            chi(x1)*eta*            chi(x2)
     $        *eta*            chi(x3)*eta*            chi(x4)
     $        *xmax**( eta-beta)*tt1**(1-beta/eta)
     $        *xmax**( eta-beta)*tt2**(1-beta/eta)
     $                          *tt3**(1-beta/eta)
     $                          *tt4**(1-beta/eta)
! zero order
      DDMZR     =  strufu(310,beta,x1)*strufu(310,beta,x2)
     $            *strufu(310,beta,x3)*strufu(310,beta,x4)
      WTSET( 1) =  DDMZR/DDCRU*WTBOR *WTfin
! first order
      DDMO1     =  strufu(311,beta,x1)*strufu(311,beta,x2)
     $            *strufu(311,beta,x3)*strufu(311,beta,x4)
      WTSET( 2) =  DDMO1/DDCRU*WTBOR *WTfin
! second order
      DDMO2     =  strufu(312,beta,x1)*strufu(312,beta,x2)
     $            *strufu(312,beta,x3)*strufu(312,beta,x4)
      WTSET( 3) =  DDMO2/DDCRU*WTBOR *WTfin
! third order
      DDMO3     =  strufu(313,beta,x1)*strufu(313,beta,x2)
     $            *strufu(313,beta,x3)*strufu(313,beta,x4)
      WTSET( 4) =  DDMO3/DDCRU*WTBOR *WTfin
! third order + pairs       
      BETR = -3D0*DLOG(1-BETA/3D0)
      DDCRR =  eta*           chi(x1)*eta*           chi(x2)
     $        *eta*           chi(x3)*eta*           chi(x4)
     $        *xmax**( eta-betr)*tt1**(1-betr/eta)
     $        *xmax**( eta-betr)*tt2**(1-betr/eta)
     $                          *tt3**(1-betr/eta)
     $                          *tt4**(1-betr/eta)
      DDMO4     =  strufu(313,betr,x1)*strufu(313,betr,x2)
     $            *strufu(313,betr,x3)*strufu(313,betr,x4)
      WTSET( 5) =  DDMO4/DDCRR*WTBOR *WTfin
!---------------------------------------------------------------------
!          Upper line only
!---------------------------------------------------------------------
!===  DDCRU =  eta*x1**(eta-1)*chi(x1) *eta*x2**(eta-1)*chi(x2)
!=== $        *eta*x3**(eta-1)*chi(x3) *eta*x4**(eta-1)*chi(x4)
      DDCRU =  eta*            chi(x1) *eta*            chi(x2)
     $        *eta*            chi(x3) *eta*            chi(x4)
     $        *xmax**( eta-beta)*tt1**(1-beta/eta)
     $        *xmax**( eta-beta)*tt2**(1-beta/eta)
     $                          *tt3**(1-beta/eta)
     $                          *tt4**(1-beta/eta)
! third order
!===  DDMO6     =  strufu(313,beta,x1) *beta*x2**(beta-1)
!=== $            *strufu(313,beta,x3) *beta*x4**(beta-1)
      DDMO6     =  strufu(313,beta,x1) *beta
     $            *strufu(313,beta,x3) *beta
      WTSET( 6) =  DDMO6/DDCRU*WTBOR *WTfin
!---------------------------------------------------------------------
      END

      SUBROUTINE mdlinu      
!     ***********************
!  Model weight for second order LL unexponentiated
!  note infrared cut k0 is present as usual
!  (stj) 14 febr 91, 
!  correction of bug xk0**eta instead of (xk0/xmax)**eta
!
! correction 5.05.92     E.Richter-Was
! soft region x1<k0, x2<k0  is not corrected to true crude
! str. function i.e. should be divided by 1/((1+(1-x1)**2)/2) *
! 1/((1+(1-x2)**2)/2)   --- effect completelly negligible numerically
! futher corrections included in fortran code
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(ALFA   =  1D0/137.03604D0,   ALF1=  ALFA/PI)
      PARAMETER(GNANOB =  389385D0) 
      PARAMETER(AMEL   =  0.511D-3)       
      PARAMETER(X0     =  1D-5)       
      COMMON / INOUT  / NINP,NOUT 
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)

      KEYBLO = MOD(KEYRAD,10) 

      XMAX=  1-ZMIN   
      T0 = (XK0/XMAX)**ETA
      XX1=0D0
      IF(TT1.GT.T0) XX1 = XMAX*DEXP(1/ETA*DLOG(TT1))
      XX2=0D0
      IF(TT2.GT.T0) XX2 = XMAX*DEXP(1/ETA*DLOG(TT2))  
      Z1 = 1D0-XX1
      Z2 = 1D0-XX2 
! Born Factor
      WTBOR = (1+(1-XI)**2)/2D0    
! Two types of the big-log definition
      IF(KEYBLO.EQ.3) THEN 
        SPRIM = CMSENE**2*Z1*Z2
        BETA  = ALF1*(DLOG(SPRIM*XI/AMEL**2) -1)
      ELSEIF(KEYBLO.EQ.4) THEN
        BETA  = ALF1* DLOG((CMSENE/AMEL)**2*XIA) 
![[[[[        BETA  = ALF1* DLOG((CMSENE/AMEL)**2*XIB) 
      ELSE
        WRITE(6,*) ' WRONG KEYBLO=', KEYBLO
      ENDIF  
!   WTSET(11)   Born
!   WTSET(12)   O(alf1)
! 
      WTSET(11) =0D0
      WTSET(12) =0D0 
      WTSET(13) =0D0 
! soft region x1<k0, x2<k0
!CCCCCCCC      IF(TT1.LT.XK0**ETA.AND.TT2.LT.XK0**ETA) THEN
      IF(TT1.LT.T0.AND.TT2.LT.T0) THEN
         WTSET(11) = 1D0/XK0**(2*ETA) * WTBOR
         WTSET(12) = (1D0+ 2D0*BETA*DLOG(XK0)+ 3D0/2D0*BETA)
     $        /XK0**(2*ETA) * WTBOR
         WTSET(13) = 
     $  ( 1D0+4D0*BETA/4D0*(3/2D0+2D0*DLOG(XK0))
     $    +    (BETA/2D0)**2 *(3/2D0+2D0*DLOG(XK0))**2
     $    +2D0*(BETA/2D0)**2 *(9/8D0+2D0*(DLOG(XK0))**2
     $           +3D0*DLOG(XK0)    -PI**2/3D0))      
     $        /XK0**(2*ETA) * WTBOR
      ENDIF       
! single bremss. one soft photon (virt.) one hard  
!CCC      IF(TT1.LT.XK0**ETA.AND.TT2.GT.XK0**ETA) THEN
      IF(TT1.LT.T0.AND.TT2.GT.T0) THEN
         WTSET(12) = 
     $     BETA/( XK0**ETA *ETA *XX2**ETA)  * WTBOR  
         WTSET(13) =
     $      BETA *(1D0+BETA/2D0*( 3/2D0+2D0*DLOG(XK0)
     $                 +2D0*DLOG(XX2)-DLOG(1D0-XX2)+3/2D0)
!c5.05$    +BETA/2D0*((2D0-XX2)/2D0*DLOG(1D0-XX2)-XX2)*XX2 )
     $     +BETA/2D0*((2D0-XX2)/2D0*DLOG(1D0-XX2)-XX2)*XX2
     $                                      /(1D0+(1D0-XX2)**2) )
     $    /( XK0**ETA *ETA *XX2**ETA)  * WTBOR
      ENDIF                 
!CCC      IF(TT1.GT.XK0**ETA.AND.TT2.LT.XK0**ETA) THEN
      IF(TT1.GT.T0.AND.TT2.LT.T0) THEN
         WTSET(12) = 
     $     BETA/( XK0**ETA *ETA *XX1**ETA)  * WTBOR
         WTSET(13) =
     $      BETA *(1D0+BETA/2D0*( 3/2D0+2D0*DLOG(XK0)
     $                 +2D0*DLOG(XX1)-DLOG(1D0-XX1)+3/2D0)
!C5.05$     +BETA/2D0*((2D0-XX1)/2D0*DLOG(1D0-XX1)-XX1)*XX1 )
     $     +BETA/2D0*((2D0-XX1)/2D0*DLOG(1D0-XX1)-XX1)*XX1 
     $                                      /(1D0+(1D0-XX1)**2) )
     $    /( XK0**ETA *ETA *XX1**ETA)  * WTBOR
      ENDIF             
! two hardreal  photons above k0    
!CCCC      IF(TT1.GT.XK0**ETA.AND.TT2.GT.XK0**ETA) THEN
      IF(TT1.GT.T0.AND.TT2.GT.T0) THEN
         WTSET(13) = 
     $     BETA**2/( ETA**2 *XX1**ETA*XX2**ETA)  * WTBOR
      ENDIF
      END

      SUBROUTINE mdlifu(wtfin3,wtfin4) 
!     ********************************
!  Model weight for second order LL un-exponentiated INITIAL+FINAL
!  note infrared cut k0 is present as usual
!  WTSET(11)   Born
!  WTSET(12)   O(alf1) without exponentiation
!  WTSET(13)   O(alf2) without exponentiation
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(ALFA   =  1D0/137.03604D0,   ALF1=  ALFA/PI)
      PARAMETER(GNANOB =  389385D0) 
      PARAMETER(AMEL   =  0.511D-3)       
      PARAMETER(X0     =  1D-5)       
      COMMON / INOUT  / NINP,NOUT 
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      LOGICAL sof1,sof2,sof3,sof4,har1,har2,har3,har4
! Inline function
      chi(x)= (1+(1-x)**2)/2

      KEYBLO = MOD(KEYRAD,10) 
! Here xk0 is infrared dummy parameter,
      xmax=  1-zmin   
      t0I = (xk0/xmax)**eta
      t0F =  xk0**eta
      xx1=0
      xx2=0
      xx3=0
      xx4=0
      IF(tt1.gt.t0I) xx1 = xmax*exp(1/eta*log(tt1))
      IF(tt2.gt.t0I) xx2 = xmax*exp(1/eta*log(tt2))  
      IF(tt3.gt.t0F) xx3 =      exp(1/eta*log(tt3))
      IF(tt4.gt.t0F) xx4 =      exp(1/eta*log(tt4))  
      z1 = 1-xx1
      z2 = 1-xx2
      z3 = 1-xx3
      z4 = 1-xx4
! Born Factor
      WTBOR = chi(xi) *wtfin3*wtfin4   
! Two types of the big-log definition
      IF(KEYBLO.EQ.3) THEN 
        sprim = cmsene**2*z1*z2
        beta  = alf1*(dlog(sprim*xi/amel**2) -1)
      ELSEIF(KEYBLO.EQ.4) THEN
        beta  = alf1* dlog((cmsene/amel)**2*xia) 
![[[[[        beta  = alf1* dlog((cmsene/amel)**2*xib) 
      ELSE
        WRITE(6,*) ' WRONG KEYBLO=', KEYBLO
      ENDIF  
      WTSET(11) =0
      WTSET(12) =0 
      WTSET(13) =0 
      sof1 = tt1.lt.t0I
      sof2 = tt2.lt.t0I
      sof3 = tt3.lt.t0F
      sof4 = tt4.lt.t0F
      har1 = tt1.ge.t0I
      har2 = tt2.ge.t0I
      har3 = tt3.ge.t0F
      har4 = tt4.ge.t0F
! all four soft
      IF(sof1.AND.sof2.AND.sof3.AND.sof4) THEN
        WTSET(11) =  (1/xk0**eta)**4 * wtbor
        WTSET(12) = (1 + 4*beta*(3./4+log(xk0))  )
     $              *(1/xk0**eta)**4 * wtbor 
        WTSET(13) = 
     $  (1 +4*beta*(3./4+log(xk0)) -4*beta**2*(1./12*pi**2)
     $     +1./2*(4*beta)**2*(3./4+log(xk0))**2 )
     $              *(1/xk0**eta)**4 * wtbor 
      ENDIF
!-----------------------------------------------------
! one hard other three soft
      IF(har1.AND.sof2.AND.sof3.AND.sof4) THEN
        WTSET(12) =  beta
     $    *1/(eta*xx1**eta) *(1/xk0**eta)**3   * wtbor
        WTSET(13) = (beta *chi(xx1) -beta**2 *(1./4)*xx1**2
     $    -beta**2 *log(1-xx1) *(1./8)*(4-6*xx1+3*xx1**2)
     $    +beta**2 *chi(xx1)   *(3./4+log(xx1)) 
     $  +3*beta**2 *chi(xx1)   *(3./4+log(xk0)) 
     $  )  *1/(eta*xx1**eta*chi(xx1)) *(1/xk0**eta)**3 * wtbor
      ENDIF
      IF(sof1.AND.har2.AND.sof3.AND.sof4) THEN
        WTSET(12) =  beta
     $    *1/(eta*xx2**eta) *(1/xk0**eta)**3   * wtbor
        WTSET(13) = (beta *chi(xx2) -beta**2 *(1./4)*xx2**2
     $    -beta**2 *log(1-xx2) *(1./8)*(4-6*xx2+3*xx2**2)
     $    +beta**2 *chi(xx2)   *(3./4+log(xx2)) 
     $  +3*beta**2 *chi(xx2)   *(3./4+log(xk0)) 
     $  )  *1/(eta*xx2**eta*chi(xx2)) *(1/xk0**eta)**3 * wtbor
      ENDIF
      IF(sof1.AND.sof2.AND.har3.AND.sof4) THEN
        WTSET(12) =  beta
     $    *1/(eta*xx3**eta) *(1/xk0**eta)**3   * wtbor
        WTSET(13) = (beta *chi(xx3) -beta**2 *(1./4)*xx3**2
     $    -beta**2 *log(1-xx3) *(1./8)*(4-6*xx3+3*xx3**2)
     $    +beta**2 *chi(xx3)   *(3./4+log(xx3)) 
     $  +3*beta**2 *chi(xx3)   *(3./4+log(xk0)) 
     $  )  *1/(eta*xx3**eta*chi(xx3)) *(1/xk0**eta)**3 * wtbor
      ENDIF
      IF(sof1.AND.sof2.AND.sof3.AND.har4) THEN
        WTSET(12) =  beta
     $    *1/(eta*xx4**eta) *(1/xk0**eta)**3   * wtbor
        WTSET(13) = (beta *chi(xx4) -beta**2 *(1./4)*xx4**2
     $    -beta**2 *log(1-xx4) *(1./8)*(4-6*xx4+3*xx4**2)
     $    +beta**2 *chi(xx4)   *(3./4+log(xx4)) 
     $  +3*beta**2 *chi(xx4)   *(3./4+log(xk0)) 
     $  )  *1/(eta*xx4**eta*chi(xx4)) *(1/xk0**eta)**3 * wtbor
      ENDIF
!-----------------------------------------------------
! two hard real  photons above k0, other two are soft  
      IF(har1.AND.har2.AND.sof3.AND.sof4) THEN
        WTSET(13) = beta**2 
     $   /(eta*xx1**eta) /(eta*xx2**eta) *(1/xk0**eta)**2 * wtbor
      ENDIF
      IF(har1.AND.sof2.AND.har3.AND.sof4) THEN
        WTSET(13) = beta**2 
     $   /(eta*xx1**eta) /(eta*xx3**eta) *(1/xk0**eta)**2 * wtbor
      ENDIF
      IF(har1.AND.sof2.AND.sof3.AND.har4) THEN
        WTSET(13) = beta**2 
     $   /(eta*xx1**eta) /(eta*xx4**eta) *(1/xk0**eta)**2 * wtbor
      ENDIF
      IF(sof1.AND.har2.AND.har3.AND.sof4) THEN
        WTSET(13) = beta**2 
     $   /(eta*xx2**eta) /(eta*xx3**eta) *(1/xk0**eta)**2 * wtbor
      ENDIF
      IF(sof1.AND.har2.AND.sof3.AND.har4) THEN
        WTSET(13) = beta**2 
     $   /(eta*xx2**eta) /(eta*xx4**eta) *(1/xk0**eta)**2 * wtbor
      ENDIF
      IF(sof1.AND.sof2.AND.har3.AND.har4) THEN
        WTSET(13) = beta**2 
     $   /(eta*xx3**eta) /(eta*xx4**eta) *(1/xk0**eta)**2 * wtbor
      ENDIF
      END        



      SUBROUTINE modlu 
!     *****************
!  This routine is not used, it is completely equivalent to mdlinu
!  but has nicer representation of all formulas
!-----------------------------------------------------------------
!  Model weight for second order LL un-exponentiated INITIAL+FINAL
!  note infrared cut k0 is present as usual
!  WTSET(11)   Born
!  WTSET(12)   O(alf1)
!  WTSET(13)   O(alf2)
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(ALFA   =  1D0/137.03604D0,   ALF1=  ALFA/PI)
      PARAMETER(GNANOB =  389385D0) 
      PARAMETER(AMEL   =  0.511D-3)       
      PARAMETER(X0     =  1D-5)       
      COMMON / INOUT  / NINP,NOUT 
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      LOGICAL sof1,sof2,har1,har2
! Inline function
      chi(x)= (1+(1-x)**2)/2

      KEYBLO = MOD(KEYRAD,10) 

      XMAX=  1-ZMIN   
      T0 = (XK0/XMAX)**ETA
      xx1=0
      IF(tt1.gt.t0) xx1 = xmax*dexp(1/eta*dlog(tt1))
      xx2=0
      IF(tt2.gt.t0) xx2 = xmax*dexp(1/eta*dlog(tt2))  
      z1 = 1-xx1
      z2 = 1-xx2 
! Born Factor
      WTBOR = chi(xi)    
! Two types of the big-log definition
      IF(KEYBLO.EQ.3) THEN 
        sprim = cmsene**2*z1*z2
        beta  = alf1*(dlog(sprim*xi/amel**2) -1)
      ELSEIF(KEYBLO.EQ.4) THEN
        beta  = alf1* dlog((cmsene/amel)**2*xia) 
![[[[[        beta  = alf1* dlog((cmsene/amel)**2*xib) 
      ELSE
        WRITE(6,*) ' WRONG KEYBLO=', KEYBLO
      ENDIF  
      WTSET(11) =0
      WTSET(12) =0 
      WTSET(13) =0 
      sof1 = tt1.lt.t0
      sof2 = tt2.lt.t0
      har1 = tt1.ge.t0
      har2 = tt2.ge.t0
! soft region x1<k0, x2<k0
      IF(sof1.AND.sof2) THEN
        WTSET(11) =  1/xk0**(2*eta) * wtbor
        WTSET(12) = (1 + 2*beta*(3./4+log(xk0))  )
     $                /xk0**(2*eta) * wtbor
        WTSET(13) = 
     $  (1 +2*beta*(3./4+log(xk0)) -2*beta**2*(1./12*pi**2)
     $     +1./2*(2*beta)**2*(3./4+log(xk0))**2 )
     $  /xk0**eta /xk0**eta * wtbor 
      ENDIF       
! single bremss. one soft photon (virt.) one hard  
      IF(sof1.AND.har2) THEN
        WTSET(12) =  beta /xk0**eta /(eta*xx2**eta)  * wtbor  
        WTSET(13) =
     $  (   beta *chi(xx2)
     $  +beta**2 *chi(xx2)   *(3./4+log(xk0) +3./4+log(xx2)) 
     $  -beta**2 *log(1-xx2) *(1./8)*(4-6*xx2+3*xx2**2)
     $  -beta**2 *(1./4)*xx2**2
     $  ) /xk0**eta  /(eta*xx2**eta*chi(xx2))  * wtbor
      ENDIF                 
      IF(har1.AND.sof2) THEN
        WTSET(12) =  beta /xk0**eta /(eta*xx1**eta)  * wtbor
        WTSET(13) =
     $  (   beta *chi(xx1)
     $  +beta**2 *chi(xx1)   *(3./4+log(xk0) +3./4+log(xx1)) 
     $  -beta**2 *log(1-xx1) *(1./8)*(4-6*xx1+3*xx1**2)
     $  -beta**2 *(1./4)*xx1**2
     $  ) /xk0**eta  /(eta*xx1**eta*chi(xx1))  * wtbor
      ENDIF             
! two hardreal  photons above k0    
      IF(har1.AND.har2) THEN
        WTSET(13) = beta**2 
     $     /(eta*xx1**eta) /(eta*xx2**eta)  * wtbor
      ENDIF
      END        

      FUNCTION STRUFU(KEYD,BETA,VV)
!-------------------------------------------------------------------C
! This originates from XSTFIG 
! non-singlet structure functions, factor x**(beta-1) removed!!!
!-------------------------------------------------------------------C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter( nout=6 )
      PARAMETER (DZ2=1.6449340668482264D0,DZ3=1.2020569031595943D0)
      PARAMETER(PI= 3.1415926535897932D0, ALFINV=137.03604D0)
      PARAMETER(CEULER =0.57721566D0)
   
      BETI=BETA
      X=VV
      Z=1D0-VV
! ....Zero Order without subleading terms
      IF(KEYD  .EQ.310)  THEN
!===     DISTR=BETI*X**(BETI-1D0)
         DISTR=BETI
     @      *EXP(BETI*(.75D0-CEULER))/DPGAMM(1D0+BETI)
! ....First Order without subleading terms
      ELSEIF(KEYD  .EQ.311)  THEN
!===     DISTR=BETI*X**(BETI-1D0)
         DISTR=BETI
     @      *EXP(BETI*(.75D0-CEULER))/DPGAMM(1D0+BETI)
     @      *(1D0 -.5D0*(1D0-(1D0-X)**2) )
! ....Second Order without subleading terms
      ELSEIF(KEYD  .EQ.312)  THEN
!===     DISTR=BETI*X**(BETI-1D0)
         DISTR=BETI
     @      *EXP(BETI*(.75D0-CEULER))/DPGAMM(1D0+BETI)
     @      *(1D0 -.5D0*(1D0-(1D0-X)**2)
     @       +BETI*( -(1D0+3D0*(1D0-X)**2)/8D0*LOG(1D0-X) -.25D0*X**2)
     @       )
! ....Third Order without subleading terms
      ELSEIF(KEYD  .EQ.313)  THEN
! redefine beta in the case of fermion pairs, 
! (this is: LL nonsinglet, electron pair only)
!===     DISTR=BETI*X**(BETI-1D0)        
         DISTR=BETI
     @      *EXP(BETI*(.75D0-CEULER))/DPGAMM(1D0+BETI)
     @      *(1D0 -.5D0*(1D0-(1D0-X)**2)
     @       +BETI*( -(1D0+3D0*(1D0-X)**2)/8D0*LOG(1D0-X) -.25D0*X**2)
     @          +BETI**2*(
     @                +(3D0*X-2D0)*X/16*LOG(1D0-X)
     @                 +(8D0-14D0*X+7D0*X**2)/96*LOG(1D0-X)**2
     @                 +X**2/8D0
     @                  +(2D0-X)*X/8*DILOGY(X)
     @                    )  )
      ELSE
            GOTO 900
      ENDIF        
      STRUFU=DISTR      
      RETURN
 900  WRITE(NOUT,*) ' ===--->  WRONG KEYDIS IN STRUFU'
      STOP
      END 




      SUBROUTINE mdlbhU(wtfin3,wtfin4)
!     ********************************
! LL emulation of BHLUMI amplitude up to O(gam**3)exp.
! Explicit emission UPPER LINE INITIAL+FINAL 
! ---------------------------------------------------------------------
! ===  In order to assure numerical stability
! ===  factors x**(beta-1) =[xmax*t**(1/eta)]**(beta-1) removed.
! ===  Similar division of x**(eta-1) is done in FUNSKI.
! ---------------------------------------------------------------------
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(ALFA   =  1D0/137.03604D0,   ALF1=  ALFA/PI)
      PARAMETER(AMEL   =  0.511D-3)       
      PARAMETER(CEULER =0.57721566D0)
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      COMMON / INOUT  / NINP,NOUT 
! ---------------------------------------------------------------------
! Inline functions
! ---------------------------------------------------------------------
      chi(x)   = (1+(1-x)**2)/2
      FF(g)    = EXP( -g*ceuler )/dpgamm( 1d0 + g )
      STR(g,v) = g* v**(g)
      YFS(g,v) = exp( -g*log(1d0-v) +0.25d0*g)
      YFSd(g) = exp(0.25d0*g)
! ---------------------------------------------------------------------
      DATA fleps / 1d-6 /
! ---------------------------------------------------------------------
!
      XMAX=  1-ZMIN
      z1  =1-x1
      z2  =1-x2
      z3  =1-x3
      z4  =1-x4
!-----------------------------------------------------------------------
! Only one type of the big-log definition alowed
      sprim = cmsene**2*z1*z2
      beta  = alf1*(dlog(sprim*xi/amel**2) -1)
      gam   = 2d0*beta
!-----------------------------------------------------------------------
! Born Factor
      WTBOR = chi(xi)
      WTfin = wtfin3*wtfin4
!-----------------------------------------------------------------------
      vI  = x1*z3
      vF  = x3
      uI  = x2*z4
      uF  = x4
      v   = vI+vF
      u   = uI+uF
      FFac = FF(beta)
!-----------------------------------------------------------------------
      bt00  =0
      bt01  =0
      bt02  =0
      bt10u =0
      bt10l =0
      bt11u =0
      bt11l =0
      bt20u =0
      bt20l =0
      bt2ul =0
!=======================================================================
!               UPPER LINE   only   INITIAL+FINAL 
!=======================================================================
!-----------------------------------------------------------------------
! Crude distribution, 
! x**(beta1-1) divided out for the sake of numerical stability!
!-----------------------------------------------------------------------
!===> CRU1=   chi(x1) *eta *x1**(eta-1)
!===> CRU3=   chi(x3) *eta *x3**(eta-1)
!===> CRU2=   chi(x2) *eta *x2**(eta-1)
!===> CRU4=   chi(x4) *eta *x4**(eta-1)
      XRU1=   chi(x1) *eta *xmax**( eta-beta)*tt1**(1-beta/eta)
      XRU3=   chi(x3) *eta                   *tt3**(1-beta/eta)
      XRU2=   chi(x2) *eta *xmax**( eta-beta)*tt2**(1-beta/eta)
      XRU4=   chi(x4) *eta                   *tt4**(1-beta/eta)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                         beta zero
!-----------------------------------------------------------------------
      DDCRU =  XRU1*XRU3*XRU2*XRU4
      bt00 =    YFS(2*beta,vF)
!===>$     *FFac *beta *(x1*z3)**(beta-1) *z3   ! Upper Initial
!===>$     *FFac *beta *   (x3)**(beta-1)       ! Upper Final
!===>$           *beta *   (x2)**(beta-1)    ! Lower Initial
!===>$           *beta *   (x4)**(beta-1)    ! Lower Final
     $     *FFac *beta *(   z3)**(beta-1) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $           *beta                       ! Lower Initial
     $           *beta                       ! Lower Final
     $     /DDCRU
      bt01  = bt00*(1 + gam/2            )
      bt02  = bt00*(1 + gam/2  +1/8d0*gam**2)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                         beta one
!-----------------------------------------------------------------------
! Contributions from beta1 UPPER INITIAL line
!-----------------------------------------------------------------------
      bt11uI = 0
      IF(vI .GT. fleps ) THEN
         CRU1   =  chi(x1) *eta *x1**(eta-1)
         DDCRU  =  CRU1*XRU3*XRU2*XRU4
         bt11uI =  YFS(2*beta,vF)
     $     *FFac *vI**beta  *RRI2u(gam,vI) *z3   ! Upper Initial
     $     *FFac *beta                           ! Upper Final
     $           *beta                       ! Lower Initial
     $           *beta                       ! Lower Final
     $        /DDCRU
      ENDIF
!-----------------------------------------------------------------------
! Contributions from beta1 UPPER FINAL line
!-----------------------------------------------------------------------
      bt11uF = 0
      IF(vF .GT. fleps ) THEN
         CRU3   =  chi(x3) *eta *x3**(eta-1)
         DDCRU  =  XRU1*CRU3*XRU2*XRU4
         bt11uF =  YFS(2*beta,vF)
     $     *FFac *beta *(   z3)**(beta-1) *z3  ! Upper Initial
     $     *FFac *vF**beta  *RRF2u(gam,vF)     ! Upper Final
     $           *beta                       ! Lower Initial
     $           *beta                       ! Lower Final
     $        /DDCRU
      ENDIF
      bt11u = bt11uI +bt11uF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                       beta two Upper Upper
!-----------------------------------------------------------------------
      bt2uuIF = 0
!-- Upper Initial * Upper Final
      IF(vI .GT. fleps .AND. vF .GT. fleps ) THEN
         CRU1   =  chi(x1) *eta *x1**(eta-1)
         CRU3   =  chi(x3) *eta *x3**(eta-1)
         DDCRU  =  CRU1*CRU3*XRU2*XRU4
         bt2uuIF =   YFS(2*beta,vF)
     $        *FFac *vI**beta                *z3 ! Upper Initial
     $        *FFac *vF**beta                    ! Upper Final
     $   *1d0/4*gam**2*( xkIF(vI,vF) +gam/2*wIF(vI,vF) )
     $           *beta                       ! Lower Initial
     $           *beta                       ! Lower Final
     $        /DDCRU
      ENDIF
      bt2uuII = 0
!-- Upper Initial * Upper Initial
      IF(vI .GT. fleps ) THEN
         CRU1   =  chi(x1) *eta *x1**(eta-1)
         DDCRU  =  CRU1*XRU3*XRU2*XRU4
         bt2uuII  =  YFS(2*beta,vF)
     $     *FFac *vI**beta  *b2II(gam,vI) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $           *beta                       ! Lower Initial
     $           *beta                       ! Lower Final
     $        /DDCRU
      ENDIF
      bt2uuFF = 0
!-- Upper Final * Upper Final
      IF(vF .GT. fleps ) THEN
         CRU3   =  chi(x3) *eta *x3**(eta-1)
         DDCRU  =  XRU1*CRU3*XRU2*XRU4
         bt2uuFF =  YFS(2*beta,vF)
     $     *FFac *beta *(   z3)**(beta-1) *z3 ! Upper Initial
     $     *FFac *vF**beta  *b2II(gam,vF)     ! Upper Final
     $           *beta                       ! Lower Initial
     $           *beta                       ! Lower Final
     $        /DDCRU
      ENDIF
      bt20u=bt2uuIF+ bt2uuII + bt2uuFF


!=======================================================================
!-----------------------------------------------------------------------
!                   Define MC weights
!-----------------------------------------------------------------------
!=======================================================================
      Facton =  WTBOR *WTfin
!           ---------------------------------
!           ! EXPONentated  UPPER line only !
!           ---------------------------------
! All beta's:  O(alf0),O(alf1),O(alf2)
      WTSET(130) =    bt00  *Facton
      WTSET(131) =   (bt01+bt10u)*Facton 
      WTSET(132) =   (bt02+bt11u+bt20u)*Facton
! Individual beta's in various orders.
! O(alf1)
      WTSET(135) =    bt01*Facton                   ! beta0
      WTSET(136) =    bt10u*Facton                  ! beta1
! O(alf2)
      WTSET(147) =    bt02*Facton                   ! beta0
      WTSET(148) =    bt11u*Facton                  ! beta1
      WTSET(149) =    bt20u*Facton                  ! beta2
!-----------------------------------------------------------------------
      END

      FUNCTION RRF2u(gam,V)
!     ***********************************
! Version for upper line only, adequate one line virtual in r_1 !!!
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(CEULER =0.57721566D0)
!-----------------------------------------------------------------------
!   [R_F^2(V)] =
!       + [N]*[1]^3*gam^3*[ln(1-V)]^2 * (  - 1/32 + 1/64*V + 1/32*V^-1 )
!       + [N]*[1]^3*gam^3*[ln(1-V)] * ( 1/32 - 1/32*V )
!       + [N]*[1]^3*gam^3*[Li2(V)] * (  - 1/16 + 1/32*V )
!       + [N]*[1]^3*gam^3 * ( 1/32*V )
!       + [N]*[1]^2*gam^2*[ln(1-V)] * (  - 1/8 + 1/16*V + 1/8*V^-1 )
!       + [N]*[1]^2*gam^2 * (  - 1/4 )
!       + [N]*[1]*gam * (  - 1/2 + 1/4*V );
!-----------------------------------------------------------------------
      RRF2u =
     $ +gam**3*log(1-V)**2   * (  - 1d0/32 + 1d0/64*V + 1d0/32/V )
     $ +gam**3*log(1-V)      * ( 1d0/32    - 1d0/32*V )
     $ +gam**3*dilogy(V)     * (  - 1d0/16 + 1d0/32*V )
     $ +gam**3               * (             1d0/32*V )
     $ +gam**2*log(1-V)      * (  - 1d0/8 + 1d0/16*V + 1d0/8/V )
     $ +gam**2               * (  - 1d0/4 )
     $ +gam                  * (  - 1d0/2 + 1d0/4*V )
      END
      FUNCTION RRI2u(gam,V)
!     ***********************************
! Version for upper line only, adequate one line virtual in r_2 !!!
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(CEULER =0.57721566D0)
!-----------------------------------------------------------------------
!   [R_I^2(V)] =
!       + [N]*[1]^3*gam^3*[ln(1-V)]^2 * ( 1/32 - 1/64*V - 1/32*V^-1 )
!       + [N]*[1]^3*gam^3*[ln(1-V)] * (  - 1/32 + 1/32*V )
!       + [N]*[1]^3*gam^3*[Li2(V)] * ( 1/16 - 1/32*V )
!       + [N]*[1]^3*gam^3 * (  - 1/32*V )
!       + [N]*[1]^2*gam^2*[ln(1-V)] * ( 1/8 - 1/16*V - 1/8*V^-1 )
!       + [N]*[1]^2*gam^2 * (  - 1/4 )
!       + [N]*[1]*gam * (  - 1/2 + 1/4*V );
!-----------------------------------------------------------------------
      RRI2u =
     $ +gam**3*log(1-V)**2   * ( 1d0/32 - 1d0/64*V - 1d0/32/V )
     $ +gam**3*log(1-V)      * (  - 1d0/32 + 1d0/32*V )
     $ +gam**3*dilogy(V)     * ( 1d0/16 - 1d0/32*V )
     $ +gam**3               * (  - 1d0/32*V )
     $ +gam**2*log(1-V)      * ( 1d0/8 - 1d0/16*V - 1d0/8/V )
     $ +gam**2               * (  - 1d0/4 )
     $ +gam                  * (  - 1d0/2 + 1d0/4*V )
      END


      SUBROUTINE mdlbhe(wtfin3,wtfin4)
!     ********************************
! LL emulation of BHLUMI amplitude up to O(gam**3)exp.
! Explicit emission INITIAL+FINAL STATE
! ---------------------------------------------------------------------
! ===  In order to assure numerical stability
! ===  factors x**(beta-1) =[xmax*t**(1/eta)]**(beta-1) removed.
! ===  Similar division of x**(eta-1) is done in FUNSKI.
! ---------------------------------------------------------------------
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(ALFA   =  1D0/137.03604D0,   ALF1=  ALFA/PI)
      PARAMETER(AMEL   =  0.511D-3)       
      PARAMETER(CEULER =0.57721566D0)
      COMMON / PARFUN / ETA,ZMIN,XIA,XIB,TT1,TT2,TT3,TT4,KEYTES
      COMMON / UTIBLO / X1,X2,X3,X4,XI1,XI2,XI
      COMMON / PARBLO / CMSENE,TMINL,RAXI,XK0,KEYRAD      
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      COMMON / INOUT  / NINP,NOUT 
! ---------------------------------------------------------------------
! Inline functions
! ---------------------------------------------------------------------
      chi(x)   = (1+(1-x)**2)/2
      FF(g)    = EXP( -g*ceuler )/dpgamm( 1d0 + g )
      STR(g,v) = g* v**(g)
      YFS(g,v) = exp( -g*log(1d0-v) +0.25d0*g)
! ---------------------------------------------------------------------
      DATA fleps / 1d-6 /
! ---------------------------------------------------------------------
!
      XMAX=  1-ZMIN
      z1  =1-x1
      z2  =1-x2
      z3  =1-x3
      z4  =1-x4
!-----------------------------------------------------------------------
! Only one type of the big-log definition alowed
      sprim = cmsene**2*z1*z2
      beta  = alf1*(dlog(sprim*xi/amel**2) -1)
      gam   = 2d0*beta
!-----------------------------------------------------------------------
! Born Factor
      WTBOR = chi(xi)
      WTfin = wtfin3*wtfin4
!-----------------------------------------------------------------------
      vI  = x1*z3
      vF  = x3
      uI  = x2*z4
      uF  = x4
      v   = vI+vF
      u   = uI+uF
      FFac = FF(beta)
!-----------------------------------------------------------------------
      bt00  =0
      bt01  =0
      bt02  =0
      bt10u =0
      bt10l =0
      bt11u =0
      bt11l =0
      bt20u =0
      bt20l =0
      bt2ul =0
!-----------------------------------------------------------------------
! Crude distribution, 
! x**(beta1-1) divided out for the sake of numerical stability!
!-----------------------------------------------------------------------
!===> CRU1=   chi(x1) *eta *x1**(eta-1)
!===> CRU3=   chi(x3) *eta *x3**(eta-1)
!===> CRU2=   chi(x2) *eta *x2**(eta-1)
!===> CRU4=   chi(x4) *eta *x4**(eta-1)
      XRU1=   chi(x1) *eta *xmax**( eta-beta)*tt1**(1-beta/eta)
      XRU3=   chi(x3) *eta                   *tt3**(1-beta/eta)
      XRU2=   chi(x2) *eta *xmax**( eta-beta)*tt2**(1-beta/eta)
      XRU4=   chi(x4) *eta                   *tt4**(1-beta/eta)
!=======================================================================
!                         beta zero
!=======================================================================
!-----------------------------------------------------------------------
      DDCRU =  XRU1*XRU2*XRU3*XRU4
      bt00 =    YFS(2*beta,vF) *YFS(2*beta,uF)
!===>$     *FFac *beta *(x1*z3)**(beta-1) *z3   ! Upper Initial
!===>$     *FFac *beta *   (x3)**(beta-1)       ! Upper Final
!===>$     *FFac *beta *(x2*z4)**(beta-1) *z4   ! Lower Initial
!===>$     *FFac *beta *   (x4)**(beta-1)       ! Lower Final
     $     *FFac *beta *(   z3)**(beta-1) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $     *FFac *beta *(   z4)**(beta-1) *z4   ! Lower Initial
     $     *FFac *beta                          ! Lower Final
     $     /DDCRU
      bt01  = bt00*(1 + gam            )
      bt02  = bt00*(1 + gam +0.5*gam**2)
!=======================================================================
!                         beta one
!=======================================================================
!-----------------------------------------------------------------------
! Contributions from beta1 UPPER INITIAL line
!-----------------------------------------------------------------------
      bt11uI = 0
      IF(vI .GT. fleps ) THEN
         CRU1   =  chi(x1) *eta *x1**(eta-1)
         DDCRU  =  CRU1*XRU2*XRU3*XRU4
         bt11uI =  YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *vI**beta  *RRI2(gam,vI) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $     *FFac *beta *(   z4)**(beta-1) *z4   ! Lower Initial
     $     *FFac *beta                          ! Lower Final
     $        /DDCRU
      ENDIF
!-----------------------------------------------------------------------
! Contributions from beta1 UPPER FINAL line
!-----------------------------------------------------------------------
      bt11uF = 0
      IF(vF .GT. fleps ) THEN
         CRU3   =  chi(x3) *eta *x3**(eta-1)
         DDCRU  =  XRU1*XRU2*CRU3*XRU4
         bt11uF =  YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *beta *(   z3)**(beta-1) *z3 ! Upper Initial
     $     *FFac *vF**beta  *RRF2(gam,vF)     ! Upper Final
     $     *FFac *beta *(   z4)**(beta-1) *z4 ! Lower Initial
     $     *FFac *beta                        ! Lower Final
     $        /DDCRU
      ENDIF
      bt11u = bt11uI +bt11uF
!-----------------------------------------------------------------------
! Contributions from beta1 LOWER INITIAL line
!-----------------------------------------------------------------------
      bt11lI = 0
      IF(uI .GT. fleps ) THEN
         CRU2   =  chi(x2) *eta *x2**(eta-1)
         DDCRU  =  XRU1*CRU2*XRU3*XRU4
         bt11lI =    YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *beta *(   z3)**(beta-1) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $     *FFac *uI**beta *RRI2(gam,uI)  *z4   ! Lower Initial
     $     *FFac *beta                          ! Lower Final
     $        /DDCRU
      ENDIF
!-----------------------------------------------------------------------
! Contributions from beta1 LOWER FINAL line
!-----------------------------------------------------------------------
      bt11lF = 0
      IF(uF .GT. fleps ) THEN
         CRU4   =  chi(x4) *eta *x4**(eta-1)
         DDCRU  =  XRU1*XRU2*XRU3*CRU4
         bt11lF =    YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *beta *(   z3)**(beta-1) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $     *FFac *beta *(   z4)**(beta-1) *z4   ! Lower Initial
     $     *FFac *uF**beta     *RRF2(gam,uF)    ! Lower Final
     $        /DDCRU
      ENDIF
      bt11l = bt11lI +bt11lF
!=======================================================================
!                   beta two   Upper Lower
!=======================================================================
!-- Upper Initial * Lower Initial 
      bt2ulII = 0
      IF(vI .GT. fleps .AND. uI .GT. fleps) THEN
         CRU1   =  chi(x1) *eta *x1**(eta-1)
         CRU2   =  chi(x2) *eta *x2**(eta-1)
         DDCRU  =  CRU1*CRU2*XRU3*XRU4
         bt2ulII = YFS(2*beta,vF) *YFS(2*beta,uF)
     $        *FFac *vI**beta *RRI1(gam,vI)  *z3 ! Upper Initial
     $        *FFac *beta                        ! Upper Final
     $        *FFac *uI**beta *RRI1(gam,uI)  *z4 ! Lower Initial
     $        *FFac *beta                        ! Lower Final
     $        /DDCRU
      ENDIF
!-- Upper Final * Lower Final
      bt2ulFF = 0
      IF(vF .GT. fleps .AND. uF .GT. fleps) THEN
         CRU3   =  chi(x3) *eta *x3**(eta-1)
         CRU4   =  chi(x4) *eta *x4**(eta-1)
         DDCRU  =  XRU1*XRU2*CRU3*CRU4
         bt2ulFF =  YFS(2*beta,vF) *YFS(2*beta,uF)
     $        *FFac *beta *(   z3)**(beta-1) *z3 ! Upper Initial
     $        *FFac *vF**beta  *RRF1(gam,vF)     ! Upper Final
     $        *FFac *beta *(   z4)**(beta-1) *z4 ! Lower Initial
     $        *FFac *uF**beta  *RRF1(gam,uF)     ! Lower Final
     $        /DDCRU
      ENDIF
!-- Upper Initial * Lower Final
      bt2ulIF = 0
      IF(vI .GT. fleps .AND. uF .GT. fleps) THEN
         CRU1   =  chi(x1) *eta *x1**(eta-1)
         CRU4   =  chi(x4) *eta *x4**(eta-1)
         DDCRU  =  CRU1*XRU2*XRU3*CRU4
         bt2ulIF =  YFS(2*beta,vF) *YFS(2*beta,uF)
     $        *FFac *vI**beta  *RRI1(gam,vI) *z3 ! Upper Initial
     $        *FFac *beta                        ! Upper Final
     $        *FFac *beta *(   z4)**(beta-1) *z4 ! Lower Initial
     $        *FFac *uF**beta  *RRF1(gam,uF)     ! Lower Final
     $        /DDCRU
      ENDIF
      bt2ulFI = 0
!-- Upper Final * Lower Initial 
      IF(vF .GT. fleps .AND. uI .GT. fleps ) THEN
         CRU3   =  chi(x3) *eta *x3**(eta-1)
         CRU2   =  chi(x2) *eta *x2**(eta-1)
         DDCRU  =  XRU1*CRU2*CRU3*XRU4
         bt2ulFI =   YFS(2*beta,vF) *YFS(2*beta,uF)
     $        *FFac *beta *(   z3)**(beta-1) *z3 ! Upper Initial
     $        *FFac *vF**beta *RRF1(gam,vF)      ! Upper Final
     $        *FFac *uI**beta *RRI1(gam,uI)  *z4 ! Lower Initial
     $        *FFac *beta                        ! Lower Final
     $        /DDCRU
      ENDIF
      bt2ul = bt2ulII+bt2ulFF+bt2ulIF+bt2ulFI
!=======================================================================
!                   beta two   Upper Upper
!=======================================================================
      bt2uuIF = 0
!-- Upper Initial * Upper Final
      IF(vI .GT. fleps .AND. vF .GT. fleps ) THEN
         CRU1   =  chi(x1) *eta *x1**(eta-1)
         CRU3   =  chi(x3) *eta *x3**(eta-1)
         DDCRU  =  CRU1*XRU2*CRU3*XRU4
         bt2uuIF =   YFS(2*beta,vF) *YFS(2*beta,uF)
     $        *FFac *vI**beta                *z3 ! Upper Initial
     $        *FFac *vF**beta                    ! Upper Final
     $   *1d0/4*gam**2*( xkIF(vI,vF) +gam/2*wIF(vI,vF) )
     $        *FFac *beta *(   z4)**(beta-1) *z4 ! Lower Initial
     $        *FFac *beta                        ! Lower Final
     $        /DDCRU
      ENDIF
      bt2uuII = 0
!-- Upper Initial * Upper Initial
      IF(vI .GT. fleps ) THEN
         CRU1   =  chi(x1) *eta *x1**(eta-1)
         DDCRU  =  CRU1*XRU2*XRU3*XRU4
         bt2uuII  =  YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *vI**beta  *b2II(gam,vI) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $     *FFac *beta *(   z4)**(beta-1) *z4   ! Lower Initial
     $     *FFac *beta                          ! Lower Final
     $        /DDCRU
      ENDIF
      bt2uuFF = 0
!-- Upper Final * Upper Final
      IF(vF .GT. fleps ) THEN
         CRU3   =  chi(x3) *eta *x3**(eta-1)
         DDCRU  =  XRU1*XRU2*CRU3*XRU4
         bt2uuFF =  YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *beta *(   z3)**(beta-1) *z3 ! Upper Initial
     $     *FFac *vF**beta  *b2II(gam,vF)     ! Upper Final
     $     *FFac *beta *(   z4)**(beta-1) *z4 ! Lower Initial
     $     *FFac *beta                        ! Lower Final
     $        /DDCRU
      ENDIF
      bt20u=bt2uuIF+ bt2uuII + bt2uuFF
!=======================================================================
!                   beta two   Lower Lower
!=======================================================================
      bt2llIF = 0
!-- Lower Initial * Lower Final
      IF(uI .GT. fleps .AND. uF .GT. fleps ) THEN
         CRU2   =  chi(x2) *eta *x2**(eta-1)
         CRU4   =  chi(x4) *eta *x4**(eta-1)
         DDCRU  =  XRU1*CRU2*XRU3*CRU4
         bt2llIF =   YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *beta *(   z3)**(beta-1) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $     *FFac *uI**beta                *z4   ! Lower Initial
     $     *FFac *uI**beta                      ! Lower Final
     $   *1d0/4*gam**2*( xkIF(uI,uF) +gam/2*wIF(uI,uF) )
     $        /DDCRU
      ENDIF
      bt2llII = 0
!-- Lower Initial * Lower Initial
      IF(uI .GT. fleps ) THEN
         CRU2   =  chi(x2) *eta *x2**(eta-1)
         DDCRU  =  XRU1*CRU2*XRU3*XRU4
         bt2llII  =  YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *beta *(   z3)**(beta-1) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $     *FFac *uI**beta *b2II(gam,uI)  *z4   ! Lower Initial
     $     *FFac *beta                          ! Lower Final
     $        /DDCRU
      ENDIF
      bt2llFF = 0
!-- Lower Final * Lower Final
      IF(uF .GT. fleps ) THEN
         CRU4   =  chi(x4) *eta *x4**(eta-1)
         DDCRU  =  XRU1*XRU2*XRU3*CRU4
         bt2llFF =  YFS(2*beta,vF) *YFS(2*beta,uF)
     $     *FFac *beta *(   z3)**(beta-1) *z3   ! Upper Initial
     $     *FFac *beta                          ! Upper Final
     $     *FFac *beta *(   z4)**(beta-1) *z4   ! Lower Initial
     $     *FFac *uF**beta     *b2II(gam,uF)    ! Lower Final
     $        /DDCRU
      ENDIF
      bt20l=bt2llIF+ bt2llII + bt2llFF


!=======================================================================
!-----------------------------------------------------------------------
!                   Define MC weights
!-----------------------------------------------------------------------
!=======================================================================
      Facton =  WTBOR *WTfin
!           ---------------------------------
!           ! EXPON  UPPER + LOWER line     !
!           ---------------------------------
! All beta's:  O(alf0),O(alf1),O(alf2)
      WTSET( 40) =    bt00  *Facton
      WTSET( 41) =   (bt01+bt10u+bt10l)*Facton 
      WTSET( 42) =   (bt02+bt11u+bt11l+bt20u+bt20l+bt2ul)*Facton 
! Individual beta's in various orders.
! O(alf1)
      WTSET( 43) =    bt01*Facton                   ! beta0
      WTSET( 44) =   (bt10u+bt10l)*Facton           ! beta1
      WTSET( 45) =    bt10u*Facton                  ! beta1 components
      WTSET( 46) =    bt10l*Facton  
! O(alf2)
      WTSET( 47) =    bt02*Facton                   ! beta0
      WTSET( 48) =   (bt11u+bt11l)*Facton           ! beta1
      WTSET( 49) =   (bt20u+bt20l+bt2ul)*Facton     ! beta2
      WTSET( 50) =    bt11u*Facton                  ! beta1 components
      WTSET( 51) =    bt11l*Facton  
      WTSET( 52) =    bt2ul*Facton  
      WTSET( 53) =    bt20u*Facton                  ! beta2 components
      WTSET( 54) =    bt20l*Facton 
! Debug weights 
      WTSET( 55) =    bt2ulII*Facton
      WTSET( 56) =    bt2ulFF*Facton
      WTSET( 57) =    bt2llIF*Facton
      WTSET( 58) =    bt2llII*Facton
      WTSET( 59) =    bt2llFF*Facton
      WTSET( 60) =    bt2ulIF*Facton
      WTSET( 61) =    bt2ulFI*Facton
      WTSET( 62) =    bt11uI*Facton
      WTSET( 63) =    bt11uF*Facton
!-----------------------------------------------------------------------
      END

      FUNCTION RRI2(gam,V)
!     ***********************************
! O(alf2) Result of convolution of beta-1 contribution with near-by
! bunch of soft photons, on the same INITIAL state fermion leg.
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(CEULER =0.57721566D0)
!-----------------------------------------------------------------------
!   2*[N] = exp(DYFS(gamb,v)) *F(gamb) *[v^gamb]*b_0
!   [R_I^2(V)] =
!       + [N]*[1]^3*gam^3*[ln(1-V)]^2 * ( 1/32 - 1/64*V - 1/32*V^-1 )
!       + [N]*[1]^3*gam^3*[ln(1-V)] * (  - 1/32 + 1/32*V )
!       + [N]*[1]^3*gam^3*[Li2(V)] * ( 1/16 - 1/32*V )
!       + [N]*[1]^3*gam^3 * (  - 3/32*V )
!       + [N]*[1]^2*gam^2*[ln(1-V)] * ( 1/8 - 1/16*V - 1/8*V^-1 )
!       + [N]*[1]^2*gam^2 * (  - 1/2 + 1/8*V )
!       + [N]*[1]*gam * (  - 1/2 + 1/4*V );
!-----------------------------------------------------------------------
      RRI2 =
     $  +gam**3*log(1-V)**2 * ( 1d0/32 - 1d0/64*V - 1d0/32/V )
     $  +gam**3*log(1-V)    * (  - 1d0/32 + 1d0/32*V )
     $  +gam**3*dilogy(V)   * ( 1d0/16 - 1d0/32*V )
     $  +gam**3             * (  - 3d0/32*V )
     $  +gam**2*log(1-V)    * ( 1d0/8 - 1d0/16*V - 1d0/8/V )
     $  +gam**2             * (  - 1d0/2 + 1d0/8*V )
     $  +gam                * (  - 1d0/2 + 1d0/4*V )
      END
      FUNCTION RRI1(gam,V)
!     ***********************************
! O(alf1) Result of convolution...
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
!-----------------------------------------------------------------------
!   [R_I^1(V)] =
!       + [N]*[1]^3*gam^3 * (         1/16*V )
!       + [N]*[1]^2*gam^2 * (        - 1/8*V )
!       + [N]*[1]*gam     * (  - 1/2 + 1/4*V )
!-----------------------------------------------------------------------
      RRI1 =
     $     +gam**3     *(        +1d0/16*V )
     $     +gam**2     *(        - 1d0/8*V )
     $     +gam        *( -1d0/2 + 1d0/4*V )
      END


      FUNCTION RRF2(gam,V)
!     ***********************************
! O(alf2) Result of convolution of beta-1 contribution with near-by
! bunch of soft photons, on the same FINAL state fermion leg.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      SAVE
      PARAMETER(PI     =  3.1415926535897932D0) 
      PARAMETER(CEULER =0.57721566D0)
!-----------------------------------------------------------------------
!   2*[N] = exp(DYFS(gamb,v)) *F(gamb) *[v^gamb]*b_0
!   [R_F^2(V)] =
!       + [N]*[1]^3*gam^3*[ln(1-V)]^2 * (  - 1/32 + 1/64*V + 1/32*V^-1 )
!       + [N]*[1]^3*gam^3*[ln(1-V)] * ( 1/32 - 1/32*V )
!       + [N]*[1]^3*gam^3*[Li2(V)] * (  - 1/16 + 1/32*V )
!       + [N]*[1]^3*gam^3 * (  - 1/32*V )
!       + [N]*[1]^2*gam^2*[ln(1-V)] * (  - 1/8 + 1/16*V + 1/8*V^-1 )
!       + [N]*[1]^2*gam^2 * (  - 1/2 + 1/8*V )
!       + [N]*[1]*gam * (  - 1/2 + 1/4*V );
!-----------------------------------------------------------------------
      RRF2 =
     $  +gam**3*log(1-V)**2 * (  - 1d0/32 + 1d0/64*V + 1d0/32/V )
     $  +gam**3*log(1-V)    * ( 1d0/32 - 1d0/32*V )
     $  +gam**3*dilogy(V)   * (  - 1d0/16 + 1d0/32*V )
     $  +gam**3             * (  - 1d0/32*V )
     $  +gam**2*log(1-V)    * (  - 1d0/8 + 1d0/16*V + 1d0/8/V )
     $  +gam**2             * (  - 1d0/2 + 1d0/8*V )
     $  +gam                * (  - 1d0/2 + 1d0/4*V )
      END
      FUNCTION RRF1(gam,V)
!     ***********************************
! O(alf1) Result of convolution...
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
!-----------------------------------------------------------------------
!   [R_F^1(V)] =
!       + [N]*[1]^3*gam^3 * ( 1/16*V )
!       + [N]*[1]^2*gam^2 * (  - 1/8*V )
!       + [N]*[1]*gam * (  - 1/2 + 1/4*V );
!-----------------------------------------------------------------------
      RRF1 =
     $     +gam**3     *(        +1d0/16*V )
     $     +gam**2     *(        - 1d0/8*V )
     $     +gam        *( -1d0/2 + 1d0/4*V )
      END


      FUNCTION wIF(v1,v2)
!     ***********************************
! Element of beta2UU calculation
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
**!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
**!%        [wIF(v1,v2)] =
**!%       + v1*[1-v2]^-1 * (  - 1/8 )
**!%       + v1*[1-v2]^-2 * (  - 1/8 )
**!%       + v2*[1-v1]^-1 * ( 1/8 )
**!%       + v2*[1-v1]^-2 * (  - 1/8 )
**!%       + ln(1 - v1)*v2*[1-v1]^-1 * ( 1/8 )
**!%       + ln(1 - v1)*v2*[1-v1]^-2 * ( 1/8 )
**!%       + ln(1 - v1)*[1-v1]^-1 * (  - 1/4 );
**!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
**      wIF=
**     $ + v1/(1-v2)    * (  - 1d0/8 )
**     $ + v1/(1-v2)**2 * (  - 1d0/8 )
**     $ + v2/(1-v1)    * ( 1d0/8 )
**     $ + v2/(1-v1)**2 * (  - 1d0/8 )
**     $ + log(1 - v1)*v2/(1-v1)    * ( 1d0/8 )
**     $ + log(1 - v1)*v2/(1-v1)**2 * ( 1d0/8 )
**     $ + log(1 - v1)   /(1-v1)    * (-1d0/4 )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%   [wIF(v1,v2)] =
!%       +            v1*[1-v2]^-2 * (-1/2 )
!%       + ln(1 - v2)   *[1-v2]^-1 * (-1/2 );
!%       + ln(1 - v2)*v1*[1-v2]^-1 * ( 1/4 )
!%       + ln(1 - v2)*v1*[1-v2]^-2 * ( 1/4 )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      wIF=
     $ +             v1/(1-v2)**2 * (-1d0/2 )
     $ + log(1 - v2)   /(1-v2)    * (-1d0/2 )
     $ + log(1 - v2)*v1/(1-v2)    * ( 1d0/4 )
     $ + log(1 - v2)*v1/(1-v2)**2 * ( 1d0/4 )
      END

      FUNCTION xkIF(v1,v2)
!     ***********************************
! Element of beta2UU calculation
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
**!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
**!%   [k(v1,v2)] =
**!%       + v1*[1-v2]^-1 * ( 1/8 )
**!%       + v1*[1-v2]^-2 * ( 1/8 )
**!%       + v2*[1-v1]^-1 * ( 1/8 )
**!%       + v2*[1-v1]^-2 * ( 1/8 )
**!%       + [1-v1]^-1 * (  - 1/4 )
**!%       + [1-v2]^-1 * (  - 1/4 )
**!%       + 1/2;
**!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
**      xkIF=
**     $ + v1/(1-v2)    * ( 1d0/8 )
**     $ + v1/(1-v2)**2 * ( 1d0/8 )
**     $ + v2/(1-v1)    * ( 1d0/8 )
**     $ + v2/(1-v1)**2 * ( 1d0/8 )
**     $ + 1/(1-v1)     * (-1d0/4 )
**     $ + 1/(1-v2)     * (-1d0/4 )
**     $ + 1d0/2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%   [kIF(v1,v2)] =
!%       +                  1/2;
!%       +    [1-v2]^-1 * (-1/2 )
!%       + v1*[1-v2]^-1 * ( 1/4 )
!%       + v1*[1-v2]^-2 * ( 1/4 )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      xkIF=
     $ +                  1d0/2
     $ +  1/(1-v2)    * (-1d0/2 )
     $ + v1/(1-v2)    * ( 1d0/4 )
     $ + v1/(1-v2)**2 * ( 1d0/4 )
      END

      FUNCTION b2II(gam,v)
!     ***********************************
! Element of beta2UU calculation
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%   [d_0^*(v)] =
!%       + ln(1 - v) * ( 1/2 - 1/4*v )
!%       + 1/2*v;
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      dd0 =  log(1 - v) * ( 1d0/2 - 1d0/4*v )+ 1d0/2*v
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%   [d_1^*(v)] =
!%       + ln(1 - v) * (  - 1/4 + 1/4*v )
!%       + Li2( - 1/(1 - v)*v) * (  - 1/2 + 1/4*v )
!%       - 3/4*v;
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      dd1 =
     $     + log(1 - v)               * (  - 1d0/4 + 1d0/4*v )
     $     + dilogy( - 1d0/(1 - v)*v) * (  - 1d0/2 + 1d0/4*v )
     $     - 3d0/4*v
!------
      b2II = 1d0/8 *gam**2 *(dd0 +gam/2*dd1)
      END
