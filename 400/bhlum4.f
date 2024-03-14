!
      SUBROUTINE BHLUM4(MODE,XPAR,NPAR)  
!     *********************************  
! 
!           ************************************************** 
!           *       **********************************       *
!           *       *      *******************       *       *
!           *       *      *                 *       *       *
!           *       *      *   B H L U M 4   *       *       *
!           *       *      *                 *       *       *
!           *       *      *******************       *       *
!           *       **********************************       *
!           **************************************************
!
! =======================     AUTHORS      =============================
! ==  S. Jadach, W. Placzek, E. Richter-Was, B.F.L. Ward  and Z. Was  ==
! =================  vers. 4.04 September 1996 =========================
! 
! Main subprogram of the  Monte Carlo multiphoton t-channel generator 
! BHLUMI version 4.04.
! It is multiphoton generator with Yennie-Frautschi-Suura second 
! order exponentiation based on refs. [1,2,3,4,5,6] and other.
! According to ref. [6] it features an OVERALL PRECISION is 0.11%, 
! see there for the validity range of the above statement.
! Z-contribution now added according to ref. [5].
! [1] S. Jadach, E. Richter-Was, B.F.L. Ward  and Z. Was,
!     Comp. Phys. Commun. 70 (1992) 305 (TH-6230, sept. 1991).
! [2] S. Jadach, E. Richter-Was, B.F.L. Ward  and Z. Was,
!     Phys.Lett. B268 (1991), 253 (TH-6118, June 1991). 
! [3] S. Jadach and B.F.L. Ward,
!     Phys. Rev. D40 (1989) 3582.
! [4] S. Jadach, E. Richter-Was, B.F.L. Ward  and Z. Was,
!     TH-95-38 (February 1995); Phys. Lett., B353 (1995) 362.
! [5] S. Jadach, W. Placzek and B.F.L. Ward
!     TH-95-74 (March 1995); Phys. Lett. B353 (1995) 349.
! [6] S. Jadach and O. Nicrosini (conveners of Bhabha Working Group),  
!     in Physics at LEP2, CERN Yellow Report 96-01, edited by 
!     G. Altarelli, T. Sjostrand, and F. Zwirner (CERN, Geneva, 1996), 
!     Vol. 2, p. 229, (hep-ph/9602393);
!     Summary paper of LEP2 Bhabha working group, 
!     (Convenors S. Jadach and O. Nicrosini)
!     A. Arbuzov et.al., Phys. Lett. B383 (1996) 238 (hep-ph/9605239)
! Postscript files for all the above papers are
! available from http://hpjmiady.ifj.edu.pl/programs/programs.html
! 
!----------------------------------------------------------------------
!                 INPUT and OUTPUT of BHLUM4
!----------------------------------------------------------------------
! All input and output goes through parameters in 
!                 CALL BHLUM4(MODE,XPAR,NPAR)
! and through /MOMSET/ and /WGTALL/ common blocks.
! In the following we shall  briefly indicate the meaning of the
! above parameters/variables.
! 
! IF( MODE =-1 ) THEN
! ===================
! Initialization is performed, all input parameters are transfered
! through XPAR and NPAR.
! In the following table we indicate the meaning of NPAR, XPAR 
! entries for LUMLOG subgenerator in the initialization mode. 
!      Table,           Input parameters of BHLUM4
!----------------------------------------------------------------------
!  Entry    Variable   Meaning
!----------------------------------------------------------------------
!  NPAR( 1)  KeyOpt =1000*KeyGen +100*KeyREM +10*KeyWGT +KeyRND 
!                    General Option switch in wchich:
!            KeyGen =3 for this sub-generator
!            KeyRem =0,1 removal/no-removal switch, both are OK,
!                   =1 no-removal technicaly simpler/safer
!                   =0 OBLIGATORY for KeyZet =1 !!!
!            KeyRnd =1,2 type of random number generator RANMAR,RANECU
!                   =1 better for parallel production on computer farm.
!            KeyWgt =0,1,2 for constant/variable weight WTMOD,
!                   =0, WTMOD =1 useful for apparatus Monte Carlo.
!                   =1, WTMOD variable, option faster/safer, RECOMMENDED
!                   =2, WTMOD variable, events below trmin generated
!  NPAR( 2)  KeyRad =1000*KeyZet+100*KeyUpd+10*KeyMod +KeyPia 
!                   Switch defining type of matrix element:
!            KeyZet =0,1 test switch,
!                   =0   Z contribution OFF 
!                   =1   Z contribution ON, DEFAULT!!!
!            KeyUpd =0,1,2 test switch, 
!                   =0 normal position DEFAULT!!!
!                   =1 only upper line bremss., =2 only lower line
!            KeyMod =1,2 type of MODEL subrogram and QED matrix element
!                   =1 version compatible with Comp. Phys. Comm. 70 (1992) 305
!                   =2 version 4.x which is now DEFAULT!
!            KeyPia =0,1,2,3 photon vacuum polarization and s-chanel photon
!                   =0 OFF, it used in semianalytical tests,
!                   =1 ON,  Burkhardt et.al. 1989, as in BHLUMI 2.0x
!                   =2 ON,  S. Eidelman, F. Jegerlehner, Z. Phys. C (1995)
!                   =3 ON,  Burkhardt and Pietrzyk 1995 (Moriond).
!  XPAR( 1)  CMSENE Total center mass energy [GeV]
!  XPAR( 2)   TRMIN Minimum transfer (positive) [GeV**2] 
!  XPAR( 3)   TRMAX Maximum transfer (positive) [GeV**2] 
!  XPAR( 4)   EPSCM Dimensionless infrared cut on CMS energy of soft
!                   photons, ( E_phot > CMSENE*EPSCM/2 )
!----------------------------------------------------------------------
! 
! ELSE IF( MODE = 0 ) THEN
! ========================
! Generation of the single Monte Carlo event. 
! The four momenta of the final state electron, positron and photon
! and of real photons are encoded in 
!      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
! where P1 and Q1 are four-momenta of positron and electron beams.
! P2 and Q2 are four-momenta of outgoing positron and electron.
! The list PHOT(100,4) four-momenta contains 
! NPHOT four-momenta of real the photons, all in GeV units.
! The principal weight WTM of the event is placed in
!      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)   
! It is often of interest to use 'parallel weights' from WTSET.
! The event weight is constructed then as WT= WTCRU1*WTCRU2*WTSET(J).
! Which J is allowed and what version of the QED matrix element 
! it represents is summarized in the tables below.
!
! One of the weights in WTSET (the best) is used for WTMOD.
! This is the DEFAULT weight used for rejection for KeyWgt =0.
! WARNING: using 'parallel weights' makes only sense for calculation
!          with variable weights, KeyWgt=1,2 !!!
!======================================================================!
!======================================================================!
! The list of actual weights corresponding to  KeyMod = 2              !
! (the obsolete weights for KeyMod = 1 listed later below              !
!======================================================================!
!  WtMod = wtset(1) = wtset(142)+wtset(12)  for Z on, Vac.Pol. on
!  WtMod = wtset(1) = wtset(142)+wtset(11)  for Z on, Vac.Pol. off
!        where wtset(142) is calculated in modl2b.f
!        and   wtset(11,12) are from m2agzi.f
!  The s-chan contribution and VacPol are introduced 
!  multiplicatively wtset(i)=wtset(i)*(1 +dels)*vpfac in model.f
!----------------------------------------------------------------------!
! (A) List of auxiliary weights from modl2a
!    WTSET( 30) =   O(alf0)exp      total    upper line only
!    WTSET( 31) =     O(alf1)exp    total    upper line only
!    WTSET( 32) =       O(alf2)exp  total    upper line only
!    WTSET( 35) =     O(alf1)exp    beta0    upper line only
!    WTSET( 36) =     O(alf1)exp    beta1    upper line only
!    WTSET( 37) =       O(alf2)exp  beta0    upper line only
!    WTSET( 38) =       O(alf2)exp  beta1    upper line only
!    WTSET( 39) =       O(alf2)exp  beta2    upper line only
!    WTSET( 40) =   O(alf0)exp       total    upper + lower
!    WTSET( 41) =     O(alf1)exp     total    upper + lower
!    WTSET( 42) =       O(alf2)exp   total    upper + lower
!    WTSET( 43) =     O(alf1)exp     beta0    upper + lower
!    WTSET( 44) =     O(alf1)exp     beta1    upper + lower
!    WTSET( 45) =     O(alf1)exp     beta1    upper component
!    WTSET( 46) =     O(alf1)exp     beta1    lower component
!    WTSET( 47) =       O(alf2)exp   beta0    upper + lower
!    WTSET( 48) =       O(alf2)exp   beta1    upper + lower
!    WTSET( 49) =       O(alf2)exp   beta2    upper+lower+upper*lower
!    WTSET( 50) =       O(alf2)exp   beta1    upper component
!    WTSET( 51) =       O(alf2)exp   beta1    lower component
!    WTSET( 52) =       O(alf2)exp   beta2    upper*lower component
!    WTSET( 53) =       O(alf2)exp   beta2    upper component
!    WTSET( 54) =       O(alf2)exp   beta2    lower component
!    WTSET( 60) =    O(alf0)       total    upper line only
!    WTSET( 61) =      O(alf1)     total    upper line only
!    WTSET( 62) =        O(alf2)   total    upper line only
!    WTSET( 70) =    O(alf0)       total    upper + lower
!    WTSET( 71) =      O(alf1)     total    upper + lower
!    WTSET( 72) =        O(alf2)   total    upper + lower
!----------------------------------------------------------------------!
! (B) List of auxiliary weights from modl2b
!    WTSET(130) =   O(alf0)exp      total    upper line only
!    WTSET(131) =     O(alf1)exp    total    upper line only
!    WTSET(132) =       O(alf2)exp  total    upper line only
!    WTSET(135) =     O(alf1)exp    beta0    upper line only
!    WTSET(136) =     O(alf1)exp    beta1    upper line only
!    WTSET(137) =       O(alf2)exp  beta0    upper line only
!    WTSET(138) =       O(alf2)exp  beta1    upper line only
!    WTSET(139) =       O(alf2)exp  beta2    upper line only
!    WTSET(140) =   O(alf0)exp       total    upper + lower
!    WTSET(141) =     O(alf1)exp     total    upper + lower
!    WTSET(142) =       O(alf2)exp   total    upper + lower
!    WTSET(143) =     O(alf1)exp     beta0    upper + lower
!    WTSET(144) =     O(alf1)exp     beta1    upper + lower
!    WTSET(145) =     O(alf1)exp     beta1    upper component
!    WTSET(146) =     O(alf1)exp     beta1    lower component
!    WTSET(147) =       O(alf2)exp   beta0    upper + lower
!    WTSET(148) =       O(alf2)exp   beta1    upper + lower
!    WTSET(149) =       O(alf2)exp   beta2    upper+lower+upper*lower
!    WTSET(150) =       O(alf2)exp   beta1    upper component
!    WTSET(151) =       O(alf2)exp   beta1    lower component
!    WTSET(152) =       O(alf2)exp   beta2    upper*lower component
!    WTSET(153) =       O(alf2)exp   beta2    upper component
!    WTSET(154) =       O(alf2)exp   beta2    lower component
!    WTSET(160) =    O(alf0)       total    upper line only
!    WTSET(161) =      O(alf1)     total    upper line only
!    WTSET(162) =        O(alf2)   total    upper line only
!    WTSET(170) =    O(alf0)       total    upper + lower
!    WTSET(171) =      O(alf1)     total    upper + lower
!    WTSET(172) =        O(alf2)   total    upper + lower
!----------------------------------------------------------------------!
!
!======================================================================!
!======================================================================!
! These weights are OBSOLETE and are still accessible for KeyMod = 1   !
! seeabove  for the actual weights for KeyMod = 2                      !
!======================================================================!
!    The table explains weights as published in BHLUMI Version 2.01    !
!    in  Comp. Phys. Commun. {\bf 70} (1992) 305, table 4 there.       !
!----------------------------------------------------------------------!
!  Entry      Type of QED calculation                                  !
!----------------------------------------------------------------------!
!  WTSET( 1)  = WTSET(11)                                              !
!  WTSET( 2)  = WTSET(12)           (principal weight)                 !
!             ---------------------------------------------------------!
!             QED order   Vacuum pol.   Z-exchange    s-chan.gamma exch!
!             ---------------------------------------------------------!
!  WTSET(11)  Zero-th      Yes             Yes            Yes          !
!  WTSET(12)  First        Yes             Yes            Yes          !
!  WTSET(51)  Zero          No              No             No          !
!  WTSET(52)  First         No              No             No          !
!                           ----- Miscelanous ----          !
!  WTSET(20)  First        Yes              No             No          !
!  WTSET(21)  = WTSET(20)-WTSET(52)                                    !
!  WTSET(22)  First        Yes             Yes             No          !
!  WTSET(23)  = WTSET(22)-WTSET(20)                                    !
!  WTSET(24)  First        Yes             Yes            Yes          !
!  WTSET(25)  = WTSET(24)-WTSET(22)                                    !
!             ---------------------------------------------------------!
!  WTSET(26)  Beta_0 component in WTSET(20)                            !
!  WTSET(27)  Beta_1 component in WTSET(20),                           !
!             i.e. WTSET(20)=WTSET(26)+WTSET(27)                       !
!             ---------------------------------------------------------!
!  WTSET(61)  Kind of  LL  component in    WTSET(51)                   !
!  WTSET(62)  Kind of  NLL component in    WTSET(52)                   !
!======================================================================!
!
! 
! ELSE IF( MODE = 1 ) THEN
! ========================
! The total cross section corresponding to generated series of event,
! i.e. resulting from MC integration is calculated and stored in XPAR
! and NPAR, see table below.
!----------------------------------------------------------------------
!  Entry    Variable   Meaning
!----------------------------------------------------------------------
!  NPAR(10)  NEVGEN  Number of generated MC events
!  NPAR(20)  NEVGEN  Number of generated MC events
!  XPAR(10)   XMCNB  Total x-section [nb]
!  XPAR(11)    EREL  The relative error of XPAR(10)
!  XPAR(12)     XMC  Total x-section [GEV**-2]
!  XPAR(20)  SIG0NB  Crude total MC x-section [nb] which is necessary
!                    for rescaling histograms in run with 
!                    weighted events.
!  XPAR(21)          =0, error of XPAR(20) is zero
!  XPAR(20)    SIG0  Crude x-sectio as XPAR(20) but in [GeV**-2]
!----------------------------------------------------------------------
! For constant weight option KeyWgt=0 (convenience in rescaling histos)
! we put XPAR(20,21,22)=XPAR(10,11,12) !
! For MODE=1 program is called upon many times in the process of 
! rescaling histograms, therefore, there is no output printed 
! in this mode.
! 
! ELSE IF( MODE = 2 ) THEN
! ========================                     
! Only in this MODE=2 in addition to filling XPAR and NPAR as for 
! MODE=1 the values of various x-sections are printed on standard 
! output file.
!                
! ENDIF
! ====
!*******************************************************************
!      For hidded dip-switches adjustable by hand in the code see:
!&&&&  definition of del in terms of epscm
!%%%%  TRMX2 = TRMAX !!! only for tests with fixed TRAN 
!*******************************************************************
!     *********************************   

      IMPLICIT REAL*8(A-H,O-Z)   
      INTEGER           m_EvtUnit, m_WriteLHE
      PARAMETER( PI = 3.1415926535897932D0, ALFINV = 137.03604D0)
      PARAMETER( ALFPI=  1D0/PI/ALFINV ,ALFA=1D0/ALFINV)
      PARAMETER( GNANOB=389.385D-30*1.D33 )
      CHARACTER*80   BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G 
      PARAMETER(
     $BXOPE =  '(//1X,15(5H=====)    )',
     $BXTXT =  '(1X,1H=,                  A48,25X,    1H=)',
     $BXL1I =  '(1X,1H=,I17,                 16X, A20,A12,A7, 1X,1H=)',
     $BXL1F =  '(1X,1H=,F17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2F =  '(1X,1H=,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXL1G =  '(1X,1H=,G17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2G =  '(1X,1H=,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXCLO =  '(1X,15(5H=====)/   )'    )   
      DIMENSION  XPAR(*), NPAR(*)                    
      COMMON / BHPAR1 / DEL,EPSCM,THMIN,XMIVIS        
      COMMON / BHPAR2 / CMSENE,AMEL       
      COMMON / BHPAR3 / KeyRad,KeyOpt
! CMONIT COMMUNICATES WITH Gmonit         
      COMMON / CMONIT/ AVERWT,ERRELA,NEVTOT,NEVACC,NEVNEG,NEVOVE,NEVZER 
      COMMON / TRANSR / TRAN,TRMIN,TRMAX
      COMMON / MOMS1  / TRANP,P1(4),P2(4),PHOT1(100,4),PHSU1(4),NPHOT1
      COMMON / MOMS2  / TRANQ,Q1(4),Q2(4),PHOT2(100,4),PHSU2(4),NPHOT2
      COMMON / MOMZ1  / AL1(50),BE1(50),WTM1(50),MK1(50)
      COMMON / MOMZ2  / AL2(50),BE2(50),WTM2(50),MK2(50)
      COMMON / MOMSET / PX1(4),QX1(4),PX2(4),QX2(4),PHOT(100,4),NPHOT
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)        
      COMMON / WGTSUP / WKP,WKQ,WTT1,WTT2,FPHS,FYFSU,FYFSD,WT3
      COMMON / INOUT  / NINP,NOUT
      SAVE   / BHPAR1 /, / BHPAR2 /, / BHPAR3 /, / CMONIT/, / TRANSR / 
      SAVE   / MOMS1  /, / MOMS2  /, / MOMZ1  /, / MOMZ2 /, / MOMSET /
      SAVE   / WGTALL /, / WGTSUP /, / INOUT  /
      SAVE   SVAR, WTMAX, TRMX2, EMIN, IDA,IDC
      SAVE   IDGEN, NEVGEN, IEVENT, SIG0, SIG0NB
      SAVE   KeyWgt, KeyRem, KeyUpd
      DOUBLE PRECISION DRVEC(100)
      
      m_EvtUnit = 77 

      
      IF(MODE.EQ.-1) THEN        
!     ===================
      CALL FILBH2(XPAR,NPAR)       
      
!This is Generator Identificator
      IDGEN = 3        
      SVAR=CMSENE**2 
!       
      KeyWgt = MOD(KeyOpt,100)/10   
      KeyRem = MOD(KeyOpt,1000)/100   
      KeyUpd = MOD(KeyRad,1000)/100   
      IF( KeyWgt .EQ. 2 .AND. abs(trmax/trmin-1).GT.0.1d0 ) THEN
!-- NEW feature: Generation of transfer below trmin
         sig0= 4d0*pi*alfa**2 *( 2d0/trmin-1d0/trmax)
      ELSE
!-- OLD method: no generation of transfer below trmin
         sig0= 4d0*pi*alfa**2 *( 1d0/trmin-1d0/trmax)
      ENDIF
!
      SIG0NB=SIG0*GNANOB       
      WTMAX= 2.8D0
      WTMAX= 3.0D0
! Important histo which remembers total x-section 
      CALL Gmonit(  -1, IDGEN,0D0,SIG0NB*WTMAX,1D0)          
!-- Set maximum transfer for photon angular distributions TRMX2 :
!-- In general TRMX2=svar is a safe choice. 
!-- However, for low  angles (i.e. thmin,thmax << 100mrad)
!-- It can be lowered, in order to speed up generation,
!-- but with great care!. 
      TRMX2 = TRMAX
! (over)conservative and safe choice is TRMX2=SVAR 
!     TRMX2 = SVAR   !%%%% For standard  low angles can be comented out!
      IF(TRMX2.GT.SVAR) TRMX2=SVAR
      EMIN = CMSENE/2D0*EPSCM
! Weight monitoring initialization
      IDA=50
      CALL Gmonit(  -1, IDA+1 ,0D0,1D0,1D0)  ! Obsolete        
      CALL Gmonit(  -1, IDA+2 ,0D0,1D0,1D0)  ! Obsolete           
      CALL Gmonit(  -1, IDA+3 ,0D0,1D0,1D0)  ! Obsolete         
      CALL Gmonit(  -1, IDA+4 ,0D0,1D0,1D0)  ! Obsolete          
      CALL Gmonit(  -1, IDA+18,0D0,1D0,1D0)          
      CALL Gmonit(  -1, IDA+19,0D0,1D0,1D0)          
      CALL Gmonit(  -1, IDA+20,0D0,1D0,1D0)          
      IDC = 90
      DO 12 I=IDC+1,IDC+6
  12  CALL Gmonit(  -1, I,0D0,1D0,1D0)  
      CALL GBOOK1(9001,' bhlum4, weight distribution $', 40, -1D0,7D0)
      ievent=0  
      nevgen=0
!     ======================     
      ELSEIF(MODE.EQ.0) THEN     
!     ======================     
      nevgen = nevgen+1
  200 CONTINUE
      ievent=ievent+1
      WT0=0D0
      WT1=0D0
      WT2=0D0
      WT3=0D0
      CALL VarRan(drvec,2)
!--------------------------------------------------------------------
! Generate t-channel transfer (the true one)
      rn1 = drvec(1)
      IF( KeyWgt .EQ. 2 .AND. abs(trmax/trmin-1).GT.0.1d0 ) THEN
!-- NEW feature: Generation of transfer below trmin
        ynorm= 2/trmin -1/trmax
        prob1= (1/trmin) /ynorm
        prob2= (1/trmin -1/trmax)/ynorm
        IF(rn1.LT.prob1) THEN
!-- flat distribution below trmin
          tran= trmin*rn1/prob1
          WT0 = trmin**2/tran**2
!-- to prevent crash of kinematics  (1d-8 effect)
          IF(tran.lt.(10*Amel)**2) GOTO 200
        ELSE
!-- normal 1/tra**2 distribution above trmin
          rnx  = (rn1-prob1)/prob2
          tran = 1d0 / (  rnx/trmin   +(1d0-rnx)/trmax )
          WT0  = 1d0
        ENDIF
      ELSE
!-- OLD method: no generation of transfer below trmin
        tran   = 1d0 / (  rn1/trmin   +(1d0-rn1)/trmax )
        WT0=1d0
      ENDIF
!--------------------------------------------------------------------
! Generate Photon multiplicity and momenta
      IF( KeyUpd .EQ. 0 ) THEN
        CALL mltibr(TRAN,TRMX2,AMEL,DEL,
     $      NPHOT1,PHOT1,PHSU1,AL1,BE1,TRANP,AMSP,MK1,WKP,WTM1)
        CALL mltibr(TRAN,TRMX2,AMEL,DEL,
     $      NPHOT2,PHOT2,PHSU2,AL2,BE2,TRANQ,AMSQ,MK2,WKQ,WTM2)
      ELSEIF( KeyUpd .EQ. 1) THEN
! Upper line bremss. only -- for special tests only!!
        CALL mltibr(TRAN,TRMX2,AMEL,DEL,
     $      NPHOT1,PHOT1,PHSU1,AL1,BE1,TRANP,AMSP,MK1,WKP,WTM1)
        NPHOT2 = 0
        TRANQ  = TRAN
        AMSQ   = AMEL**2
        WKQ    = 1D0
      ELSEIF( KeyUpd .eq. 2) THEN
! Lower line bremss. only -- for special tests only!!
        CALL mltibr(TRAN,TRMX2,AMEL,DEL,
     $      NPHOT2,PHOT2,PHSU2,AL2,BE2,TRANQ,AMSQ,MK2,WKQ,WTM2)
        NPHOT1 = 0
        TRANP  = TRAN
        AMSP   = AMEL**2
        WKP    = 1D0
      ELSE
        WRITE(nout,*) '+++++BHLUM4:  wrong KeyUpd=',KeyUpd
        STOP
      ENDIF
!
      IF(WKP*WKQ.EQ.0D0) GOTO 140
!--------------------------------------------------------------------
!-- Construct fermions, transform photons and fermions to CMS frame
      CALL kino4(SVAR,TRAN,AMEL,AMSP,AMSQ,WT3)
      CALL Gmonit(   0,IDA+3,  WT3, 1D0,5D-4) 
      IF(WT3.EQ.0D0) GOTO 140      
!-- Beyond this point events DO CONSERVE total four-momentum !!!
      WTKIN=1D0      
!-- Manipulations on mass weights, removal of soft photons
      CALL PIATEK(CMSENE,TRMX2,AMEL,EMIN,DEL,
     $      NPHOT1,P1,P2,PHOT1,PHSU1,WTM1,WTT1,WTMR1,WCTA1,WCTB1)
      CALL PIATEK(CMSENE,TRMX2,AMEL,EMIN,DEL,
     $      NPHOT2,Q1,Q2,PHOT2,PHSU2,WTM2,WTT2,WTMR2,WCTA2,WCTB2)
      IF( KeyRem .EQ. 0 ) THEN
!-- Removing photons < epsCM from the record 
!-- Mass weight WTMR1,2 is product of mass weights for ENE>EminCM times
!-- average weight for photons with  ENE<EminCM.
        CALL REMPHO(EMIN,NPHOT1,PHOT1,AL1,BE1,WTM1,MK1)
        CALL REMPHO(EMIN,NPHOT2,PHOT2,AL2,BE2,WTM2,MK2)
! In future the KeyUpd -logics below has to be mooved into routine Piatek
        IF( KeyUpd .EQ. 0 ) THEN
          WT1= WTMR1*WKP
          WT2= WTMR2*WKQ
        ELSEIF( KeyUpd .EQ. 1) THEN
          WT1= WTMR1*WKP
          WT2= WKQ
          WCTA2 =1D0
          WCTB2 =1D0
        ELSEIF( KeyUpd .EQ. 2) THEN
          WT1= WKP
          WT2= WTMR2*WKQ
          WCTA1 =1D0
          WCTB1 =1D0
        ENDIF
      ELSEIF( KeyRem .EQ. 1) THEN
!-- No removal, valid but non-default option. 
        WT1= WTT1*WKP
        WT2= WTT2*WKQ 
      ELSE
        WRITE(nout,*) '+++++BHLUM4:  wrong KeyRem=',KeyRem
        STOP
      ENDIF
!--------------------------------------------------------------------
!  Monitoring control weights
      CALL Gmonit(   0,IDC+1,       WCTA1, 1D0,5D-4)   
      CALL Gmonit(   0,IDC+2,       WCTA2, 1D0,5D-4)   
      CALL Gmonit(   0,IDC+3, WCTA1*WCTA2, 1D0,5D-4)        
      CALL Gmonit(   0,IDC+4,       WCTB1, 1D0,5D-4)   
      CALL Gmonit(   0,IDC+5,       WCTB2, 1D0,5D-4)   
      CALL Gmonit(   0,IDC+6, WCTB1*WCTB2, 1D0,5D-4)        
      WTM1T2 = WTMR1*WTMR2
      CALL Gmonit(   0,IDA+1,      WTM1T2,  2D0,5D-4) 
!--------------------------------------------------------------------
!  Calculate Formfactor  
      PDEL = DEL*BCUD(P1,P2,PHSU1)
      QDEL = DEL*BCUD(Q1,Q2,PHSU2) 
      FYFSU= EXP( ALFPI*(
     $  -2.D0*(DLOG(TRANP/AMEL**2)-1)*DLOG(1D0/PDEL)
     $  +0.5D0*DLOG(TRANP/AMEL**2)  -1D0
     $    ))
      FYFSD= EXP( ALFPI*(
     $  -2.D0*(DLOG(TRANQ/AMEL**2)-1)*DLOG(1D0/QDEL)
     $  +0.5D0*DLOG(TRANQ/AMEL**2)  -1D0
     $    ))
      FPHS  = EXP( 2D0*ALFPI* LOG(TRMX2/AMEL**2)* LOG(1D0/DEL ) )
      IF(     KeyUpd .EQ. 0 ) THEN
        WT4 = FPHS**2  *FYFSU *FYFSD
      ELSEIF( KeyUpd .EQ. 1) THEN
        WT4 = FPHS  *FYFSU 
      ELSEIF( KeyUpd .EQ. 2) THEN
        WT4 = FPHS  *FYFSD
      ENDIF
      CALL Gmonit(   0,IDA+4,WT4,  1D0,5D-4)  
!--------------------------------------------------------------------
!  Calculate QED Matrix element
      CALL Model(1,WT5)
!--------------------------------------------------------------------
!  Define Total Principal Weight
  140 CONTINUE      
      WT  = WT0*WT1*WT2*WT3*WT4   *WT5
!-- Monitoring model weight       
      CALL Gmonit(   0,IDA+20,WT,WTmax,RN)
      WTovr = MAX(0D0,WT-WTMAX)
      CALL Gmonit(   0,IDA+18,  WTovr,0D0,0D0)
      WTneg = MIN(WT,0D0)
      CALL Gmonit(   0,IDA+19,  WTneg,0D0,0D0)
      CALL GF1( 9001, WT,1D0)
!--------------------------------------------------------------------
!  Optional Rejection according to principal weight  
      IF(KeyWgt .EQ. 0) THEN  
!-- Constant weihgt events WT=1 
            RN = DRVEC(2)
            CALL Gmonit(  0, IDGEN, SIG0NB*WTMAX, WTMAX,1D0)
            IF(WT.LT.RN*WTmax) GOTO 200
            WTMOD=1.D0 
            WTCRU1=1D0  
            WTCRU2=1D0
!-- Precaution measure
            DO i=1,300
               WTSET(i)=0
            ENDDO
         ELSE 
!-- Weighted events  
            WTCRU1= WT0*WT1*WT2
            WTCRU2=     WT3*WT4
            WTMOD  = WTCRU1*WTCRU2*WT5
            CALL Gmonit(  0, IDGEN, SIG0NB, WTMAX,1D0)
         ENDIF
!
!--------------------------------------------------------------------
!  Merge photons/fermions into single common block
      CALL Mergik
!
!     ===========
      ELSE 
!     ===========
      npar(10)= nevgen  
      npar(20)= nevgen   
      CALL Gmonit(   1,ida+20,awtot,dwtot,wwmx )
      xsmc   = sig0  *averwt 
      xsmcnb = sig0nb*averwt 
      erel   = errela 
      ermc   = xsmcnb*erel 
      xpar(10) = xsmcnb
      xpar(11) = erel
      xpar(12) = xsmc 
      IF(KeyWgt .EQ. 0) THEN 
!-- WT=1  events, normal option... 
         xpar(20)=xsmcnb 
         xpar(21)=erel 
         xpar(22)=xsmc
      ELSE
!-- Weighted events, additional information on x-sections 
         xpar(20)= sig0nb
         xpar(21)= 0d0  
         xpar(22)= sig0
      ENDIF    

! printout only for MODE=2
      IF(MODE.EQ.1) RETURN    
!
      CALL Gprint(9001) 
!
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '  BHLUM4:        WINDOW A        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXL1I) NEVGEN,     ' Accepted total    ','NEVGEN','A1'
      WRITE(NOUT,BXL1I) IEVENT,     ' Raw prior reject. ','IEVENT','A2'
      WRITE(NOUT,BXL2G) XSMCNB,ERMC,' Xsec M.C. [nb]    ','XSECMC','A3'
      WRITE(NOUT,BXL1F) EREL,       ' relat. error      ','ERELMC','A4'
      WRITE(NOUT,BXL2F) AWTOT,DWTOT,' weight  M.C.      ','AWT   ','A5'
      WRITE(NOUT,BXL1I) NEVNEG,     ' WT<0              ','NEVNEG','A6'
      WRITE(NOUT,BXL1I) NEVOVE,     ' WT>WTMAX          ','NEVOVE','A7'
      WRITE(NOUT,BXL1F) WTMAX ,     ' Maximum WT        ','WWMX  ','A8'
      WRITE(NOUT,BXCLO)  
! Print additional infos                  
!------------------------------------------------------------       
      CALL Gmonit(   1,IDA+1, AWT1 ,DWT1 ,DUMM3)                  
      CALL Gmonit(   1,IDA+2, AWT2 ,DWT2 ,DUMM3)         
      CALL Gmonit(   1,IDA+3, AWT3 ,DWT3 ,DUMM3)         
      CALL Gmonit(   1,IDA+4, AWT4 ,DWT4 ,DUMM3) 
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '  BHLUM4:        WINDOW B        '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXL2F) AWT1,DWT1  ,'WT1*WT2*T/TP*T/TQ  ','      ','B1'
      WRITE(NOUT,BXL2F) AWT3,DWT3  ,'WT3 from KINO4     ','      ','B2'
      WRITE(NOUT,BXL2F) AWT4,DWT4  ,'YFS formfac        ','WT    ','B4'
      WRITE(NOUT,BXL2F) AWTOT,DWTOT,'TOTAL              ','      ','B5'
      CALL Gmonit(   1,IDA+18, AWT18 ,RWT18 ,DUMM3)     
      XWT18 = AWT18/AWTOT
      DWT18 = XWT18*RWT18
      WRITE(NOUT,BXL2F) XWT18,DWT18,'xsec/xtot: WT>WTMAX','WT    ','B6'
      CALL Gmonit(   1,IDA+19, AWT19 ,RWT19 ,DUMM3)     
      XWT19 = AWT19/AWTOT
      DWT19 = XWT19*RWT19
      WRITE(NOUT,BXL2F) XWT19,DWT19,'xsec/xtot: WT<0    ','WT    ','B7'
      WRITE(NOUT,BXCLO)  
! ---------------------------------------------------------------
      CALL Gmonit( 1,IDC+1,AWT1,DWT1,DUMM3)                            
      CALL Gmonit( 1,IDC+2,AWT2,DWT2,DUMM3)          
      CALL Gmonit( 1,IDC+3,AWT3,DWT3,DUMM3)          
      CALL Gmonit( 1,IDC+4,AWT4,DWT4,DUMM3)          
      CALL Gmonit( 1,IDC+5,AWT5,DWT5,DUMM3)          
      CALL Gmonit( 1,IDC+6,AWT6,DWT6,DUMM3)          
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '            WINDOW C             '
      WRITE(NOUT,BXTXT) 'Built-in average control weights.'
      WRITE(NOUT,BXTXT) 'Should equal one +- statist. err.'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXL2F) AWT1,DWT1  ,'<WCTA1>            ','      ','C1'
      WRITE(NOUT,BXL2F) AWT2,DWT2  ,'<WCTA2>            ','      ','C2'
      WRITE(NOUT,BXL2F) AWT3,DWT3  ,'<WCTA1*WCTA2>      ','      ','C3'
      WRITE(NOUT,BXL2F) AWT4,DWT4  ,'<WCTB1>            ','      ','C4'
      WRITE(NOUT,BXL2F) AWT5,DWT5  ,'<WCTB2>            ','      ','C5'
      WRITE(NOUT,BXL2F) AWT6,DWT6  ,'<WCTB1*WCTB2>      ','      ','C6'
      WRITE(NOUT,BXCLO)  
      ENDIF 
!     =====
      END  

      SUBROUTINE FILBH2(XPAR,NPAR) 
!     **************************** 
      IMPLICIT REAL*8 (A-H,O-Z)   
      DIMENSION  XPAR(*), NPAR(*) 
      PARAMETER( PI = 3.1415926535897932D0, ALFINV = 137.03604D0)
      CHARACTER*80   BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G 
      PARAMETER(
     $BXOPE =  '(//1X,15(5H=====)    )',
     $BXTXT =  '(1X,1H=,                  A48,25X,    1H=)',
     $BXL1I =  '(1X,1H=,I17,                 16X, A20,A12,A7, 1X,1H=)',
     $BXL1F =  '(1X,1H=,F17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2F =  '(1X,1H=,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXL1G =  '(1X,1H=,G17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2G =  '(1X,1H=,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXCLO =  '(1X,15(5H=====)/   )'    )   
      COMMON / TRANSR / TRAN,TRMIN,TRMAX
      COMMON / BHPAR3 / KeyRad,KeyOpt
      COMMON / BHPAR2 / CMSENE,AMEL  
      COMMON / BHPAR1 / DEL,EPSCM,THMIN,XMIVIS 
      COMMON / BHPARZ / AMAZ,GAMMZ,SINW2,GV,GA
! Communicates with VarRan
      COMMON / RANPAR / KeyRnd
      COMMON / INOUT  / NINP,NOUT
      real *4 enq,st2,der,errder,deg,errdeg
      SAVE   / INOUT  /, / RANPAR /
      SAVE   / TRANSR /, / BHPAR3 /, / BHPAR2 / , / BHPAR1 /, / BHPARZ /
!--
      CMSENE=XPAR(1)
      TRMIN =XPAR(2)
      TRMAX =XPAR(3)
      EPSCM =XPAR(4)
      KeyOpt=NPAR(1)
      KeyRnd = MOD(KeyOpt,10)   
      KeyWgt = MOD(KeyOpt,100)/10   
      KeyRem = MOD(KeyOpt,1000)/100   
      KeyRad=NPAR(2)
      KeyPia = MOD(KeyRad,10)
      KeyMod = MOD(KeyRad,100)/10 
      KeyUpd = MOD(KeyRad,1000)/100 
      KeyZet = MOD(KeyRad,10000)/1000 
!--
      AMEL  =  0.5111D-3
      SVAR  = CMSENE**2
      XIMIN = TRMIN/SVAR
      XIMAX = TRMAX/SVAR
!&&&&&&& varius choices adjusted manualy
      DEL   = EPSCM*0.01D0  ! This one seems to be optimal
!&    DEL   = EPSCM*0.1D0
!&    DEL   = EPSCM*XIMIN
!&&&&&&&
      IF(TRMAX.GT.SVAR) TRMAX=SVAR
! Inputs for Z correction   
      AMAZ  = 91.187d0
      GAMMZ =  2.490d0
      SINW2 = 0.2319d0
      GA = -1/(4*DSQRT(SINW2*(1-SINW2)))
      GV = GA*(1-4*SINW2)
! Test of vacuum correction 
      CALL vacpol(KeyPia,-trmin,SINW2,RePi1,dRePi1)
      CALL vacpol(KeyPia,-trmax,SINW2,RePi2,dRePi2)

!--
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '  BHLUM4: INPUT PARAMETRES       '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXL1I) KeyOpt,     ' OPTIONS   switch  ','KeyOpt','  '
      WRITE(NOUT,BXL1I) KeyRnd,     ' rand. numb. switch','KeyRnd','  '
      WRITE(NOUT,BXL1I) KeyWgt,     ' weighting switch  ','KeyWgt','  '
      WRITE(NOUT,BXL1I) KeyRem,     ' photon removal  sw','KeyRem','  '
      WRITE(NOUT,BXL1I) KeyRad,     ' RADIATION switch  ','KeyRad','  '
      WRITE(NOUT,BXL1I) KeyPia,     ' vac_pol   switch  ','KeyPia','  '
      WRITE(NOUT,BXL1I) KeyMod,     ' QED mat. elm. type','KeyMod','  '
      WRITE(NOUT,BXL1I) KeyUpd,     ' Test switch, def=0','KeyUpd','  '
      WRITE(NOUT,BXL1I) KeyZet,     ' Z contribution    ','KeyZet','  '
      WRITE(NOUT,BXL1F) CMSENE,     ' CMS energy  [GeV] ','CMSENE','  '
      WRITE(NOUT,BXL1G) SVAR  ,     ' CMSENE^2  [GeV^2] ','SVAR  ','  '
      WRITE(NOUT,BXL1G) TRMIN ,     ' trasf_min [GeV^2] ','TRMIN ','  '
      WRITE(NOUT,BXL1G) TRMAX ,     ' trasf_max [GeV^2] ','TRMAX ','  '
      WRITE(NOUT,BXL1G) XIMIN ,     ' xi_min=TRMIN/SVAR ','XIMIN ','  '
      WRITE(NOUT,BXL1G) XIMAX ,     ' xi_max=TRMAX/SVAR ','XIMAX ','  '
      THMIN  = ACOS(1-2*XIMIN)*1000D0
      THMAX  = ACOS(1-2*XIMAX)*1000D0
      WRITE(NOUT,BXL1F) THMIN ,     ' theta_min  [mrad] ','THMIN ','  '
      WRITE(NOUT,BXL1F) THMAX ,     ' theta_max  [mrad] ','THMAX ','  '
      THMIN  = ACOS(1-2*XIMIN)*180/PI
      THMAX  = ACOS(1-2*XIMAX)*180/PI
      WRITE(NOUT,BXL1F) THMIN ,     ' theta_min   [deg] ','THMIN ','  '
      WRITE(NOUT,BXL1F) THMAX ,     ' theta_max   [deg] ','THMAX ','  '
      WRITE(NOUT,BXL1G) EPSCM ,     ' eps_CM infr. cut  ','EPSCM ','  '
      WRITE(NOUT,BXL1G) DEL   ,     ' delta  infr. cut  ','DEL   ','  '
      WRITE(NOUT,BXL1F) REPI1 ,     ' RePi(transf_min)  ','REPI1 ','  '
      WRITE(NOUT,BXL1F) dREPI1,     ' error             ','dREPI1','  '
      WRITE(NOUT,BXL1F) REPI2 ,     ' RePi(transf_max)  ','REPI2 ','  '
      WRITE(NOUT,BXL1F) dREPI2 ,     'error             ','dREPI2','  '
      WRITE(NOUT,BXL1F) AMAZ  ,     ' Z-mass GeV        ','AMAZ  ','  '
      WRITE(NOUT,BXL1F) GAMMZ ,     ' Z-width GeV       ','GAMMZ ','  '
      WRITE(NOUT,BXL1F) SINW2 ,     ' weak mixing angle ','SINW2 ','  '
      WRITE(NOUT,BXCLO)  
      END 

      SUBROUTINE REMPHO(EMIN,NPHOT,PHOT,ALF,BET,WTM,MK)
!     *************************************************
      IMPLICIT REAL*8(A-H,O-Z)       
      DIMENSION PHOT(100,4),ALF(50),BET(50),WTM(50),MK(50)
!
      if(nphot.eq.0) return
      NPH=NPHOT 
      DO 100 J=NPHOT,1,-1 
      IF(PHOT(J,4).LT.EMIN) THEN 
         DO 60 I=J+1,NPH  
         ALF(I-1)=ALF(I)
         BET(I-1)=BET(I)
         WTM(I-1)=WTM(I)
         MK( I-1)=MK( I)
         DO 60 K=1,4      
   60    PHOT(I-1,K)=PHOT(I,K) 
         NPH=NPH-1 
      ENDIF   
  100 CONTINUE
!.....................................................
! Correction by Alex Read and Tiziano Camporesi DELPHI
! Date: Fri, 25 Nov 94 17:50:37 WET
!      Code added by ALR 22.11.94 to fix problem with
!      photon handling. Have to erase the discarded
!      photons or they cause occasional problems in
!      MERGIK when merging the PHOT1 and PHOT2 arrays
!      (REMPHO operates on these here).
!
      DO J=NPH+1,NPHOT
         DO K=1,4
           PHOT(J,K) = 0.D0
         ENDDO
      ENDDO
!.....................................................
      NPHOT=NPH
      END

      SUBROUTINE PIATEK(CMSENE,TRMAX,AMEL,EMIN,DELTA,
     $         NPHOT,P1,P2,PHOT,PHSU,WMAT,WTAL,WTMRE,WCTR1,WCTR2)
!     ***************************************************************
! Input:
!        CMSENE         CMS energy
!        TRMAX          maximum transfer (dummy)
!        AMEL           electron mass (GeV)
!        EMIN           CMS minimum photon energy (GeV)
!        DELTA          MC minimum photon energy (dimensionless)
!        NPHOT          photon number
!        P1,P2(4)       fermion momenta
!        PHOT(100,4)    photon four-momenta
!        PHSU(50)       sum of photon four-momenta
!        WMAT(50)       mass weights from MLTIBR
! Output:
!        WTAL       mass weight for all photons
!        WTMRE      In the case of removal the new mass weight
!        WCTR1      Control weight for delta-->epsilon rejection
!        WCTR2      control weight for photons below EMIN removal
!     ************************
      IMPLICIT REAL*8(A-H,O-Z)  
      PARAMETER(PI=3.1415926535897932D0,ALFINV=137.03604D0)
      PARAMETER(ALF1=1/ALFINV/PI )
      DIMENSION P1(4),P2(4),WMAT(50),PHOT(100,4),PHSU(4)
      DATA ICONT /0/
      ICONT=ICONT+1
! Calculate mass weight for all photos and separately for
! photons with energy below/above EMIN
      EPSCM = 2*EMIN/CMSENE
      WTAL = 1D0   
      WTM1 = 1D0
      WTM2 = 1D0
      WTEPSP=1D0
      DO 100 I=1,NPHOT 
      WTAL = WTAL*WMAT(I)
      IF(WTAL.LT.1D-15) WTAL =0D0
      IF(PHOT(I,4).LT.EMIN) THEN
        WTM1 = WTM1*WMAT(I)
        IF(WTM1.LT.1D-15) WTM1=0D0
        WTEPSP = 0D0
      ELSE
        WTM2=WTM2*WMAT(I)
        IF(WTM2.LT.1D-15) WTM2=0D0
      ENDIF
  100 CONTINUE   
! Control weight for delta-->epsilon  'REJECTION'
      DELT1 = DELTA*BCUD(P1,P2,PHSU)
      CALL WFORM(TRMAX,P1,P2,AMEL,DELT1,EMIN,PDYFS)
      WCTR1 = WTEPSP*PDYFS   
! control weight for photons ENE<EMIN  'REMOVAL'
      TRANP = 2D0*(P1(4)*P2(4)-P1(3)*P2(3)-P1(2)*P2(2)-P1(1)*P2(1))
      EPS1  =  SQRT(EMIN**2/P1(4)/P2(4))     
      DELB2 = -2*ALF1*(DLOG(TRMAX/TRANP)+1) *DLOG(EPS1/DELT1)
      WCTR2 = WTM1*EXP(-DELB2)   
! In the case of removal the new mass weight is this
      WTMRE = WTM2*EXP(DELB2)
      END

      SUBROUTINE WFORM(TRMAX,Q1,Q2,AMF,DELTA,EMIN,DYFS)
!     *************************************************  
! For tests only.     
! Yennie-Frautschi-Suura Formfactors for the single fermion pair 
! This is for crude distribition before mass wights
! The triangle effect included (pi**2/6)
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI=3.1415926535897932D0,ALFINV=137.03604D0)
      PARAMETER(ALF1=1/ALFINV/PI )
      DIMENSION Q1(4),Q2(4)       
! ...Momenta q1,q2 should be in CMS
      Q1Q2  = Q1(4)*Q2(4)-Q1(3)*Q2(3)-Q1(2)*Q2(2)-Q1(1)*Q2(1) 
      E1 = Q1(4) 
      E2 = Q2(4)
      BETF2 = 2*ALF1* DLOG(TRMAX /AMF**2) 
      DELB  = BETF2*DLOG(EMIN/SQRT(E1*E2)/DELTA)
      EP    = E1+E2
      EM    = E1-E2  
      DL    = SQRT( 2*Q1Q2 +EM**2 )     
! Note that approximately REMN= +(1./6.)*PI**2 for t-channel
      REMN  = PI**2/2 
     $        -0.50*DLOG(E1/E2)**2  
     $        -0.25*DLOG((DL+EM)**2/(4*E1*E2))**2 
     $        -0.25*DLOG((DL-EM)**2/(4*E1*E2))**2    
     $        - DILOGY((DL+EP)/(DL+EM)) -DILOGY((DL-EP)/(DL-EM))
     $        - DILOGY((DL-EP)/(DL+EM)) -DILOGY((DL+EP)/(DL-EM)) 
! This (alf/pi)*pi**2/6 is related to replacing (y+z)>epsilon
! by max(y,z)>epsilon.   (Rejection delta=> epsilon over-estimated)
      TRIANG = -PI**2/6D0 
      DYFS   = EXP( DELB +ALF1*REMN +ALF1*TRIANG) 
      END

      FUNCTION BCUD(P1,P2,SF)
!     ************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P1(4),P2(4),SF(4)
      XPP  = P2(4)*P1(4)-P2(3)*P1(3)-P2(2)*P1(2)-P2(1)*P1(1)
      XPR  = P1(4)*(P2(4)+SF(4)) - P1(3)*(P2(3)+SF(3))
     $     - P1(2)*(P2(2)+SF(2)) - P1(1)*(P2(1)+SF(1))
      BCUD= XPR/XPP
      END

      SUBROUTINE MLTIBR(TRAN,TRMAX,AMEL,DEL,
     $      NPH,PHOT,PHSU,ALF1,BET1,TRANP,AMSP,MK,WT1,WTM)   
!     ****************************************************   
! This provides momenta of photons in a fermion proper frame 
! Input : TRAN    = principal t-channel transfer     (GEV**2)
!         TRMAX   = max. transf. (>0) for angular phot. dist. [GEV**2]
!         AMEL    = electron energy         (GEV)
!         DEL     = low energy photon limit   (dimensionless)
! Output: NPH     = photon multiplicity
!         PHOT    = list of photon four-momenta
!         PHSU    = sum of photon momenta
!         ALF1,BET1   = Sudakov variables
!         TRANP   = (P2-P1)**2
!         AMSP    = (P2+PHSU)**2  
!         MK      = marked photons
!         WT1     = TRANP/TRAN is Jacobian, =0 outside ph.sp.
!         WTM     = list of mass weights
!     ************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0)
      PARAMETER (NMAX=20)
      DIMENSION PHOT(100,4),PHSU(4),PH(4)
      DIMENSION MK(50),WTM(50)   
      DIMENSION ALF(50),BET(50),ALF1(50),BET1(50),Y(50),Z(50)
      DIMENSION RR(100),P2(4)
      DOUBLE PRECISION DRVEC(100)
      DATA ICONT /0/ 
      ICONT=ICONT+1 
    
      DELS  = AMEL**2/TRMAX
      BILGS = LOG(TRMAX/AMEL**2)
      DELL  = LOG(DEL) 
      WT1    = 1D0
      DO  11 I=1,50 
      DO  11 K=1,4 
   11 PHOT(I,K)=0D0
      AVERG=2D0/(PI*ALFINV)*BILGS*LOG(1D0/DEL)
      CALL POISSG(AVERG,NMAX,NPH,RR) 
      IF(NPH.GT.NMAX) GOTO 900  
! No photon        
      DO 45 K=1,4
   45 PHSU(K)=0D0
      IF(NPH.EQ.0) THEN 
        TRANP=TRAN  
      ELSE  
! One or more photons 
   50   CALL VarRan(DRVEC,NPH) 
        BSUM=0D0 
        DO 80 I=1,NPH  
! We define R=LOG(MAX(YGR,ZET)) 
        R=DELL*RR(I) 
! Photons close to lower infrared boundry are marked for tests 
        MK(I)=0 
        IF(EXP(R).LT.DEL*3D0) MK(I)=1  
        T= 2D0*DRVEC(I)
        IF(T.GT.1D0) THEN  
           YGR=R  
           ZET=R-(T-1D0)*BILGS  
        ELSE  
           ZET=R     
           YGR=R-T*BILGS 
        ENDIF  
        YGR=EXP(YGR) 
        ZET=EXP(ZET) 
        Y(I)=YGR 
        Z(I)=ZET 
! Define ALPHA and BETA (prim)
        ALF1(I)=YGR-ZET*DELS
        BET1(I)=ZET-YGR*DELS 
        IF(ALF1(I).LE.0D0.OR.BET1(I).LE.0D0) GOTO 50
   80   BSUM=BSUM+BET1(I)   
        IF(BSUM.GT.1D0) GOTO 800
! Rescale ALPHA and BETA        
        CALL VarRan(DRVEC,NPH) 
        DO 90 I=1,NPH 
        ALF(I)=ALF1(I)/(1D0-BSUM) 
   90   BET(I)=BET1(I)/(1D0-BSUM)  
! Define photon four momenta in SQRT(TRANP)/2 units 
        DO 100 I=1,NPH 
        PHOT(I,4)= ALF(I)+BET(I) 
        PHOT(I,3)=-ALF(I)+BET(I)
        R1 = DRVEC(I)
        PHI=2D0*PI*R1
        PTRAN=2D0*DSQRT(DABS(ALF(I)*BET(I)))
        PHOT(I,1)=PTRAN*COS(PHI) 
        PHOT(I,2)=PTRAN*SIN(PHI) 
        DO 100 K=1,4 
  100   PHSU(K)=PHSU(K)+PHOT(I,K)
! Define factor for rescaling photon momenta    
        XMK2=PHSU(4)**2-PHSU(3)**2-PHSU(2)**2-PHSU(1)**2 
        YY2=1D0/(1D0+PHSU(3)-.25D0*XMK2)      
! YY2 negative when outside phase space (abs(t)>abs(s))
        IF(YY2.LE.0D0) GOTO 900 
        TRANP=TRAN*YY2      
        ENER =SQRT(TRANP)/2D0        
! RESCALE ALL PHOTON MOMENTA         
        DO 120 K=1,4        
        PHSU(K)=PHSU(K)*ENER         
        DO 120 I=1,NPH      
  120   PHOT(I,K)=PHOT(I,K)*ENER     
! This rotation makes PHSU(2)=0 
! (we get rid her of dummy angle, see 'poprawka' in the notes)
        PSIT=ANGFI(PHSU(1),PHSU(2))  
        CALL ROTOD3(-PSIT, PHSU,PHSU)         
        DO 140 I=1,NPH      
        DO 135 K=1,4        
  135   PH(K)=PHOT(I,K)     
        CALL ROTOD3(-PSIT, PH,PH)    
        DO 140 K=1,4        
  140   PHOT(I,K)=PH(K)     
      ENDIF         
!+++      IF(TRANP.EQ.0D0) GO TO 900      
      IF(TRANP.LE. 4*AMEL**2) GO TO 900      
      BETEL=SQRT(1D0-4D0*AMEL**2/(TRANP+4D0*AMEL**2))          
      P2(3)=SQRT(TRANP)/2D0          
      P2(4)=SQRT(P2(3)*P2(3)+AMEL**2)         
      P2(2)=0D0     
      P2(1)=0D0     
      AMSP=(P2(4)+PHSU(4))**2-(P2(3)+PHSU(3))**2      
     $    -(P2(2)+PHSU(2))**2-(P2(1)+PHSU(1))**2      
! And the weight finally    
      WT1 = TRANP/TRAN  
      DELT=AMEL**2/TRANP    
      DO 200 I=1,NPH        
! Generated distribution         
! here some numerical regularization  
      DIST0 = 1D0/((ALF1(I)+DELS*BET1(I))*(BET1(I)+DELS*ALF1(I)))
      YGR=ALF1(I)+DELT*BET1(I)
      ZET=BET1(I)+DELT*ALF1(I)       
! Desired distribution = soft factor
      DIST1 = ALF1(I)*BET1(I)/(YGR*ZET)**2
      WTM(I)= DIST1/DIST0     
  200 CONTINUE
      RETURN 
! Event outside phase space
! Note that distinction is made (TRANP=-2,-1) to facilitate tests 
! event dropped due to: sum(beta) > 1
 800  CONTINUE
      TRANP = -1D0
      WT1   =  0D0
      RETURN
! event dropped due to: tranp < m^2, or earlier because YY2 < 0
 900  CONTINUE
      TRANP = -2D0
      WT1   =  0D0     
      END  
  
      SUBROUTINE POISSG(AVERG,NMAX,MULT ,RR)      
!     **************************************
! DIFFERS FROM THAT IN EXPAND DEC. 87         
! THIS GENERATES PHOTON MULTIPLICITY ACCORDING TO POISSON DISTRIBUTION
! INPUT:  AVERG = AVERAGE MULTIPLICITY        
!         NMAX  = MAXIMUM MULTIPLICITY        
! OUTPUT: MULT  = GENERATED MULTIPLICITY      
!         RR(1:100) LIST OF ORDERED UNIFORM RANDOM NUMBERS,    
!         A BYPRODUCT RESULT, TO BE EVENTUALLY USED FOR SOME FURTHER  
!         PURPOSE (I.E.  GENERATION OF PHOTON ENERGIES).       
!     **************************************      
      IMPLICIT REAL*8(A-H,O-Z)       
      REAL*8 RR(100)        
      DOUBLE PRECISION DRVEC(100)
      SAVE NFAIL
      DATA NFAIL/0/         
   50 NN=0    
      IF(NMAX.GT.100) GOTO 900       
      CALL VarRan(DRVEC,NMAX)
      SUM=0D0         
      DO 100 IT=1,NMAX      
      RN = DRVEC(IT)
      Y= LOG(RN)    
      SUM=SUM+Y     
      NN=NN+1       
      IF(SUM.LT.-AVERG) GOTO 130     
      RR(NN)=SUM/(-AVERG)   
  100 CONTINUE      
      NFAIL=NFAIL+1         
      IF(NFAIL.GT.100) GOTO 900      
      GOTO 50       
  130 MULT =NN-1    
      RETURN        
  900 WRITE(6,*) ' POISSG: TO SMALL OR TO BIG NMAX',NMAX  
      STOP 
      END  
  
      SUBROUTINE KINO4(SVAR,TRAN,AMEL,AMSP,AMSQ,WTKK)   
!     ************************************************
! Kinematics, cnstruction of momenta in CMS   
!     ************************************************
      IMPLICIT REAL*8(A-H,O-Z)      
      PARAMETER( PI =3.1415926535897932D0)
      COMMON / MOMS1  / TRANP,P1(4),P2(4),PHOT1(100,4),PHSU1(4),NPHOT1
      COMMON / MOMS2  / TRANQ,Q1(4),Q2(4),PHOT2(100,4),PHSU2(4),NPHOT2
      COMMON / MOMZ1  / AL1(50),BE1(50),WTM1(50),MK1(50)
      COMMON / MOMZ2  / AL2(50),BE2(50),WTM2(50),MK2(50)
      COMMON / PSIPHI / TH1,EXT1,EXB1,PSI1,EXW1,
     $                  TH2,EXT2,EXB2,PSI2,EXW2,EXE2,FIF,PHI
      SAVE   / MOMS1  /, / MOMS2  /, / MOMZ1  /, / MOMZ2  /, /PSIPHI/
      DOUBLE PRECISION DRVEC(100)
      REAL*8 PH(4)  
      DIMENSION QCM(4)
        
      BTEL=DSQRT(1D0-4D0*AMEL**2/SVAR)        
      WTKK=1D0
! Three azimuthal angles        
      CALL VarRan(DRVEC,3)      
      PSI1= 2D0*PI*DRVEC(1)
      PSI2= 2D0*PI*DRVEC(2)
      PHI = 2D0*PI*DRVEC(3)
! Upper vertex: transf. from P2-P1 proper frame       
      CALL KLIPER(TRANP,AMEL,PHSU1,P2,TH1,EXT1,EXB1)  
! Lower vertex: trascf. from q2-q1 proper frame       
      CALL KLIPER(TRANQ,AMEL,PHSU2,Q2,TH2,EXT2,EXB2)  
! Define P1, Q1  in central QMS      
      P1(3)= -(TRAN+AMSP-AMEL**2)/SQRT(TRAN)/2D0      
      Q1(3)=  (TRAN+AMSQ-AMEL**2)/SQRT(TRAN)/2D0      
      RPQK=(Q1(3)+P1(3))/DSQRT(SVAR) 
! Correcting for electron mass       
!     PX2=SVAR*(SVAR+4D0*P1(3)*Q1(3))/((Q1(3)+P1(3))**2+SVAR)/4D0
!     PX2=PX2-AMEL**2       
      GPQK= P1(3)-Q1(3)     
      PX2=(BTEL**2*SVAR*(1D0+RPQK*RPQK)-GPQK*GPQK)/(1D0+RPQK*RPQK)/4D0
      IF(PX2.LE.0D0)  GOTO 900       
      PX=SQRT(PX2)  
      P1(2)=  0D0   
      Q1(2)=  0D0   
      P1(1)=  -PX   
      Q1(1)=   PX   
      P1(4)=  SQRT(P1(1)**2+P1(2)**2+P1(3)**2+AMEL**2)         
      Q1(4)=  SQRT(Q1(1)**2+Q1(2)**2+Q1(3)**2+AMEL**2)         
! Correcting for electron mass       
!     BETP = SQRT(1D0-(AMEL/P1(4))**2)        
!     BETQ = SQRT(1D0-(AMEL/Q1(4))**2)        
!     DO 7 K=1,3    
!     P1(K)=BETP* P1(K)     
!   7 Q1(K)=BETQ* Q1(K)     
      EXW1=SQRT((P1(4)+P1(1))/(P1(4)-P1(1)))  
      EXW2=SQRT((Q1(4)+Q1(1))/(Q1(4)-Q1(1)))  
! Construct momentum transfer Q in CMS        
      QCM(4)=(AMSP-AMSQ)/SQRT(SVAR)/2D0       
      QMOD=SQRT(TRAN+QCM(4)**2)      
      QCM(3)=(-TRAN-AMSP/2D0-AMSQ/2D0+AMEL**2)/SQRT(SVAR-4D0*AMEL**2)
      QCM(2)=0D0    
      QCM(1)=SQRT(QMOD**2-QCM(3)**2) 
      FIF =ANGFI(QCM(3),QCM(1))      
      EXE2=SQRT((QMOD+QCM(4))/(QMOD-QCM(4)))  
  
! Final set of transformations from QMSP and QMSQ to CMS       
! First branch, tranformed are P2, PHSU1, PHOT1
      CALL  PTRAL(TH1,EXT1,EXB1,PSI1,EXW1,EXE2,FIF,PHI,P2)     
      IF(NPHOT1.NE.0) THEN  
       CALL PTRAL(TH1,EXT1,EXB1,PSI1,EXW1,EXE2,FIF,PHI,PHSU1)  
       DO 20 I=1,NPHOT1     
       DO 15 K=1,4  
   15  PH(K)=PHOT1(I,K)     
       CALL PTRAL(TH1,EXT1,EXB1,PSI1,EXW1,EXE2,FIF,PHI,PH)     
       DO 16 K=1,4  
   16  PHOT1(I,K)=PH(K)     
   20  CONTINUE     
      ENDIF         
! Second branch, tranformed are Q2, PHSU2, PHOT2
      CALL  QTRAL(TH2,EXT2,EXB2,PSI2,EXW2,EXE2,FIF,PHI,Q2)     
      IF(NPHOT2.NE.0) THEN  
       CALL QTRAL(TH2,EXT2,EXB2,PSI2,EXW2,EXE2,FIF,PHI,PHSU2)  
       DO 30 I=1,NPHOT2     
       DO 25 K=1,4  
   25  PH(K)=PHOT2(I,K)     
       CALL QTRAL(TH2,EXT2,EXB2,PSI2,EXW2,EXE2,FIF,PHI,PH)     
       DO 26 K=1,4  
   26  PHOT2(I,K)=PH(K)     
   30  CONTINUE     
      ENDIF         
! Finally, beams P1 and Q1   
      CALL BOSTD3(EXE2,P1,P1)        
      CALL ROTOD2( FIF,P1,P1)        
      CALL BOSTD3(EXE2,Q1,Q1)        
      CALL ROTOD2( FIF,Q1,Q1)        
      RETURN        
! Event outside phase space          
  900 WTKK=0D0      
      END  
  
      SUBROUTINE PTRAL(TH,EXT,EXB,PSI,EXW,EXE,FIF,PHI,P)       
!     **************************************************       
      IMPLICIT REAL*8(A-H,O-Z)       
      REAL*8 P(4)   
      CALL ROTOD2( -TH, P, P)        
      CALL BOSTD3( EXT, P, P)        
      CALL BOSTD1( EXB, P, P)        
      CALL ROTOD3( PSI, P, P)        
      CALL BOSTD1( EXW, P, P)        
      CALL BOSTD3( EXE, P, P)        
      CALL ROTOD2( FIF, P, P)        
      CALL ROTOD3( PHI, P, P)        
      END  
  
      SUBROUTINE QTRAL(TH,EXT,EXB,PSI,EXW,EXE,FIF,PHI,P)       
!     **************************************************       
      IMPLICIT REAL*8(A-H,O-Z)       
      PARAMETER( PI =3.1415926535897932D0) 
      REAL*8 P(4)   
      CALL ROTOD2( -TH, P, P)        
      CALL BOSTD3( EXT, P, P)        
      CALL BOSTD1( EXB, P, P)        
      CALL ROTOD3( PSI, P, P)        
      CALL ROTOD2(  PI, P, P)        
      CALL BOSTD1( EXW, P, P)        
      CALL BOSTD3( EXE, P, P)        
      CALL ROTOD2( FIF, P, P)        
      CALL ROTOD3( PHI, P, P)        
      END  
  
  
      SUBROUTINE KLIPER(TRANP,AMEL,PHSUM,P2,TH,EXT,EXB)
!     **************************************************
! Deals with Lorentz transf. from QQ1 to QQ frame
! where QQ1=P2-P1, QQ=P2+PHSUM-P1, TRANP=QQ1**2, P1**2=P2**2=AMEL**2 
! Input: TRANP,AMEL,PHSUM
! Output: P2,TH,EXT,EXB,PHSUM
! Here, TH, EXT, EXB are transformation params.
!     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)       
      REAL*8 PHSUM(4),P2(4)          
      REAL*8 P1(4),QQ1(4)   
  
      BETEL=SQRT(1D0-4D0*AMEL**2/(TRANP+4D0*AMEL**2))          
! No photon         
      IF(PHSUM(4).EQ.0D0) THEN       
        P2(3)= SQRT(TRANP)/2D0       
        P2(4)= SQRT(P2(3)*P2(3)+AMEL**2)      
        P2(2)=0D0   
        P2(1)=0D0           
        TH =0D0         
        EXT=1D0         
        EXB=1D0         
      ELSE  
! One photon or more    
        ENER1=SQRT(TRANP)/2D0    
   
        P1(1)=0D0       
        P1(2)=0D0       
        P1(3)=-ENER1    
        P1(4)= SQRT(P1(3)*P1(3)+AMEL**2)  
   
        P2(1)=0D0       
        P2(2)=0D0       
        P2(3)= ENER1    
        P2(4)= SQRT(P2(3)*P2(3)+AMEL**2)  
   
        DO 33 I=1,4     
  33    QQ1(I)=P2(I)+PHSUM(I)-P1(I)       
   
! Rotation 2 puts QQ1 paralel to axis 3   
! Note that  PHSUM(2)=0 is already assured in MLTIBR!
        TH  =ANGFI(QQ1(3),QQ1(1))         
        CALL ROTOD2(-TH ,QQ1,QQ1)         
        CALL ROTOD2(-TH ,P1,P1)  
! Boost 3 puts QQ1(4)=0          
        EXT = SQRT((QQ1(3)-QQ1(4))/(QQ1(3)+QQ1(4)))   
        CALL BOSTD3( EXT ,QQ1,QQ1)        
        CALL BOSTD3( EXT , P1, P1)        
        EXB = SQRT((P1(4)-P1(1))/(P1(4)+P1(1)))       
!  Testing obsolete appendix
!  Boost 1 puts P1 antiparallel to axis 3  
!       CALL ROTOD2( -TH , P2, P2)
!       CALL BOSTD3( EXT , P2, P2)
!       CALL BOSTD1( EXB , P2, P2)
      ENDIF    
      END   

      SUBROUTINE MERGIK
!     *****************
! Transfer momenta and mark into proper commons        
! photons ordered according to cms energy 
! (the hardest in the first position)
!     ************************
      IMPLICIT REAL*8(A-H,O-Z)       
      COMMON / MOMS1  / TRANP,P1(4),P2(4),PHOT1(100,4),PHSU1(4),NPHOT1
      COMMON / MOMS2  / TRANQ,Q1(4),Q2(4),PHOT2(100,4),PHSU2(4),NPHOT2
      COMMON / MOMZ1  / AL1(50),BE1(50),WTM1(50),MK1(50)
      COMMON / MOMZ2  / AL2(50),BE2(50),WTM2(50),MK2(50)
      COMMON / MOMSET / PX1(4),QX1(4),PX2(4),QX2(4),PHOT(100,4),NPHOT
      COMMON / MARPKP / MARKP(100) 
      SAVE   / MOMS1/,/ MOMS2/,/ MOMZ1/,/ MOMZ2/,/ MOMSET/,/ MARPKP /
      NPHOT=NPHOT1+NPHOT2   
      I1=1 
      I2=1 
      DO 207 I=1,NPHOT      
      IF(PHOT1(I1,4).GT.PHOT2(I2,4)) THEN     
         DO 205 K=1,4       
  205    PHOT( I,K)=PHOT1(I1,K)      
         MARKP(I)  =  MK1(I1)        
         I1=I1+1    
      ELSE 
         DO 206 K=1,4       
  206    PHOT( I,K)=PHOT2(I2,K)      
         MARKP(I)  =  MK2(I2)        
         I2=I2+1    
      ENDIF         
  207 CONTINUE
      DO 300 K=1,4
      PX1(K)=P1(K)
      PX2(K)=P2(K)
      QX1(K)=Q1(K)
      QX2(K)=Q2(K)
  300 CONTINUE
      END      


      SUBROUTINE DUMPR(NUNIT,WORD,PP,QQ)    
!     **********************************
! 15 Jan 90 (SJ)
! prints twice dot-products of two four momentum PP and QQ
! more precisely:   2*PP.QQ  and  (PP+QQ).(PP+QQ)
!     ************************   
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8 WORD          
      REAL*8 PP(4),QQ(4)  
      DOT1=2*(PP(4)*QQ(4)-PP(3)*QQ(3)-PP(2)*QQ(2)-PP(1)*QQ(1))
      DOT2=(PP(4)+QQ(4))**2-(PP(3)+QQ(3))**2
     $    -(PP(2)+QQ(2))**2-(PP(1)+QQ(1))**2
      WRITE(NUNIT,'(1X,A8,5(1X,F20.10))') WORD,DOT1,DOT2        
      END  


      SUBROUTINE DUMPS(NOUT)     
!     **********************     
! THIS PRINTS OUT FOUR MOMENTA OF PHOTONS 
! ON OUTPUT UNIT NOUT
      IMPLICIT REAL*8(A-H,O-Z)   
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT  
      SAVE   / MOMSET /
      REAL*8 SUM(4)     
      WRITE(NOUT,*) '=====================DUMPS====================' 
      WRITE(NOUT,3100) ' P2',(P2(K),K=1,4)   
      WRITE(NOUT,3100) ' Q2',(Q2(K),K=1,4)   
      DO 100 I=1,NPHOT  
  100 WRITE(NOUT,3100) 'PHO',(PHOT(I,K),K=1,4)        
      DO 200 K=1,4      
  200 SUM(K)=P2(K)+Q2(K)         
      DO 210 I=1,NPHOT  
      DO 210 K=1,4      
  210 SUM(K)=SUM(K)+PHOT(I,K)    
      WRITE(NOUT,3100) 'SUM',(SUM(K),K=1,4)           
 3100 FORMAT(1X,A3,1X,5F18.13)   
      END
      
      SUBROUTINE INIT_LHE(m_EvtUnit, CMSENE)
      ! subrutine that inits the LHE file
      DOUBLE PRECISION CMSENE
      
      OPEN(m_EvtUnit, file = "hepevt.txt")
!1    WRITE(m_EvtUnit, '(a)')'<LesHouchesEvents version="1.0">'
!1    WRITE(m_EvtUnit, '(a)')'<!--'
!1    WRITE(m_EvtUnit, '(a)')'   File Created with BHLUMI'
!1    WRITE(m_EvtUnit, '(a)')'-->'
!1    WRITE(m_EvtUnit, '(a)')'<init>'
!1    WRITE(m_EvtUnit, '(A, es17.8, A, es17.8, A)' )'  11  -11  ', CMSENE/2.,'  ', CMSENE/2. , '  0  0  0  0  3  1'  
!1    WRITE(m_EvtUnit, '(a)')'  0.1  1.0e-06  1.000000e+00   9999'                                                   
!1    WRITE(m_EvtUnit, '(a)')'</init>'                                                                               
!1    WRITE(m_EvtUnit, '(a)')'  '                                                                                    

      END
      
      SUBROUTINE END_LHE(m_EvtUnit)
      ! subrutine that ends the LHE file   
!1    WRITE(m_EvtUnit, '(a)') '</LesHouchesEvents>'

      END
      
      SUBROUTINE WRITE_LHE(m_EvtUnit, CMSENE )
!     **********************
! THIS PRINTS OUT FOUR MOMENTA OF PHOTONS
! ON OUTPUT UNIT NOUT
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER nhep
      DOUBLE PRECISION AlphaQED, AlphaQCD, mass_e, CMSENE, MZstar
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      SAVE   / MOMSET /
      REAL*8 SUM(4)

      DO J=1,5
         SUM(J)=0.
      ENDDO
      
      mass_e=5.1099906999999998E-004

      
      DO I=1,NPHOT
         DO J=1,5
            SUM(J)= SUM(J)+ PHOT(I,J)
         ENDDO
      ENDDO

      MZstar=SQRT( (CMSENE-SUM(4))**2 - (SUM(1))**2 - (SUM(2))**2-(SUM(3))**2)

      AlphaQCD=0.1201789
      AlphaQED=137.03604D0
      

!1    WRITE(m_EvtUnit,'(a)')  '<event>'
      nhep=NPHOT+5
      
      WRITE(m_EvtUnit,*) nhep
      WRITE(m_EvtUnit,*)  3,  11, 0, 0, Q1(1), Q1(2), Q1(3), mass_e
      WRITE(m_EvtUnit,*)  3, -11, 0, 0, P1(1), P1(2), P1(3), mass_e
      WRITE(m_EvtUnit,*)  2,  23, 0, 0, -SUM(1), -SUM(2), -SUM(3),MZstar
      DO I=1,NPHOT
         WRITE(m_EvtUnit,*) 1, 22, 0 ,0, PHOT(I,1), PHOT(I,2), PHOT(I,3), 0
      ENDDO
      WRITE(m_EvtUnit,*) 1,  11, 0, 0, Q2(1), Q2(2), Q2(3), mass_e
      WRITE(m_EvtUnit,*) 1, -11, 0, 0, P2(1), P2(2), P2(3), mass_e


      
      ! nphotons+5, 2 incomming electrons 2 outgoing electrons and 1 photon/z
!1    WRITE(m_EvtUnit,*) nhep, 9999, 1., CMSENE, 1./AlphaQED, AlphaQCD 
      !WRITE(m_EvtUnit,*) 11, -1, 0, 0, 0, 0., 0., p_e, CMSENE/2., mass_e, 0 , 0
      !WRITE(m_EvtUnit,*) -11, -1, 0, 0, 0, 0., 0., -p_e, CMSENE/2., mass_e, 0 , 0 
!1    WRITE(m_EvtUnit,*)  11, -1, 0, 0, 0, 0, Q1(1), Q1(2), Q1(3), Q1(4), mass_e, 0 , 0  
!1    WRITE(m_EvtUnit,*) -11, -1, 0, 0, 0, 0, P1(1), P1(2), P1(3), P1(4), mass_e, 0 , 0 
!1    WRITE(m_EvtUnit,*)  23,  2, 1, 2, 0, 0, -SUM(1), -SUM(2), -SUM(3), CMSENE-SUM(4),MZstar, 0, 0
! what is written above is a dummy, it's there just to make LHE file, but it's not actuall simulation
!1    DO I=1,NPHOT
!1       WRITE(m_EvtUnit,*) 22, 1 , 1, 2, 0 ,0, PHOT(I,1), PHOT(I,2), PHOT(I,3), PHOT(I,4), 0., 0, 0 
!1    ENDDO
!1    WRITE(m_EvtUnit,*) 11,  1, 3, 3, 0, 0, Q2(1), Q2(2), Q2(3), Q2(4), mass_e, 0 , 0
!1    WRITE(m_EvtUnit,*) -11, 1, 3, 3, 0, 0, P2(1), P2(2), P2(3), P2(4), mass_e, 0 , 0
      
      
!1    WRITE(m_EvtUnit,'(a)')  '</event>'

      END
      
