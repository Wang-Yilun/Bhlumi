      SUBROUTINE BHLUMI(MODE,XPAR,NPAR)        
*     *********************************
!
!
!   BBBBBBB    BBB   BBB  BBB      BBB  BBB  BBB     BBB   BBB
!   BBB  BBB   BBB   BBB  BBB      BBB  BBB  BBBB   BBBB   BBB
!   BBB  BBB   BBB   BBB  BBB      BBB  BBB  BBBBB BBBBB   BBB
!   BBBBBB     BBBBBBBBB  BBB      BBB  BBB  BBB BBB BBB   BBB
!   BBBBBBBBB  BBBBBBBBB  BBB      BBB  BBB  BBB  B  BBB   BBB
!   BBB  BBBB  BBB   BBB  BBB  BB  BBB  BBB  BBB     BBB   BBB
!   BBBBBBBBB  BBB   BBB  BBBBBBB  BBB  BBB  BBB     BBB   BBB
!   BBBBBBBB   BBB   BBB  BBBBBBB   BBBBBB   BBB     BBB   BBB
!
!
!======================================================================
!======================================================================
!======================================================================
!===============             B H L U M I            ===================
!======================================================================
!======================================================================
!=============== MONTE CARLO FOR SMALL ANGLE BHABHA ===================
!===============            VERSION 4.04            ===================
!======================================================================
!===============    July 1995 - Sept. 1996          ===================
!======================================================================
!======================================================================
!=======================     AUTHORS      =============================
!== S. Jadach, W. Placzek, E. Richter-Was, B.F.L. Ward  and Z. Was   ==
!======================================================================
!======================================================================
!           
! The description of the usage of the program can be found
! in refs. [1] and [3].
!
! The program contains three subgenerators:
!   (1)  BHLUM4 
!    Multiphoton generator with Yennie-Frautschi-Suura SECOND order 
!    exponentiation based on refs. [1-7]
!   (2)  LUMLOG
!    Leading-Logarithmic generator with collinear emission of photons,
!    QED matrix element up to third order described in refs. [1-2],[3],[8]
!   (3)  OLDBIS
!    Modified version of the first order clasical generator OLDBAB,
!    see refs. [9-10]
!
! [1] S. Jadach, W. Placzek, E. Richter-Was, B.F.L. Ward  and Z. Was,
!     CERN preprint TH-96-158, June. 1996, submitted to Comp. Phys. Commun.
! [2] S. Jadach and B.F.L. Ward,
!     CERN preprint TH-96-156, June. 1996,  Phys. Lett. in print.
! [3] S. Jadach, E. Richter-Was, B.F.L. Ward  and Z. Was,
!     Comp. Phys. Commun. 70 (1992) 305 (TH-6230, sept. 1991).
! [4] S. Jadach, E. Richter-Was, B.F.L. Ward  and Z. Was,
!     Phys.Lett. B268 (1991), 253 (TH-6118, June 1991). 
! [5] S. Jadach and B.F.L. Ward,
!     Phys. Rev. D40 (1989) 3582.
! [6] S. Jadach, E. Richter-Was, B.F.L. Ward  and Z. Was,
!     TH-95-38; Phys. Lett., B353 (1995) 362.
! [7] S. Jadach, W. Placzek and B.F.L. Ward
!     TH-95-74; Phys. Lett. B353 (1995) 349.
! [8] S. Jadach, E. Richter-Was, B.F.L. Ward  and Z. Was
!     Phys. Lett. B260 (1991) 438, TH-5995.
! [9] S. Jadach, E. Richter-Was, B.F.L. Ward  and Z. Was
!     Phys. Lett. B253 (1991) 469, TH-5888.
! [10] F. Berends and R. Kleiss Nucl. Phys. B228 (1983) 537.  
!
! Postscript files for the above and other papers 
! available from URL http://hpjmiady.ifj.edu.pl/programs/programs.html
!
!                  IMPORTANT NOTE
!                  --------------
! The user is kindly requested to cite at least refs. [1-7]
! and any other ones from the above list if his study depends strongly 
! on the particular subgenerator.
!  
!----------------------------------------------------------------------
!                        INPUT and OUTPUT
!----------------------------------------------------------------------
! All input and output goes through parameters in 
!                 CALL BHLUMI(MODE,XPAR,NPAR)
! and through /MOMSET/ and /WGTALL/ common blocks.
! In the following we shall  briefly indicate the meaning of the
! above parameters/variables.
!
! IF( MODE =-1 ) THEN
! ===================
! Initialisation is performed, all input parameters are transfered
! through XPAR and NPAR.
! The meaning of the entries in XPAR and NPAR depends on the type of 
! the subgenerator: 
!          see tables in the source of BHLUM4.f
!          ------------------------------------ 
! and also in OLDBIS, LUMLOG and the tables in the Long-write-up. 
! In the following table we indicate very briefly
! parameters which have the same meaning for all three subgenerators
!----------------------------------------------------------------------
!  Entry    Variable   Meaning
!----------------------------------------------------------------------
!  NPAR( 1)  KeyOpt =1000*KeyGen +100*KeyREM +10*KeyWGT +KeyRND 
!                   General option switch 
!            KeyGen =1,2,3 for OLDBIS,LUMLOG,BHLUM4 sub-generators
!            KeyRem =0,1 removal/no-removal switch, see BHLUM4
!            KeyRND =1,2 type of random number generator RANMAR,RANECU
!            KeyWGT =0,1,2 for constant/variable weight WTM,
!                    see more details in tables in BHLUM4 and LUMLOG
!  NPAR( 2)  KeyRAD is QED option switch defining WTMOD weight
!                   see tables in BHLUM4, LUMLOG and OLDBIS
!  XPAR( 1)  CMSENE Total center mass energy [GeV]
!  XPAR( 2)         see tables in BHLUM4, LUMLOG and OLDBIS 
!  XPAR( 3)         see tables in BHLUM4, LUMLOG and OLDBIS 
!  XPAR( 4)         see tables in BHLUM4, LUMLOG and OLDBIS 
!  XPAR( 5)         see tables in BHLUM4, LUMLOG and OLDBIS 
!  XPAR( 6)         see tables in BHLUM4, LUMLOG and OLDBIS 
!----------------------------------------------------------------------
!
! ELSE IF( MODE = 0 ) THEN
! ========================
! Generation of the single Monte Carlo event
! The four momenta of the final state electron positron and photon
! are encoded in 
!      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
! where P1 and Q1 are four-momenta of positron and electron beams, 
! P2 and Q2 are four-momenta of outgoing positron and electron [GeV],  
! PHOT(100,4) contains list of photon four-momenta
! and NPHOT is the number of real photons in PHOT.  
! For weighted events it may be profitable to use 'paralel weights' 
! from
!      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)   
! where WTMOD is the principal model weight and another useful weights
! representing some interesting version of the QED matrix element
! can be constructed as WT= WTCRU1*WTCRU2*WTSET(J). 
! Which J is alowed and what version of the matrix element it 
! represents depends on the type of subgenerator (BHLUM4,OLDBIS,LUMLOG)
! and may be found in corresponding Tables of this source 
! code  and in Long-write-up.
! 
! ELSE IF( MODE = 1 ) THEN
! ========================
! The total cross section corresponding to generated series of event,
! i.e. resulting from MC integrartion is calculated and stored in XPAR
! and NPAR. 
! In the table below we describe their most essential entries.
! For describtion of the other entries see tables in the source code
! of the subgenerators BHLUM4, LUMLOG and OLDBIS and in Long-write-up.
!----------------------------------------------------------------------
!  Entry    Variable   Meaning
!----------------------------------------------------------------------
!  NPAR(10)  NEVGEN  Number of generated MC events
!  NPAR(20)  NEVGEN  Number of generated MC events
!  XPAR(10)    XSEC  Total x-section [nb]
!  XPAR(11)   RXSEC  The relative (statistical) error of XSEC
!  XPAR(20)          Crude total MC x-section [nb] which is necessary
!                    for rescaling histograms in run 
!                    with weighted events.
!  XPAR(21)          =0, error of XPAR(20) is zero
! ----------------------------------------------------------------------
! For constant weight option KEYWGT=0 (convevience in rescaling histos)
! we put XPAR(20,21)=XPAR(10,11).
! For MODE=1 program is typicaly called upon many times in the process of 
! rescaling histograms and therefore no output is printed.
! 
! ELSE IF( MODE = 2 ) THEN
! ========================   
! Only in this MODE=2 in addition to filling XPAR and NPAR 
! (as for MODE=1)
! the values of various x-sections are printed on the standard 
! output file.
!       
! ENDIF
! ====
!     ****************************************
      IMPLICIT REAL*8(A-H,O-Z)   
      CHARACTER*80    BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G 
      PARAMETER(
     $BXOPE =  '(//1X,15(5H=====)    )',
     $BXTXT =  '(1X,1H=,                  A48,25X,    1H=)',
     $BXL1I =  '(1X,1H=,I17,                 16X, A20,A12,A7, 1X,1H=)',
     $BXL1F =  '(1X,1H=,F17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2F =  '(1X,1H=,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXL1G =  '(1X,1H=,G17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2G =  '(1X,1H=,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXCLO =  '(1X,15(5H=====)/   )'    )   
      DIMENSION  XPAR(*),NPAR(*)
      COMMON / INOUT  / NINP,NOUT    
      COMMON / BHPAR3 / KEYRAD,KEYOPT
      SAVE   / INOUT  /, / BHPAR3 /
      SAVE   NEVG, KEYGEN
   
      IF(MODE.EQ.-1) THEN        
*     ===================        
      NINP= 15   
      NOUT= 16   
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '   '
      WRITE(NOUT,BXTXT) 'BBB   B  B  B    B   B  B    B  B'
      WRITE(NOUT,BXTXT) 'B  B  B  B  B    B   B  BB  BB  B'
      WRITE(NOUT,BXTXT) 'BBB   BBBB  B    B   B  B BB B  B'
      WRITE(NOUT,BXTXT) 'B  B  B  B  B    B   B  B    B  B'
      WRITE(NOUT,BXTXT) 'BBB   B  B  BBBB  BBB   B    B  B'
      WRITE(NOUT,BXTXT) '   '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '*   BHLUMI version 4.04         *'
      WRITE(NOUT,BXTXT) '*   June           1991  (2.01) *'
      WRITE(NOUT,BXTXT) '*   Sept           1992  (2.02) *'
      WRITE(NOUT,BXTXT) '*   January        1995  (4.00) *'
      WRITE(NOUT,BXTXT) '*   Febuary        1995  (4.01) *'
      WRITE(NOUT,BXTXT) '*   May            1995  (4.02) *'
      WRITE(NOUT,BXTXT) '*   July           1995  (4.02a)*'
      WRITE(NOUT,BXTXT) '*   June           1996  (4.03) *'
      WRITE(NOUT,BXTXT) '*   September      1996  (4.04) *'
      WRITE(NOUT,BXTXT) '*         AUTHORS               *'
      WRITE(NOUT,BXTXT) '* S. Jadach,      W. Placzek,   *'
      WRITE(NOUT,BXTXT) '* E. Richter-Was, B.F.L. Ward,  *'
      WRITE(NOUT,BXTXT) '*        and Z. Was             *'
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXCLO)  
!
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) 'This program is based on papers '
      WRITE(NOUT,BXTXT) '--------------------------------'
      WRITE(NOUT,BXTXT) 'Phys. Lett. B353 (1995) 362     '
      WRITE(NOUT,BXTXT) 'Phys. Lett. B353 (1995) 349     '
      WRITE(NOUT,BXTXT) 'Comp. Phys. Comm. 70 (1992) 305 '
      WRITE(NOUT,BXTXT) 'Phys. Lett. B268 (1991) 253     '
      WRITE(NOUT,BXTXT) 'Phys. Rev.  D40  (1989) 3582.   '
      WRITE(NOUT,BXTXT) 'Phys. Lett. B260 (1991) 438,    '
      WRITE(NOUT,BXTXT) 'Phys. Lett. B253 (1991) 469,    '
      WRITE(NOUT,BXTXT) 'Nucl. Phys. B228 (1983) 537.    '
      WRITE(NOUT,BXCLO)  
!
      NEVG=0 
      KEYOPT=NPAR(1)
      KEYGEN = MOD(KEYOPT,10000)/1000
      IF(KEYGEN.EQ.0) CALL BHLUM4(  -1,XPAR,NPAR)   
      IF(KEYGEN.EQ.1) CALL OLDBIS(  -1,XPAR,NPAR)   
      IF(KEYGEN.EQ.2) CALL LUMLOG(  -1,XPAR,NPAR)   
      IF(KEYGEN.EQ.3) CALL BHLUM4(  -1,XPAR,NPAR)   
      IF(KEYGEN.EQ.4) CALL BHLOG4(  -1,XPAR,NPAR)   
      IF(KEYGEN.EQ.5) CALL BHLOG5(  -1,XPAR,NPAR)   

      ELSEIF(MODE.EQ.0) THEN
*     ======================
      NEVG=NEVG+1  
      IF(KEYGEN.EQ.0) CALL DUMGEN
      IF(KEYGEN.EQ.1) CALL OLDBIS(   0,XPAR,NPAR)   
      IF(KEYGEN.EQ.2) CALL LUMLOG(   0,XPAR,NPAR)   
      IF(KEYGEN.EQ.3) CALL BHLUM4(   0,XPAR,NPAR)   
      IF(KEYGEN.EQ.4) CALL BHLOG4(   0,XPAR,NPAR)   
      IF(KEYGEN.EQ.5) CALL BHLOG5(   0,XPAR,NPAR)   
! clean final state common blocks if necessary (safety reason)
      CALL BHCLEN

      ELSE 
!     ==== 
      IF(KEYGEN.EQ.0) CALL DUMGEN
      IF(KEYGEN.EQ.1) CALL OLDBIS(MODE,XPAR,NPAR)   
      IF(KEYGEN.EQ.2) CALL LUMLOG(MODE,XPAR,NPAR)   
      IF(KEYGEN.EQ.3) CALL BHLUM4(MODE,XPAR,NPAR)   
      IF(KEYGEN.EQ.4) CALL BHLOG4(MODE,XPAR,NPAR)   
      IF(KEYGEN.EQ.5) CALL BHLOG5(MODE,XPAR,NPAR)   
      ENDIF 
!     =====
      END

      SUBROUTINE BHCLEN
!     *****************
! This routine prevents user from using zero weight events 
! and parellel weights when they should not be used!
!     ************************
      IMPLICIT REAL*8(A-H,O-Z)   
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300) 
      SAVE   / MOMSET /, / WGTALL /

! Parallel weights should not be used for constant weight events.
      IF(WTMOD.EQ.1D0) THEN
       DO I=1,300
         WTSET(I)=0D0
       ENDDO
! Clean final state momenta for events outside phase space 
      ELSEIF(WTCRU1*WTCRU2 .EQ.0D0 )  THEN
       DO K=1,4
         P2(K)=0D0
         Q2(K)=0D0    
       ENDDO
       NPHOT=0
       DO J=1,100
         DO K=1,4
           PHOT(J,K)=0D0
         ENDDO
       ENDDO
      ENDIF
      END

      SUBROUTINE DUMGEN
!     *****************
! This routine is dummy generator for keygen=0.
! It may be usefull if generation is for some tests done directly
! in booking routine but input is the same as in BHLUM4
      END
