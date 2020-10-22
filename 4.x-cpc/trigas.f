      SUBROUTINE TRIGAS0
     $   (NPHI,TH1N,TH2N,TH1W,TH2W,XWIDE,XNARR,XMIX1,XMIX2)
C     **********************************************************
C Idealized simple trigger on BARE/DRESSED final electrons.
C It is made entirely on P2, Q2 momenta -- 
C for BHLUMI and OLDBIS they are normal 'bare' electrons 
C while for LUMLOG they are by definition 'dressed'.
C This trigger is used in TH.5888 and TH.9555.
C Input:  TH1N,TH2N,TH1W,TH2W  theta limits of the Narrow/Wide cones
C         NPHI        number of phi-sectors (dummy)
C Output: XWIDE,XNARR,XMIX1,XMIX2 = 0 or V=1-S'/S, depending 
C         whether theta's fall into corresponding enery range
C         RANGE = N-N, W-W, N-W, W-N or does not.
C This looks a bit like puting right hand to left pocket but
C it is done in order to be in compatibility with more complicated
C realistic calorimetric triggers,  see TRIGAS1.
C What is important to remember is that to accept an event we shall 
C require 0 X>V_min in the routine calling TRIGAS0. 
C The condition X>V_min is in fact the logical multiplication
C (V>V_min) .AND. (theta's in the given RANGE) i.e. definition
C of our present trigger!         
C N.B. This trigger is identical to trigger in CUTCUT in
C version of LUMLOG distributed in January 91.
C     ****************************************** 
      IMPLICIT REAL*8(A-H,O-Z) 
      SAVE
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT 
      LOGICAL LFORW,LFORN,LBACW,LBACN
C wide/narrow sectors forward
      THET=ANGFI(P2(3),DSQRT(P2(1)**2+P2(2)**2))
      LFORW= THET.GT.TH1W.AND.THET.LT.TH2W
      LFORN= THET.GT.TH1N.AND.THET.LT.TH2N
C wide/narrow sectors backward
      THET=ANGFI(-Q2(3),DSQRT(Q2(1)**2+Q2(2)**2))
      LBACW= THET.GT.TH1W.AND.THET.LT.TH2W
      LBACN= THET.GT.TH1N.AND.THET.LT.TH2N
C ....
      SV1 = (P2(4)+Q2(4))**2 -(P2(3)+Q2(3))**2
     $     -(P2(2)+Q2(2))**2 -(P2(1)+Q2(1))**2 
      SV  = (P1(4)+Q1(4))**2 -(P1(3)+Q1(3))**2
     $     -(P1(2)+Q1(2))**2 -(P1(1)+Q1(1))**2  
      V = 1 -SV1/SV   
C ....
      XWIDE= 0D0
      XNARR= 0D0
      XMIX1= 0D0
      XMIX2= 0D0
      IF(LFORW.AND.LBACW) XWIDE= 1D0-V
      IF(LFORN.AND.LBACN) XNARR= 1D0-V
      IF(LFORW.AND.LBACN) XMIX1= 1D0-V
      IF(LFORN.AND.LBACW) XMIX2= 1D0-V
      END


      SUBROUTINE TRIGAS1
     $   (NPHI,TH1N,TH2N,TH1W,TH2W,XWIDE,XNARR,XMIX1,XMIX2)
C     **********************************************************
C Idealized exper. CALORIMETRIC trigger on dressed final electrons.
C Electrons and photons not distinguished!
C It is described in detail in TH.6118.
C Input:  TH1N,TH2N,TH1W,TH2W  theta limits of the Narrow/Wide cone
C         NPHI                 number of phi-sectors
C Output: XWIDE,XNARR,XMIX1,XMIX2, each variable X parametrizes
C         enery 'cought' in the corresponding N-N, W-W, N-W, W-N 
C         angular range, X>X_min is the energy cut contition
C     ****************************************** 
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      PARAMETER( PI = 3.1415926535897932D0 )
      PARAMETER( NMX1 = 20)        
      LOGICAL LANGW,LANGN,LPHI
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      DIMENSION PC(100,4)
      DIMENSION PCW1(NMX1,4),PCN1(NMX1,4)
      DIMENSION PCW2(NMX1,4),PCN2(NMX1,4)
      DATA ICONT /0/

      IF(icont.eq.0) THEN
        icont=icont+1
        WRITE(6,*) '*******************************'
        WRITE(6,*) '***     TH-2128 Trigger     ***'
        WRITE(6,*) '*******************************'
        WRITE(6,'(5A12/4F12.9,I12)')
     $  ' TH1W',' TH2W',' TH1N',' TH2N','NPHI',
     $    TH1W ,  TH2W ,  TH1N ,  TH2N , NPHI
      ENDIF

C Beam energy
      ENE = P1(4)
C Final electrons and photons not distinguished
      DO 10 K=1,4
      PC(1,K)=P2(K)
   10 PC(2,K)=Q2(K)
      DO 20 I=1,NPHOT
      DO 20 K=1,4
   20 PC(2+I,K)=PHOT(I,K)
      NP = NPHOT+2
      DO 40 I=1,NMX1
      DO 40 K=1,4
      PCW1(I,K)=0D0
      PCN1(I,K)=0D0
      PCW2(I,K)=0D0
   40 PCN2(I,K)=0D0       
C
C Collecting energies in calorimeter sectors
C
      DO 100 I=1,NP
C wide/narrow sectors forward
      THET=ANGFI(PC(I,3),DSQRT(PC(I,1)**2+PC(I,2)**2))
      PHI =ANGFI(PC(I,1),PC(I,2))
      LANGW= THET.GT.TH1W.AND.THET.LT.TH2W
      LANGN= THET.GT.TH1N.AND.THET.LT.TH2N
      DO 60 JPHI=1,NPHI
      PHI1 = (JPHI-1)*(2D0*PI/NPHI)
      PHI2 =    JPHI *(2D0*PI/NPHI)
      LPHI= PHI .GT.PHI1.AND. PHI.LT.PHI2
      IF(LANGW.AND.LPHI) THEN
        DO 50 K=1,4
   50   PCW1(JPHI,K)=PCW1(JPHI,K)+ PC(I,K)
      ENDIF
      IF(LANGN.AND.LPHI) THEN
        DO 51 K=1,4
   51   PCN1(JPHI,K)=PCN1(JPHI,K)+ PC(I,K)
      ENDIF
   60 CONTINUE
C wide/narrow sectors backward
      THET=ANGFI(-PC(I,3),DSQRT(PC(I,1)**2+PC(I,2)**2))
      PHI =ANGFI(-PC(I,1),-PC(I,2))
      LANGW= THET.GT.TH1W.AND.THET.LT.TH2W
      LANGN= THET.GT.TH1N.AND.THET.LT.TH2N
      DO 80 JPHI=1,NPHI
      PHI1 = (JPHI-1)*(2D0*PI/NPHI)
      PHI2 =    JPHI *(2D0*PI/NPHI)
      LPHI= PHI .GT.PHI1.AND. PHI.LT.PHI2
      IF(LANGW.AND.LPHI) THEN
        DO 70 K=1,4
   70   PCW2(JPHI,K)=PCW2(JPHI,K)+ PC(I,K)
      ENDIF
      IF(LANGN.AND.LPHI) THEN
        DO 71 K=1,4
   71   PCN2(JPHI,K)=PCN2(JPHI,K)+ PC(I,K)
      ENDIF  
   80 CONTINUE
  100 CONTINUE   
     
C at least one coincidences in a pair of opposite calorimetric blocks
      XWIDE= 0D0
      XNARR= 0D0
      XMIX1= 0D0
      XMIX2= 0D0

      DO 150 JPHI=1,NPHI
C minimum of energies deposed in the opposite blocks,
C we will require energies in opposite block above certain minimum
      EWW   = DMIN1(PCW1(JPHI,4),PCW2(JPHI,4))
      ENN   = DMIN1(PCN1(JPHI,4),PCN2(JPHI,4))
      EWN   = DMIN1(PCW1(JPHI,4),PCN2(JPHI,4))
      ENW   = DMIN1(PCN1(JPHI,4),PCW2(JPHI,4))
C maximum over pairs of blocks, 
C we will ask at least one pair of blocks excited above certain min.
      XWIDE = DMAX1(XWIDE,EWW/ENE)
      XNARR = DMAX1(XNARR,ENN/ENE)
      XMIX1 = DMAX1(XMIX1,EWN/ENE)
      XMIX2 = DMAX1(XMIX2,ENW/ENE)
  150 CONTINUE
      END 




************************************************************************
*                   ALEPH LCAL
*
* Narrow range approximately   3.3-6.3 deg.
* Energy default cut
*      min(Ecl_1+Ecl_2)/E_beam  > 0.44 and
*      (Ecl_1+Ecl_2)/(2*E_beam) > 0.6
************************************************************************
      SUBROUTINE TRIGAS2
     $   (TH1N,TH2N,TH1W,TH2W,NPHI,XWIDE,XNARR,XMIX1,XMIX2)
C     **********************************************************
C MODIFIED TRIGAS1 FOR ALEPH GEOMETRICAL CUTS
C Idealized exper. CALORIMETRIC trigger on dressed final electrons.
C Electrons and photons not distinguished!
C     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0 )
      LOGICAL LANGW,LANGN,LPHI
      LOGICAL LNARROW,LWIDE,KEEPIT
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      REAL *4 X4,Y4
      DIMENSION PC(100,4)
      DIMENSION PCW1(4),PCN1(4)
      DIMENSION PCW2(4),PCN2(4)
      DATA ICONT /0/

      IF(icont.eq.0) THEN
        icont=icont+1
        WRITE(6,*) '*******************************'
        WRITE(6,*) '***     ELCAL  Trigger     ***'
        WRITE(6,*) '*******************************'
        WRITE(6,'(5A12/4F12.9,I12)')
     $  ' TH1W',' TH2W',' TH1N',' TH2N','NPHI',
     $    TH1W ,  TH2W ,  TH1N ,  TH2N , NPHI
      ENDIF

C Beam energy
      ENE = P1(4)
C Final electrons and photons not distinguished
      DO 10 K=1,4
      PC(1,K)=P2(K)
   10 PC(2,K)=Q2(K)
      DO 20 I=1,NPHOT
      DO 20 K=1,4
   20 PC(2+I,K)=PHOT(I,K)
      NP = NPHOT+2
      DO 40 K=1,4
      PCW1(K)=0D0
      PCN1(K)=0D0
      PCW2(K)=0D0
   40 PCN2(K)=0D0
C
C Collecting energies in calorimeter sectors
C
      Z0 = 2800.D0
      DO 100 I=1,NP
C wide/narrow sectors forward
C Staszek's angels
      THETA=ANGFI(ABS(PC(I,3)),DSQRT(PC(I,1)**2+PC(I,2)**2))
      PHI  =ANGFI(PC(I,1),PC(I,2))

C Bolek's angels
c****      THETA = ACOS(ABS(PC(I,3))/PC(I,4))
C***      PHIA  = ATG( PC(I,1),PC(I,2)) - PI/2.
      X     = Z0*PC(I,1)/ABS(PC(I,3))
      Y     = Z0*PC(I,2)/ABS(PC(I,3))

      LNARROW = .FALSE.
      LWIDE   = .TRUE.
      IF( KEEPIT(X,Y,5) ) LNARROW = .TRUE.
      CALL LOOSEC(THETA,PHI,LWIDE)
      LANGW = .NOT. LWIDE
      LANGN = LNARROW
C**      X4 = X
C**      Y4 = Y
C**      IF(LANGN) CALL HF2(100,X4,Y4,1.)
C**      IF(LANGW) CALL HF2(200,X4,Y4,1.)
C**      TYPE *,I,PC(I,1),PC(I,2),THETA,THET,PHIA,PHI,X,Y,
C**     &       LNARROW,LWIDE,langw,LANGN

      IF(PC(I,3) .GT. 0D0) THEN

      IF(LANGW) THEN
        DO 50 K=1,4
   50   PCW1(K)=PCW1(K)+ PC(I,K)
      ENDIF
      IF(LANGN) THEN
        DO 51 K=1,4
   51   PCN1(K)=PCN1(K)+ PC(I,K)
      ENDIF

      ELSE

C wide/narrow sectors backward
      IF(LANGW) THEN
        DO 70 K=1,4
   70   PCW2(K)=PCW2(K)+ PC(I,K)
      ENDIF
      IF(LANGN) THEN
        DO 71 K=1,4
   71   PCN2(K)=PCN2(K)+ PC(I,K)
      ENDIF

      ENDIF
C**      TYPE *,PCW1(4),PCW2(4),PCN1(4),PCN2(4)
  100 CONTINUE

C at least one coincidences in a pair of opposite calorimetric blocks
      XWIDE= 0D0
      XNARR= 0D0
      XMIX1= 0D0
      XMIX2= 0D0

c[[[[[[[[[ Replaced by S.J.
c      IF(PCW1(4)/ENE .GT. 0.44D0  .AND.
c     &   PCW2(4)/ENE .GT. 0.44D0  .AND.
c     &  (PCW1(4)+PCW2(4))/(2D0*ENE) .GT. 0.6D0) XWIDE = PCW1(4)+PCW2(4)
c      IF(PCN1(4)/ENE .GT. 0.44D0  .AND.
c     &   PCN2(4)/ENE .GT. 0.44D0  .AND.
c     &  (PCN1(4)+PCN2(4))/(2D0*ENE) .GT. 0.6D0) XNARR = PCN1(4)+PCN2(4)
c      IF(PCW1(4)/ENE .GT. 0.44D0  .AND.
c     &   PCN2(4)/ENE .GT. 0.44D0  .AND.
c     &  (PCW1(4)+PCN2(4))/(2D0*ENE) .GT. 0.6D0) XMIX1 = PCW1(4)+PCN2(4)
c      IF(PCN1(4)/ENE .GT. 0.44D0  .AND.
c     &   PCW2(4)/ENE .GT. 0.44D0  .AND.
c     &  (PCN1(4)+PCW2(4))/(2D0*ENE) .GT. 0.6D0) XMIX2 = PCN1(4)+PCW2(4)
!
! Replaced by S.J.
!
      IF(PCW1(4)/ENE .GT. 0.44D0  .AND. PCW2(4)/ENE .GT. 0.44D0 
     $        )  XWIDE = (PCW1(4)+PCW2(4))/(2*Ene)
      IF(PCN1(4)/ENE .GT. 0.44D0  .AND. PCN2(4)/ENE .GT. 0.44D0  
     &        )  XNARR = (PCN1(4)+PCN2(4))/(2*Ene)
      IF(PCW1(4)/ENE .GT. 0.44D0  .AND. PCN2(4)/ENE .GT. 0.44D0  
     $        )  XMIX1 = (PCW1(4)+PCN2(4))/(2*Ene)
      IF(PCN1(4)/ENE .GT. 0.44D0  .AND. PCW2(4)/ENE .GT. 0.44D0 
     $        )  XMIX2 = (PCN1(4)+PCW2(4))/(2*Ene)
c]]]]]]]]]

      END



C***************** ADDED BY BOLEK 20 - DEC - 1990 *****************

      LOGICAL FUNCTION KEEPIT(X,Y,METH)
C----------------------------------------------------------------------
C!  -
C!
C!   Author   :- John Renner Hansen    23-JAN-1990
C!
C!   Inputs:Method number:  METHOD = (5,6)
C!          X,Y of particle
C!        -
C!
C!   Outputs: KEEPIT
C!        -
C!
C!   Libraries required:
C!
C!   Description
C!   ===========
C!
C?
      IMPLICIT REAL*8(A-H,O-Z)
      REAL *8 PADS/29.75D0/
      REAL *8 X,Y,XP,YP
      REAL *8 LXP5(9)/4.5D0,4.5D0,3.5D0,2.5D0,5*1.5D0/
      REAL *8 HXP(9)/10.5D0,9.5D0,9.5D0,9.5D0,8.5D0,7.5D0
     &              ,6.5D0,5.5D0,3.5D0/
      REAL *8 LYP(10)/0D0,3.D0,4.D0,5.D0,6.D0,7.D0
     &               ,8.D0,9.D0,10.D0,11.D0/
      REAL *8 LXP6(9)/5.5D0,5.5D0,4.5D0,2.5D0,5*1.5D0/
      INTEGER*4 METH
      LOGICAL LOG1,LOG2,LOG3,LOG4
C!======================================================================
      KEEPIT = .FALSE.
C
C     DEFINE ACTIVE REGION
      XP  =  ABS(X/PADS)
      YP  =  ABS(Y/PADS)
C
      IF(METH.EQ.5) THEN
      DO I = 1,9
        LOG1 = XP.GT.LXP5(I)
        LOG2 = XP.LT.HXP(I)
        LOG3 = YP.GE.LYP(I)
        LOG4 = YP.LE.LYP(I+1)
        IF(LOG1.AND.LOG2.AND.LOG3.AND.LOG4) THEN
          KEEPIT = .TRUE.
          GOTO 999
        ENDIF
      ENDDO
      ELSE
      DO I = 1,9
        LOG1 = XP.GT.LXP6(I)
        LOG2 = XP.LT.HXP(I)
        LOG3 = YP.GT.LYP(I)
        LOG4 = YP.LT.LYP(I+1)
        IF(LOG1.AND.LOG2.AND.LOG3.AND.LOG4) THEN
          KEEPIT = .TRUE.
          GOTO 999
        ENDIF
      ENDDO
      ENDIF
  999 RETURN
      END
c********************************************************************
      SUBROUTINE LOOSEC(THETA,PHI,REJFL)
C----------------------------------------------------------------------
C!  - Bhabha selection method 5 on loose side
C!
C!   Author   :- Peter H. Hansen       27-NOV-1990
C!
C!   Inputs: THETA,PHI of cluster centroid in global system
C!   Output: REJFL = .TRUE. if rejected
C!                 = .FALSE. if accepted
C?
C!======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DXN(4),DYN(4),DZN(4),OMXN(4),OMYN(4),OMZN(4)
      DIMENSION DIST(2)
      LOGICAL REJFL
C-----------------------------------------------------------------------
      REJFL = .TRUE.
C
C Sign z and module number
      IF(THETA.GT.1.D0) THEN
        SIG  =-1.D0
        MODU = 2
      ELSE
        SIG  = 1.D0
        MODU = 4
      ENDIF
      IF(COS(PHI).LT.0.D0) MODU=MODU-1
C
C X and Y in global system
C
      X = SIG*280.D0*COS(PHI)*TAN(THETA)
      Y = SIG*280.D0*SIN(PHI)*TAN(THETA)
C
C Local system
      ZL = SIG*17.5D0
      XLOC = X - DXN(MODU) - OMZN(MODU)*Y + OMYN(MODU)*ZL
      YLOC = Y - DYN(MODU) - OMXN(MODU)*ZL + OMZN(MODU)*X
C
C Fold into first quadrant
      XLOC = ABS(XLOC)
      YLOC = ABS(YLOC)
C
C Find distance to edges (as in LCLUTW, but database hardwired)
      DIST(1) = XLOC - 1.9D0
      DIST(2) = 24.5D0
      IF(YLOC.GT.8.4D0+1.2D0) THEN
        DIS = XLOC-(17.5D0-YLOC)/0.75D0
        IF(DIS.LT.DIST(1)) DIST(1)=DIS
        DIST(2) = YLOC-(17.5D0-XLOC*0.75D0)
      ELSE
        DIST(1) = XLOC-11.9D0
      ENDIF
C

C Now make the cut on the inner boundary
C and on the outer boundary
      IF(THETA.GT.1.D0) THETA = 3.14159D0-THETA
      IF(DIST(1).GT.1.D0  .AND.
     &   DIST(2).GT.0.75D0  .AND.
     &   THETA.LT.0.125D0)  REJFL = .FALSE.
C**      TYPE *,DIST(1),DIST(2),THETA,REJFL
C
  999 CONTINUE
      END



************************************************************************
*                   DELPHI VSAT
************************************************************************
*Resent-Date:  Sat, 26 Nov 94 22:39:21 WET
*Resent-From: JADACH@crnvma.cern.ch
*Resent-To: Stanislaw Jadach <JADACH@hephp01.phys.utk.edu>
*Return-Path: <read@VUOEPA.UIO.NO>
*Date: Wed, 23 Nov 1994 13:39:43 +0200
*From: read@VUOEPA.UIO.NO (Alex Read, U. Oslo, (47)22.85.5062)
*To: jadach@crnvma.cern.ch
*Subject: Trigger Code DELPHI
*X-Vms-To: smtp%'jadach@cernvm.cern.ch'
*---------------------------Original message----------------------------
*Staszek,
*I have test the new 'stand-alone' version of the DELPHI SAT
*acceptance versus my standard software and find identical
*results for 100k events. This should be a good guarentee
*that the code I send you is OK. The two routines follow.
*DSATFILL is called by LOGICAL FUNCTION DSATACC(XMIN).
*XMIN is the returned value of Emin/Ebeam (REAL*4).
*
*To test the DELPHI SAT acceptance on an event in /MOMSET/:
*      .
*      .
*      REAL*4 xmin
*      LOGICAL  dsatacc,ltemp
*      EXTERNAL dsatacc
*      .
*      .
*C LTEMP is .TRUE. for events inside the DELPHI acceptance.
*C xmin is returned to you.
*      ltemp = dsatacc(xmin)
*      .
*      .
*
*Let me know if I need to do anything more or provide
*some information you need.
*
*Cheers, Alex
************************************************************************


      SUBROUTINE dsatfil
*     ******************
*  Fill something which resembles the PAW common block for
*  the standard NTUPLES. We have to play some tricks to
*  allow for muti-photon final states (BABAMC has <=1 photon).
*
*.--------------------------------------------------------
*   MOMSET as defined in BHL2CPC.FOR
*
      INTEGER nphot,nphmax
      PARAMETER (nphmax=100)
      DOUBLE PRECISION p1,p2,q1,q2,phot
      COMMON /momset/ p1(4),q1(4),p2(4),q2(4),phot(nphmax,4),nphot
*.--------------------------------------------------------
*.--------------------------------------------------------
*   Similar to common block used in PAW analysis of DELPHI
*   SAT NTUPLES. The photon arrays are the same
*   size as for /MOMSET/ in contrast to the limit of the
*   2 hardest photons written to the NTUPLE.
      REAL
     +FILL    ,NRUN    ,FILE    ,IEVT    ,DATE    ,TIME    ,
     +EBEA    ,BBIT    ,SIAR    ,B1BT    ,BABA    ,NCL1    ,
     +NCL2    ,  E1    ,  R1    ,Phi1    ,Siz1    ,E1R1    ,
     +  E2    ,  R2    ,Phi2    ,Siz2    ,E2R1    ,NLUN    ,
     +RLU1    ,TLU1    ,PLU1    ,ELU1    ,RLU2    ,TLU2    ,
     +PLU2    ,ELU2    ,RLU3    ,TLU3    ,PLU3    ,ELU3
      COMMON/PAWIDN/IDNEVT,VIDN1,VIDN2,VIDN3,VIDN(10),
     +FILL    ,NRUN    ,FILE    ,IEVT    ,DATE    ,TIME    ,
     +EBEA    ,BBIT    ,SIAR    ,B1BT    ,BABA    ,NCL1    ,
     +NCL2    ,  E1    ,  R1    ,Phi1    ,Siz1    ,E1R1    ,
     +  E2    ,  R2    ,Phi2    ,Siz2    ,E2R1    ,NLUN    ,
     +RLU1    ,TLU1    ,PLU1    ,ELU1    ,RLU2    ,TLU2    ,
     +PLU2    ,ELU2    ,
     +RLU3(nphmax) ,TLU3(nphmax) ,PLU3(nphmax) ,ELU3(nphmax)
*.--------------------------------------------------------
*
      DOUBLE PRECISION wtmodc
      PARAMETER ( PI = 3.141592654 )
      PARAMETER ( RADIAN = PI/180. )
      REAL       zcal
      PARAMETER (zcal = 231.8)
*
*......Initialisation
*
        ebea = p1(4)
        elu1 = 0.
        rlu1 = 0.
        plu1 = 0.
        tlu1 = 0.
        elu2 = 0.
        rlu2 = 0.
        plu2 = 0.
        tlu2 = 0.
        DO i=1,nphot
          elu3(i) = 0.
          rlu3(i) = 0.
          plu3(i) = 0.
          tlu3(i) = 0.
        ENDDO
*
*......Electron (hits side A/1 instead of C as in reality)
*
      PX     = Q2(1)
      PY     = Q2(2)
      PZ     = Q2(3)
      ELU1   = Q2(4)
      TANT1  = SQRT(PX*PX+PY*PY)/PZ
      TLU1   = ATAN(TANT1)/RADIAN
      IF ( TLU1.LT.0 ) TLU1 = TLU1 + 180.
      RLU1   = zcal*ABS(TANT1)
      PLU1   = PROXIM( ATAN2(PY,PX), PI )/RADIAN
*
*......Positron (hits side C/2 instead of A as in reality)
*
      PX     = P2(1)
      PY     = P2(2)
      PZ     = P2(3)
      ELU2   = P2(4)
      TANT2  = SQRT(PX*PX+PY*PY)/PZ
      TLU2   = ATAN(TANT2)/RADIAN
      IF ( TLU2.LT.0 ) TLU2 = TLU2 + 180.
      RLU2   = zcal*ABS(TANT2)
      PLU2   = PROXIM( ATAN2(PY,PX), PI )/RADIAN
*
*......Photons
*
      DO 100 i=1,nphot
        NLUN = NLUN + 1
        PX     = PHOT(I,1)
        PY     = PHOT(I,2)
        PZ     = PHOT(I,3)
        ELU3(i)= PHOT(I,4)
        TANT3  = SQRT(PX*PX+PY*PY)/PZ
        TLU3(i)= ATAN(TANT3)/RADIAN
        IF ( TLU3(i).LT.0 ) TLU3(i) = TLU3(i) + 180.
        RLU3(i)= zcal*ABS(TANT3)
        PLU3(i)= PROXIM( ATAN2(PY,PX), PI )/RADIAN
  100 CONTINUE
*
 999  RETURN
      END

      LOGICAL FUNCTION DSATACC(XMIN)
*     ******************************
*--   Emulation of MC using 4-vectors
*--   New cuts                          26-SEP-1990   Mogens Dam
*--   Update to give better description 05-FEB-1991   Mogens Dam
*--   This is version to treat multi-photon final states (BHLUMI)
*     24.09.91, Alex Read
*
* In each arm of the calorimeter the primary cluster must be located.
* Photons inside the acceptance are either added to the electron (if
* it is inside the acceptance itself) or to a secondary cluster which
* will be composed of only photons. This is an approximation since it
* is possible that there are two photons in addition to the electron
* and they may not be close to each other.
*
* The cuts are emulated once the primary cluster in each arm has been
* chosen.
*
*.--------------------------------------------------------
*   MOMSET as defined in BHL2CPC.FOR
*
      INTEGER nphot,nphmax
      PARAMETER (nphmax=100)
      DOUBLE PRECISION p1,p2,q1,q2,phot
      COMMON /momset/ p1(4),q1(4),p2(4),q2(4),phot(nphmax,4),nphot
*.--------------------------------------------------------
*.--------------------------------------------------------
*   Similar to common block used in PAW analysis of DELPHI
*   SAT NTUPLES. The photon arrays are the same
*   size as for /MOMSET/ in contrast to the limit of the
*   2 hardest photons written to the NTUPLE.
      REAL
     +FILL    ,NRUN    ,FILE    ,IEVT    ,DATE    ,TIME    ,
     +EBEA    ,BBIT    ,SIAR    ,B1BT    ,BABA    ,NCL1    ,
     +NCL2    ,  E1    ,  R1    ,Phi1    ,Siz1    ,E1R1    ,
     +  E2    ,  R2    ,Phi2    ,Siz2    ,E2R1    ,NLUN    ,
     +RLU1    ,TLU1    ,PLU1    ,ELU1    ,RLU2    ,TLU2    ,
     +PLU2    ,ELU2    ,RLU3    ,TLU3    ,PLU3    ,ELU3
      COMMON/PAWIDN/IDNEVT,VIDN1,VIDN2,VIDN3,VIDN(10),
     +FILL    ,NRUN    ,FILE    ,IEVT    ,DATE    ,TIME    ,
     +EBEA    ,BBIT    ,SIAR    ,B1BT    ,BABA    ,NCL1    ,
     +NCL2    ,  E1    ,  R1    ,Phi1    ,Siz1    ,E1R1    ,
     +  E2    ,  R2    ,Phi2    ,Siz2    ,E2R1    ,NLUN    ,
     +RLU1    ,TLU1    ,PLU1    ,ELU1    ,RLU2    ,TLU2    ,
     +PLU2    ,ELU2    ,
     +RLU3(nphmax) ,TLU3(nphmax) ,PLU3(nphmax) ,ELU3(nphmax)
*.--------------------------------------------------------
      INTEGER npar(100)
      DOUBLE PRECISION xpar(100)
*
      REAL xmin
      REAL xcal1(2),ycal1(2),ecal1(2),xcal2(2),ycal2(2),ecal2(2)
      REAL DISTCT,pf1(3),pf2(3)
      REAL xdummy
      DOUBLE PRECISION wtri,wtmodc,wtsel
      DATA DISTCT / 10.0 /
*
      PARAMETER ( PI = 3.141592654 )
      PARAMETER ( RADIAN = PI/180. )
*
      DSATACC = .FALSE.
*
*......Fill the common block used in the DELPHI PAW-based
*      analysis.
*
      CALL dsatfil
*
*=========================================================
*--  Masked calorimeter; electron
*=========================================================
*
*......Initialize sums for photon cluster
*
      xcal2(2) = 0.
      ycal2(2) = 0.
      ecal2(2) = 0.
*
*......Initialize sums for electron cluster. The curvature due to the
*      magnetic field is corrected for and the acceptance is applied.
*
      xcal2(1) = 0.
      ycal2(1) = 0.
      ecal2(1) = 0.
      PHIR2    = PLU2-0.55*45.55/MAX(ELU2,0.001)
      PHIR2R   = PHIR2*RADIAN
      IF ( (PHIR2.GT. 75.7.AND.PHIR2.LT.104.3) .OR.
     &     (PHIR2.GT.255.7.AND.PHIR2.LT.284.3) .OR.
     &     (RLU2.LT.13.0 .OR. RLU2.GT.34.0)         ) THEN
        ECAL2(1) = 0.
        X2 = 0.
        Y2 = 0.
      ELSE
        xcal2(1) = RLU2*COS(PHIR2R)*elu2
        ycal2(1) = RLU2*SIN(PHIR2R)*elu2
        ecal2(1) = elu2
        X2 = RLU2*COS(PHIR2R)
        Y2 = RLU2*SIN(PHIR2R)
      ENDIF
*
*......Now loop over the photons. Check that they are in the
*      same calorimeter arm, and if they are inside the calorimeter
*      acceptance add them to the electron (if it is also inside the
*      acceptance) or make a second 'neutral only' cluster.
*
      DO 100 i=1,nphot
        PHIR3 = PLU3(i)
        PHIR3R = PHIR3*RADIAN
        X3     = RLU3(i)*COS(PHIR3R)
        Y3     = RLU3(i)*SIN(PHIR3R)
        IF(TLU3(i).LT.90.) THEN
          IF ( .NOT.( (PHIR3.GT. 75.7.AND.PHIR3.LT.104.3) .OR.
     &         (PHIR3.GT.255.7.AND.PHIR3.LT.284.3) .OR.
     &         (RLU3(i).LT.13.0 .OR. RLU3(i).GT.34.0))     ) THEN
            IF(ECAL2(1).GT.0.0) THEN
              DIST = SQRT((X2-X3)**2+(Y2-Y3)**2)
              IF( DIST.LE.DISTCT ) THEN
                ECAL2(1) = ECAL2(1) + ELU3(i)
                XCAL2(1) = XCAL2(1) + X3*ELU3(i)
                YCAL2(1) = YCAL2(1) + Y3*ELU3(i)
              ELSE
                ECAL2(2) = ECAL2(2) + ELU3(i)
                XCAL2(2) = XCAL2(2) + X3*ELU3(i)
                YCAL2(2) = YCAL2(2) + Y3*ELU3(i)
              ENDIF
            ELSE
              ECAL2(2) = ECAL2(2) + ELU3(i)
              XCAL2(2) = XCAL2(2) + X3*ELU3(i)
              YCAL2(2) = YCAL2(2) + Y3*ELU3(i)
            ENDIF
          ENDIF
        ENDIF
  100 CONTINUE
*
*......Choose the highest energy cluster and give energy, phi, radius
*      (ECAL2(1), PCAL2, RCAL2) for this 'primary cluster' to be used
*      in the cuts later.
*
      If(ecal2(1).GT.ecal2(2)) THEN
        xcal2(1) = xcal2(1)/ecal2(1)
        ycal2(1) = ycal2(1)/ecal2(1)
      ELSE IF(ecal2(2).GT.0.0) THEN
        xcal2(1) = xcal2(2)/ecal2(2)
        ycal2(1) = ycal2(2)/ecal2(2)
        ecal2(1) = ecal2(2)
      ELSE
        xcal2(1) = 0.
        ycal2(1) = 0.
      ENDIF
      IF(ecal2(1).GT.0.0) THEN
        RCAL2 = SQRT(XCAL2(1)*XCAL2(1)+YCAL2(1)*YCAL2(1))
        PCAL2 = ATAN2(YCAL2(1),XCAL2(1))/RADIAN
        IF ( PCAL2.LT.0 )  PCAL2=PCAL2+360.
      ELSE
        rcal2 = 0.
        pcal2 = 0.
      ENDIF
*
*=========================================================
*--  Now repeat this for non-masked calorimeter; positron
*=========================================================
*
*......Initialize sums for photon cluster
*
      xcal1(2) = 0.
      ycal1(2) = 0.
      ecal1(2) = 0.
*
*......Initialize sums for electron cluster. The curvature due to the
*      magnetic field is corrected for and the acceptance is applied.
*
      xcal1(1) = 0.
      ycal1(1) = 0.
      ecal1(1) = 0.
      PHIR1  = PLU1+0.55*45.55/MAX(0.001,ELU1)
      PHIR1R = PHIR1*RADIAN
      IF ( RLU1.LT.10.0 .OR. RLU1.GT.34.0 ) THEN
        ECAL1(1) = 0.
        X1 = 0.
        Y1 = 0.
      ELSE
        xcal1(1) = RLU1*COS(PHIR1R)*elu1
        ycal1(1) = RLU1*SIN(PHIR1R)*elu1
        ecal1(1) = elu1
        X1 = RLU1*COS(PHIR1R)
        Y1 = RLU1*SIN(PHIR1R)
      ENDIF
*
*......Now loop over the photons. Check that they are in the
*      same calorimeter arm, and if they are inside the calorimeter
*      acceptance add them to the electron (if it is also inside the
*      acceptance) or make a second 'neutral only' cluster.
*
      DO 200 i=1,nphot
        PHIR3 = PLU3(i)
        PHIR3R = PHIR3*RADIAN
        X3     = RLU3(i)*COS(PHIR3R)
        Y3     = RLU3(i)*SIN(PHIR3R)
        IF(TLU3(i).GT.90.) THEN
          IF( RLU3(i).GE.10.0 .AND. RLU3(i).LE.34.0) THEN
            IF(ECAL1(1).GT.0.0) THEN
              DIST = SQRT((X1-X3)**2+(Y1-Y3)**2)
              IF( DIST.LE.DISTCT ) THEN
                ECAL1(1) = ECAL1(1) + ELU3(i)
                XCAL1(1) = XCAL1(1) + X3*ELU3(i)
                YCAL1(1) = YCAL1(1) + Y3*ELU3(i)
              ELSE
                ECAL1(2) = ECAL1(2) + ELU3(i)
                XCAL1(2) = XCAL1(2) + X3*ELU3(i)
                YCAL1(2) = YCAL1(2) + Y3*ELU3(i)
              ENDIF
            ELSE
              ECAL1(2) = ECAL1(2) + ELU3(i)
              XCAL1(2) = XCAL1(2) + X3*ELU3(i)
              YCAL1(2) = YCAL1(2) + Y3*ELU3(i)
            ENDIF
          ENDIF
        ENDIF
  200 CONTINUE
*
*......Choose the highest energy cluster and give energy, phi, radius
*      for this 'primary cluster' to be used in the cuts later.
*
      If(ecal1(1).GT.ecal1(2)) THEN
        xcal1(1) = xcal1(1)/ecal1(1)
        ycal1(1) = ycal1(1)/ecal1(1)
      ELSE IF(ecal1(2).GT.0.0) THEN
        xcal1(1) = xcal1(2)/ecal1(2)
        ycal1(1) = ycal1(2)/ecal1(2)
        ecal1(1) = ecal1(2)
      ENDIF
      IF(ecal1(1).GT.0.0) THEN
        RCAL1 = SQRT(XCAL1(1)*XCAL1(1)+YCAL1(1)*YCAL1(1))
        PCAL1 = ATAN2(YCAL1(1),XCAL1(1))/RADIAN
        IF ( PCAL1.LT.0 )  PCAL1=PCAL1+360.
      ELSE
        rcal1 = 0.
        pcal1 = 0.
      ENDIF
*
*......Make the final acceptance cuts on the primary 'clusters'
*
      ACOP = PCAL1-PCAL2-SIGN(180.,PCAL1-PCAL2)
c[[[[[[[ S. Jadach: remove energy cut from dsatacc
c      IF (  MIN(ECAL1(1),ECAL2(1)).GT.0.65*EBEA     .AND.
c     &      ABS(ACOP).LT.20.                        .AND.
c     &      RCAL2.LT.29.7                           .AND.
c     &      RCAL1.GT.12.25 .AND. RCAL1.LT.33.0       ) THEN
c        dsatacc = .TRUE.
c      ENDIF
      dsatacc =
     &      ABS(ACOP).LT.20.                        .AND.
     &      RCAL2.LT.29.7                           .AND.
     &      RCAL1.GT.12.25 .AND. RCAL1.LT.33.0
c]]]]]]]
*
*......Return Emin/Ebeam
*
      xmin = MIN(ECAL1(1),ECAL2(1))/EBEA
*
  999 END


      FUNCTION PROXIM (ANGL,STAND)
C     ****************************
C CERN PROGLIB# B102    PROXIM          .VERSION KERNFOR  4.21  890323
C ORIG. 15/03/68 JZ, re-write 8/03/89 K.S.Koelbig
C
      PARAMETER (PI2 = 6.28318 53071 79586D0, RPI2 =1/PI2)

      PROXIM = ANGL + PI2*ANINT(RPI2*(STAND-ANGL))
      RETURN
      END
