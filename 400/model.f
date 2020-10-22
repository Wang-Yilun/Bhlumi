! ==================================================================
! ===================== Matrix Element =============================
! ==================================================================
      SUBROUTINE MODEL(MODE,WTM)
! ***********************************************
! Interface to various models for matrix element
! ***********************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      COMMON / BHPAR3 / KEYRAD,KEYOPT
      COMMON / INOUT  / NINP,NOUT
      SAVE   / INOUT  /, / BHPAR3 /

      KeyMod = MOD(KeyRad,100)/10
      IF (KeyMod .EQ. 1) THEN
! Version as in CPC 2.01 (1991,92)
! It is kept for backward compatibility
        CALL MODEL1(MODE)
        WTM = WTSET(2)
      ELSEIF(KeyMod .EQ. 2) THEN
! New version (1993, 1994, 1995)
        CALL MODEL2(MODE)
        WTM = WTSET(1)
      ELSE
        WRITE(NOUT,*)   ' +++++ MODEL: wrong keyMod=',KeyMod
        STOP
      ENDIF
      END


      SUBROUTINE model2(mode)
!     ***********************
!=================================================================!
!                                                                 !
!    New version  of O(alf2) matrix element started May 92        !
!    Now, (Jan. 94) with two examples of O(alf2) matrix element   !
!    Version Jan 95 features also matrix elements without expon.  !
!                                                                 !
!=================================================================!
!     ************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / TRANSR / TRAN,TRMIN,TRMAX
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR3 / KEYRAD,KEYOPT
      COMMON / BHPARZ / AMAZ,GAMMZ,SINW2,GV,GA
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      COMMON / INOUT  / NINP,NOUT
      SAVE   / BHPAR2 /, / BHPAR3 /, / BHPARZ /, / TRANSR /
      SAVE   / WGTALL /, / INOUT  /

      KeyZet = MOD(KeyRad,10000)/1000
      s= cmsene**2

!------------------------------------------------------------------
!     Older type (A) with and without exponentiation
!------------------------------------------------------------------
      CALL modl2a(mode)

!------------------------------------------------------------------
!     New  type  (B) with and without exponentiation
!------------------------------------------------------------------
      CALL modl2b(mode)
!
!------------------------------------------------------------------
!     Calculation od Z contribution as in TH.95-74
!------------------------------------------------------------------
      CALL m2agzi(MODE)

!****************************************************************!
!****************************************************************!
!  Model weight, normaly the best... <=====
!****************************************************************!
!  The simple multiplicative ansatz is used to implement         !
!  Vacuum polarization in the LL approximation                   !
!  and s-chanel gamma in the Born approximation                  !
!*****************************************************************
      KeyPia = MOD(KeyRad,10)
      CALL vacpol(KeyPia,-tran,SINW2,RePiE,dRePiE)
!----
      vpfac  = 1/(1 + RePiE)**2
      dels   =  -tran/s +1.5d0*(-tran/s)**2
      sgfac  = 1 + dels
!----
      IF(KeyPia .NE. 0) THEN
!     VacPol and s-chan ON
         IF(KeyZet .EQ. 0) THEN
            WtSet(1) = WtSet(142)*vpfac*sgfac           ! Z off
         ELSEIF(KeyZet .EQ. 1) THEN
            WtSet(1) = WtSet(142)*vpfac*sgfac+WtSet(12) ! Z on
         ELSE
            WRITE(NOUT,*) ' +++++ MODEL: wrong keyZet=',KeyZet
            STOP
         ENDIF
      ELSE
!     VacPol and s-chan OFF
         IF(KeyZet .EQ. 0) THEN
            WtSet(1) = WtSet(142)           ! Z off
         ELSEIF(KeyZet .EQ. 1) THEN
            WtSet(1) = WtSet(142)+WtSet(11) ! Z on
         ELSE
            WRITE(NOUT,*) ' +++++ MODEL: wrong keyZet=',KeyZet
            STOP
         ENDIF
      ENDIF
! This is for special tests with parallel weights
      KeyPia2 = 3
      CALL vacpol(KeyPia2,-tran,SINW2,RePiE2,dRePiE)
      WtSet(2) = 1/(1 + RePiE2)**2
      WtSet(3) = sgfac
!*****************************************************************
      END


      SUBROUTINE MODEL1(MODE)
!     ***********************
!--------------------------------------------------!
!     Version as in CPC 2.01 (1991,92)             !
!     It is kept for backward compatibility tests  !
!     Set KEYMOD=1 if you want to use it           !
!--------------------------------------------------!
! written:     21 jan. 91 (S.Jadach)
! last update:  5 MAY. 91 (S.J.)
! Almost identical to 2.01 version but in addition to
! original one there are two additional delta_Z
! New model weight without the up-down interference
! WT0,1,2 are zero,first,second order exponentiation weights
! The normalisation of distributions is without (1/s)**2 factor
! and without phase space jacobians (like tp*tq/t**2)
! i.e. just matrix element squared
! in beta's S-tilde*d3k/k0 divided out as usual.
!     ************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER( PI = 3.1415926535897932D0, ALFINV=137.03604D0)
      PARAMETER( ALFA= 1D0/ALFINV, ALF1 = 1D0/PI/ALFINV)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR3 / KEYRAD,KEYOPT
      COMMON / TRANSR / TRAN,TRMIN,TRMAX
      COMMON / MOMS1  / TRANP,P1(4),P2(4),PHOT1(100,4),PHSU1(4),NPHOT1
      COMMON / MOMS2  / TRANQ,Q1(4),Q2(4),PHOT2(100,4),PHSU2(4),NPHOT2
      COMMON / MOMZ1  / AL1(50),BE1(50),WTM1(50),MK1(50)
      COMMON / MOMZ2  / AL2(50),BE2(50),WTM2(50),MK2(50)
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      COMMON / BHPARZ / AMAZ,GAMMZ,SINW2,GV,GA
      SAVE   / BHPAR2 /, / BHPAR3 /, / TRANSR /
      SAVE   / MOMS1  /, / MOMS2  /, / MOMZ1  /, / MOMZ2  /
      SAVE   / WGTALL /, / BHPARZ /
      DIMENSION RP(4),RQ(4)
      SAVE ICONT
      DATA ICONT/0/
      SDOT(X1,Y1,Z1,E1,X2,Y2,Z2,E2)
     $            = (E1+E2)**2-(X1+X2)**2-(Y1+Y2)**2-(Z1+Z2)**2
      TDOT(X1,Y1,Z1,E1,X2,Y2,Z2,E2)
     $            = (E1-E2)**2-(X1-X2)**2-(Y1-Y2)**2-(Z1-Z2)**2
! ........................................
      DO 10 K=1,4
      RP(K)= P2(K)+PHSU1(K)
   10 RQ(K)= Q2(K)+PHSU2(K)
      TR = TDOT(P1(1),P1(2),P1(3),P1(4) ,RP(1),RP(2),RP(3),RP(4))
      T  = TDOT(P1(1),P1(2),P1(3),P1(4) ,P2(1),P2(2),P2(3),P2(4))
      T1 = TDOT(Q1(1),Q1(2),Q1(3),Q1(4) ,Q2(1),Q2(2),Q2(3),Q2(4))
      S  = SDOT(P1(1),P1(2),P1(3),P1(4) ,Q1(1),Q1(2),Q1(3),Q1(4))
      S1 = SDOT(P2(1),P2(2),P2(3),P2(4) ,Q2(1),Q2(2),Q2(3),Q2(4))
      U  = TDOT(P1(1),P1(2),P1(3),P1(4) ,Q2(1),Q2(2),Q2(3),Q2(4))
      U1 = TDOT(P2(1),P2(2),P2(3),P2(4) ,Q1(1),Q1(2),Q1(3),Q1(4))
! ........................................
! ..  Crude MC distribution (S-tilde factors omitted)
      CRUDE  = S**2/TR**2
! ........................................
! ..  Contribution from beta0
      BILG   = DLOG(ABS(TR)/AMEL**2)
      BET    = 2D0*ALF1 *( BILG-1D0 )
! .. regular photonic corr. to beta0
      DEL0   = BET
! .. vacuum polarization factor
      VPFAC  = 1D0/(1D0+REPI(TR))**2
      VPFACS = 1D0/(1D0+REPI( S))**2
! .. Z correction
! .. old version, no s-shift
      DELZ =
     $  (GV**2+GA**2)
     $ *(1+TR/S)**3 *2*S**2/(S**2+(S+TR)**2)
     $ *TR*(S-AMAZ**2)/((S-AMAZ**2)**2+(S/AMAZ*GAMMZ)**2)
! .. Z-gamma intrf. has one power of photon and one of Z coup_cons.
     $ *VPFACS/VPFAC
! .. full s-shift (alex read)
      DELZ1 =
     $  (GV**2+GA**2)
     $ *(1+TR/S)**3 *2*S**2/(S**2+(S+TR)**2)
     $ *TR*(S1-AMAZ**2)/((S1-AMAZ**2)**2+(S1/AMAZ*GAMMZ)**2)
     $ *VPFACS/VPFAC
!  It was painfully checked that
!  Chips=0 for initial state and chips=1 for final state!! (s.j.)
!  There are a bit over/under-flows (chips>0, chips<0) at 2% level only.
      CHIPS= 0.5D0
      IF(NPHOT1+NPHOT2.NE.0)
     $  CHIPS = (S**2 *TR**2 - T*T1*S*S1)/(S**2-S1**2)/TR**2
! the cutting below is not really necessary nor relevant numerically
!    CHIPS = MAX(0D0,MIN(1D0,CHIPS))
! half-shift: in the initial state only
      DELZIN = DELZ*CHIPS +DELZ1*(1D0-CHIPS)
! .. s-channel gamma correction
      DELS   =  TR/S +1.5D0*(TR/S)**2
! .......................
      DIS0   = (S**2+U**2+S1**2+U1**2)/4D0 /(T*T1)
      BT00   =  DIS0
      SUM00  =  BT00    /CRUDE
! ..  Leading-Log version ........
      ZETA   = DABS(TRAN)/S
      DISLL0 = 0.5D0*(1 + (1-ZETA)**2)*S**2 /(T*T1)
      BTLL00 =  DISLL0
      SULL00 =  BTLL00    /CRUDE
! ........................................
! ..  Contributions from beta1 upper line
      SULL11  =  0D0
      XPP  = P2(4)*P1(4)-P2(3)*P1(3)-P2(2)*P1(2)-P2(1)*P1(1)
      XPR  = P1(4)*RP(4)-P1(3)*RP(3)-P1(2)*RP(2)-P1(1)*RP(1)
      DELT = AMEL**2/(2*XPP)
      SUM11U=0D0
      DO 150 I=1,NPHOT1
! ..  Numerically safe variables
      A   = AL1(I)
      B   = BE1(I)
      Y   = A +B*DELT
      Z   = B +A*DELT
      P1K = XPR*Z
      P2K = XPR*Y
! ..  soft factor
      WMS0 = A*B/(Y*Z)
      SFC  = 2D0*XPP/(P1K*P2K)*WMS0
! ..  one photon bremss. distribution
      WM = A*B/(Y*Z) +DELT*(Y**2+Z**2)/((1-Y)**2+(1-Z)**2)*(Y/Z+Z/Y)
      DIS1 = (S**2 +U**2 +S1**2+U1**2)
     $       /4D0/(-T1*P1K*P2K)   *WM
! ..  beta1
      BT11U   = DIS1/SFC-BT00
      SUM11U  = SUM11U  + BT11U /CRUDE
! ..  Leading-Log version ........
! Note that in DIS1/SFC we get 1/(T*T1) which after multiplication
! by jacobians (T/TRAN)(T1/TRAN) results in 1/TRAN**2, LL is OK!
      ZETA = DABS(TRAN)/(S*(1-A))
      XX  = MAX(A,B)
      DISLL1 = S**2 *(1 + (1-ZETA)**2) * (1 +(1-XX)**2 )
     $       /4D0/(-T1*P1K*P2K)  *WMS0
      BTLL11   = DISLL1/SFC-BTLL00
      SULL11  = SULL11  + BTLL11 /CRUDE
  150 CONTINUE
! ........................................
! ..  Contributions from beta1 lower line
      XQQ  = Q2(4)*Q1(4)-Q2(3)*Q1(3)-Q2(2)*Q1(2)-Q2(1)*Q1(1)
      XQR  = Q1(4)*RQ(4)-Q1(3)*RQ(3)-Q1(2)*RQ(2)-Q1(1)*RQ(1)
      DELT = AMEL**2/(2*XQQ)
      SUM11L=0D0
      DO 250 I=1,NPHOT2
! ..  numerically safe variables
      A   = AL2(I)
      B   = BE2(I)
      Y   = A+B*DELT
      Z   = B+A*DELT
      Q1K = XQR*Z
      Q2K = XQR*Y
! ..  soft factor
      WMS0 = A*B/(Y*Z)
      SFC  = 2D0*XQQ/(Q1K*Q2K)*WMS0
! ..  one photon bremss. distribution
      WM = A*B/(Y*Z) +DELT*(Y**2+Z**2)/((1-Y)**2+(1-Z)**2)*(Y/Z+Z/Y)
      DIS1 = (S**2 +U**2 +S1**2+U1**2)
     $       /4D0/(-T *Q1K*Q2K)   *WM
! ..  beta1 contrib.
      BT11L   = DIS1/SFC-BT00
      SUM11L  = SUM11L  + BT11L /CRUDE
! ..  Leading-Log version ........
      ZETA = DABS(TRAN)/(S*(1-A))
      XX  = MAX(A,B)
      DISLL1 = S**2 *(1 + (1-ZETA)**2) * (1 +(1-XX)**2 )
     $       /4D0/(-T*Q1K*Q2K)  *WMS0
      BTLL11   = DISLL1/SFC-BTLL00
      SULL11   = SULL11  + BTLL11 /CRUDE
  250 CONTINUE
!-----------------------
! the best and most complete
!-----------------------
! zero order
      WTSET( 11) = SUM00*(1+DELZ+DELS)
! first order
      SUM01      = SUM00*(1D0 + DEL0+DELZ+DELS)
      WTSET( 12) = (SUM01 + SUM11U + SUM11L)*VPFAC
! second order not implemented
      WTSET( 13) = 0D0
!-----------------------
! model weight PURE BREMSS. only, no vac.pol, no Z
! zero order
      WTSET( 51) = SUM00
! first order
      SUM01 = SUM00*(1D0 + DEL0)
      WTSET( 52) = (SUM01 + SUM11U + SUM11L)
! second order not implemented
      WTSET( 53) = 0D0
!-----------------------
! Leading-Log version, pure bremss. only
      WTSET( 61) = SULL00
      SULL01 = SULL00*(1D0 + DEL0)
      WTSET( 62) = (SULL01 + SULL11)
      WTSET( 63) = (0D0            )
!------------------------------------------
! Miscelanous, for tests in TH-6118
! Vacuum polarization effect ONN/OFF
      SUM01      = SUM00*(1D0 + DEL0 )
      WTSET( 20) = (SUM01 + SUM11U + SUM11L)* VPFAC
      WTSET( 21) = (SUM01 + SUM11U + SUM11L)* (VPFAC-1)
! Z-exchange  ON/OFF
      SUM01      = SUM00*(1D0 +DEL0 +DELZ  )
      WTSET( 22) = (SUM01 + SUM11U + SUM11L)*VPFAC
      WTSET( 23) = WTSET(22)-WTSET(20)
! s-schannel exchange
      SUM01      = SUM00*(1D0 +DEL0 +DELZ +DELS )
      WTSET( 24) = (SUM01 + SUM11U + SUM11L)*VPFAC
      WTSET( 25) = WTSET(24)-WTSET(22)
! beta_0,1 contribs.
      WTSET( 26) = (SUM01                  )*VPFAC
      WTSET( 27) = (      + SUM11U + SUM11L)*VPFAC
! new exercises on Z-exchange
      WTSET( 30) = SUM00* DELZ1  *VPFAC
      WTSET( 31) = SUM00* DELZIN *VPFAC
!---------------------------------------------------
! Model weight, normaly the best...
!-----------------------
      WTSET(1) = WTSET(11)
      WTSET(2) = WTSET(12)
      WTSET(3) = WTSET(13)
      KEYPIA = MOD(KEYRAD,10)
! pure bremsstr. as an option:
      IF(KEYPIA.EQ.0) THEN
      WTSET(1) = WTSET(51)
      WTSET(2) = WTSET(52)
      WTSET(3) = WTSET(53)
      ENDIF
      END
