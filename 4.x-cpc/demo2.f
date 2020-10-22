**********************************************************************
*                                                                    *  
*     BBBBBBB    BBB   BBB  BBB      BBB  BBB  BBB     BBB   BBB     *
*     BBB  BBB   BBB   BBB  BBB      BBB  BBB  BBBB   BBBB   BBB     *
*     BBB  BBB   BBB   BBB  BBB      BBB  BBB  BBBBB BBBBB   BBB     *
*     BBBBBB     BBBBBBBBB  BBB      BBB  BBB  BBB BBB BBB   BBB     *
*     BBBBBBBBB  BBBBBBBBB  BBB      BBB  BBB  BBB  B  BBB   BBB     *
*     BBB  BBBB  BBB   BBB  BBB  BB  BBB  BBB  BBB     BBB   BBB     *
*     BBBBBBBBB  BBB   BBB  BBBBBBB  BBB  BBB  BBB     BBB   BBB     *
*     BBBBBBBB   BBB   BBB  BBBBBBB   BBBBBB   BBB     BBB   BBB     *
*                                                                    *
**********************************************************************

! Sophisticated demonstration program

      PROGRAM MAIN
*     **********************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
* histograms in labeled common   
      COMMON / cglib / b(50000)
* Initialization of histograming package --
* here we use double precision HBOOK-like histograming/plotting
* package GLIBK written by S. Jadach, (1990-96),  unpublished.
      CALL glimit(50000) 
!
!
! Sophisticated demonstration program
      CALL Bhldem2

      END


      SUBROUTINE Bhldem2
*     **********************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
* Input/output files
      COMMON / inout  / ninp,nout
      CHARACTER*4 Semaph
      CHARACTER*5 tname,dname
!
!---------------
! General Output for everybody including Glibk
      nout =16
      OPEN(nout,FILE='./bhl.output')
      REWIND(nout)
      CALL goutpu(nout)          
!
      WRITE(nout,*) '   '
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '==========***    Bhldem2    ***==============='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '   '
!---------------
! Standard Input          
      ninp =5
! This OPEN on unix machine is not necessary
*     OPEN(ninp)
! Read test_name and data_set_name         
      READ( ninp,'(A5,1X,A5)') tname,dname
      WRITE(   6,'(A5,1X,A5)') tname,dname
! Read semaphore flag
      CALL GIVSEM(Semaph)

!---------------
      IF(Semaph.eq.'STAR') THEN
         write(6,*) ' ------- Starting from the scratch ----------'
! Read initial (root) random number seed          
         NINP3=3
         OPEN(NINP3,FILE='./iniseed')
         READ(NINP3,'(I10)') ijklin
         READ(NINP3,'(I10)') ntotin
         READ(NINP3,'(I10)') ntot2n
         CALL marini(ijklin,ntotin,ntot2n)
      ELSEIF(Semaph.eq.'CONT') THEN
         write(6,*) ' ------- Restoring from the disk   ----------'
! Restore histograms from the disk            
         ninph=10      
         OPEN(NINPH,file='./'//'bhl.hst')
         CALL GRFILE(ninph,' ',' ')      !Transfer file number
         CALL GRIN(   0,9999,0)          !Read from the disk
         CALL GREND(' ')                 !Close file
! Read random number seed stored in semaphore file          
         NINP2=2
         OPEN(NINP2,FILE='./semaphore')
         READ(NINP2,'(A4)') Semaph
         READ(NINP2,'(I10)') ijklin
         READ(NINP2,'(I10)') ntotin
         READ(NINP2,'(I10)') ntot2n
         CALL MARINI(ijklin,ntotin,ntot2n)
         CLOSE(ninp2)
      ELSEIF(Semaph.eq.'STOP') THEN
         write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(6,*) '++++ STOP: Please Change Semaph to CONT or START !'
         write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         STOP
      ELSE
         write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(6,*) '++++ STOP: Wrong key Semaph = ', Semaph
         write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         STOP
      ENDIF

! *******************************************

            CALL cpcpro

! ********************************************

      END

      SUBROUTINE GIVSEM(Semaph)
C     ************************
      IMPLICIT REAL*8(A-H,O-Z) 
      CHARACTER*4 Semaph
! ------------------------------------------------------
! Read Semaphore flag          
! ------------------------------------------------------
      NINP2=2
      OPEN(NINP2,FILE='./semaphore')
      READ(NINP2,'(A4)') Semaph
      CLOSE(NINP2)
      END


      SUBROUTINE DUMPEH(NEV)
C     ************************
      IMPLICIT REAL*8(A-H,O-Z) 
! ------------------------------------------------------
! Write histos on the disk
! ------------------------------------------------------
      NOUTH=11
      OPEN(NOUTH,FILE='./bhl.hst')
      CALL GRFILE(NOUTH,' ','N')   !Transfer file number
      CALL GROUT( 0,ICY,' ')       !Write on the disk
      CALL GREND(' ')              !Close file
! ------------------------------------------------------
! Overwrite Semaphore file flag
      NINP2=2
      OPEN(NINP2,FILE='./semaphore')
      WRITE(NINP2,'(A4)') 'CONT'
! Append semaphore file with new random number seed in          
      CALL MAROUT(IJKLIN,NTOTIN,NTOT2N)
      WRITE(NINP2,'(I10,A)') IJKLIN, ' = IJKLIN '
      WRITE(NINP2,'(I10,A)') NTOTIN, ' = NTOTIN '
      WRITE(NINP2,'(I10,A)') NTOT2N, ' = NTOT2N '
      WRITE(NINP2,'(I10,A)') NEV,    ' =    NEV '
      CLOSE(NINP2)
! ------------------------------------------------------
      write( 6,*) ' DUMPEH: Histos Dumped into disk, NEV= ', NEV
      CALL DUMPS( 6)
      END

      SUBROUTINE cpcpro
!     *****************
! -------------------------------------------------------------------
!  Trigger types: ALEPH/Sical and OPAL/SiW
! -------------------------------------------------------------------
!     ************************
      IMPLICIT REAL*8(A-H,O-Z) 
      PARAMETER( PI = 3.1415926535897932D0 )
      SAVE
      CHARACTER*4 Semaph
      COMMON / INOUT  / NINP,NOUT 
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI
      DIMENSION NPAR(100),XPAR(100)
!
      WRITE(NOUT,*) '   '
      WRITE(NOUT,*) '========================================='
      WRITE(NOUT,*) '==========****************==============='
      WRITE(NOUT,*) '==========***  cpcpro  ***==============='
      WRITE(NOUT,*) '==========****************==============='
      WRITE(NOUT,*) '========================================='
      WRITE(NOUT,*) '   '
!=======================================================
      READ( NINP,'(8I2)')   KAT1,KAT2,KAT3,KAT4,KAT5,KAT6
      READ( NINP,'(I10)')   NEVT,KEYOPT,KEYRAD,KEYTRI
      READ( NINP,'(F10.0)') CMSENE,TMING,TMAXG,VMAXG,XK0
      READ( NINP,'(F10.0)') TMINW,TMAXW,TMINN,TMAXN,VMAXE  
      READ( NINP,'(I10)')   NPHI,NTHE
!=======================================================
! control output
      WRITE(NOUT,'(6A6/6I6)')
     $ 'KAT1','KAT2','KAT3','KAT4','KAT5','KAT6',
     $  KAT1 , KAT2 , KAT3 , KAT4 , KAT5 , KAT6
      WRITE(NOUT,'(4A12/4I12)') 
     $  'NEVT','KEYRAD','KEYOPT','KEYTRI',
     $   NEVT,  KEYRAD , KEYOPT , KEYTRI
      WRITE(NOUT,'(5A12/5F12.6)')
     $ 'CMSENE','TMING','TMAXG','VMAXG','XK0',
     $  CMSENE , TMING , TMAXG , VMAXG , XK0
      WRITE(NOUT,'(5A12/5F12.6)')
     $  'TMINW','TMAXW','TMINN','TMAXN','VMAXE',
     $   TMINW , TMAXW , TMINN , TMAXN , VMAXE
      WRITE(NOUT,'(6A12/6F12.6)')
     $  'TMING','TMAXG','TMINW','TMAXW','TMINN','TMAXN',
     $   TMING , TMAXG , TMINW , TMAXW , TMINN , TMAXN   
      WRITE(NOUT,'(2A12/2I12)') 
     $  'NPHI','NTHE',
     $   NPHI,  NTHE
!=======================================================
      KEYGEN = MOD(KEYOPT,10000)/1000
      IF(KEYGEN.EQ.3) THEN
! input data for --- BHLUM2 ---
! Born limiting values for the transfer
       TRMINB =CMSENE**2*(1D0-COS(TMING))/2D0
       TRMAXB =CMSENE**2*(1D0-COS(TMAXG))/2D0
       TRMIN = TRMINB
       TRMAX = TRMAXB
       EPSCM = XK0
       NPAR(1)=KEYOPT
       NPAR(2)=KEYRAD
       XPAR(1)=CMSENE
       XPAR(2)=TRMIN
       XPAR(3)=TRMAX
       XPAR(4)=EPSCM
      ELSEIF(KEYGEN.EQ.2) THEN
! input data for --- LUMLOG ---
       NPAR(1) = KEYOPT
       NPAR(2) = KEYRAD
       XPAR(1) = CMSENE
       XPAR(2) = TMING*180/PI
       XPAR(3) = TMAXG*180/PI
       XPAR(4) = XK0
       XPAR(5) = VMAXG
      ELSEIF(KEYGEN.EQ.1) THEN 
! input data for --- OLDBIS ---
       NPAR(1)= KEYOPT
       NPAR(2)= KEYRAD
       XPAR(1)= CMSENE
       XPAR(2)= TMING*180/PI
       XPAR(3)= TMAXG*180/PI
       XPAR(4)= XK0
       XPAR(5)= VMAXG
       XPAR(6)= 0D0
      ELSE
       WRITE(6,*) '++++ WRONG KEYOPT=',KEYOPT
      STOP
      ENDIF
!-------------------------------------------------------!
!                 Initialization                        !
!-------------------------------------------------------!
      CALL BHLUMI(  -1,XPAR,NPAR)  
*     ================================== 
      IF(KAT1.EQ.1) CALL ROBOL1(-1)      
      IF(KAT2.EQ.1) CALL ROBOL2(-1)      

      WRITE(6,'(F10.2,A)') NEVT/1.E6,' Mega-events requested'
      WRITE(6,*)  ' =======> Generation starts...'
!-------------------------------------------------------!
!                 Main MC loop                          !
!-------------------------------------------------------!
      NGROUP = 100000
      IEV=0  
      DO LOOP=1,10000000 
        DO IGROUP =1,NGROUP
          IEV=IEV+1 
          IF(MOD(IEV, NGROUP).EQ.1) WRITE( 6,*)  'IEV= ',IEV 
          CALL BHLUMI(   0,XPAR,NPAR)  
          IF(IEV.LE. 4) CALL DUMPS( 6)       
          IF(IEV.LE. 4) CALL DUMPS(16) 
!         ============================================
!         Histograming 
          IF(KAT1.EQ.1) CALL ROBOL1( 0) 
          IF(KAT2.EQ.1) CALL ROBOL2( 0) 
!         ============================================
!         Check on equested no. of events 
          IF(IEV.EQ.NEVT)     GOTO 300 
        ENDDO
!       Check on semaphore flag
        CALL GIVSEM(Semaph)
        IF(Semaph.EQ.'STOP') GOTO 300
!       Dump partial results on the disk after every NGROUP
        CALL  DUMPEH(IEV)
      ENDDO
 300  CONTINUE   
      CALL  DUMPEH(IEV)
!-------------------------------------------------------!
!                 Final activity                        !
!-------------------------------------------------------!
      IF(KAT1.EQ.1) CALL ROBOL1( 1)  
      IF(KAT2.EQ.1) CALL ROBOL2( 1)  
      END


      SUBROUTINE robol1(mode)   
!     ***********************
!--------------------------------------------------------!
! This is for calculation of the Z contribution (prod1)  !
!--------------------------------------------------------!
*     ***********************    
      IMPLICIT REAL*8(A-H,O-Z)   
      SAVE
      COMMON / INOUT  / NINP,NOUT 
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300) 
      DIMENSION NPAR(100),XPAR(100)
      DIMENSION Zsic1(4),Zsic2(4),Zsic3(4),Zsic4(4)

      IF(mode.eq.-1) THEN  
*     *******************
!=============================================================
!.... initialize histos
      KEYGEN = MOD(KEYOPT,10000)/1000
      JDA = 10000*KEYGEN
      NBIV =    1
      VMAX =  1d0 -0.43d0
      VMIN =  0D0
      DO K=1,3
       CALL GBOOK1(JDA+100+K,'Z-Born      $',NBIV,VMIN,VMAX)
       CALL GBOOK1(JDA+200+K,'Z O(alf)    $',NBIV,VMIN,VMAX)
       CALL GBOOK1(JDA+300+K,'Z O(alf)exp $',NBIV,VMIN,VMAX)
      ENDDO
      NEVGEN=0 

      ELSEIF(MODE.EQ.0) THEN
*     **********************

      NEVGEN=NEVGEN+1 
      wtcrud = wtcru1*wtcru2 

      DO i=1,4
       Zsic1(i)=5
      ENDDO
      IF(wtcrud .NE. 0D0 ) THEN
!---------------------------------------------------------------
! Aleph Sical trigger general input:
        THsic1 =.024d0     ! Tmin (rad)  theta_min sical
        THsic2 =.058d0     ! Tmax (rad)  theta_max sical
        Nthe = 16
        Nphi = 32
        CALL TRISIC2(Nphi,Nthe,THsic1,THsic2,Zsic1,Zsic2,Zsic3)

! O(alf2)eB, Energy cut versus Asimetricity
        CALL Filix(JDA+100,Zsic1,wtcrud*wtset(80) ) ! Z-Born 
        CALL Filix(JDA+200,Zsic1,wtcrud*wtset(82) ) ! Z O(alf)
        CALL Filix(JDA+300,Zsic1,wtcrud*wtset(12) ) ! Z O(alf)exp
      ENDIF
!---------------------------------------------------------------
      ELSEIF(MODE.EQ.1) THEN 
*     ***********************

C control printout
      CALL gminim(0,-1d0)
      DO i=1,3
       CALL gidopt(JDA+100+i,'ERRO')
       CALL gidopt(JDA+200+i,'ERRO')
       CALL gidopt(JDA+300+i,'ERRO')
       CALL GPRINT(JDA+100+i)
       CALL GPRINT(JDA+200+i)
       CALL GPRINT(JDA+300+i)
      ENDDO


      ENDIF      
*     *****

      END 


      SUBROUTINE ROBOL2(MODE)
!     ***********************
!--------------------------------------------------------!
!   SMALLER angles and NEWER triggers                    !
!   New Triggers:  SICAL, OSiW                           !
!   Any of three sub-generator BHLUM2, OLDBIS, LUMLOG.   !
!--------------------------------------------------------!
! This for any type of generator BHLUM2, OLDBIS, LUMLOG.
! Using these results, in the later analysis,
! the inter-generator differences are calculated.
*     ***********************
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      PARAMETER( PI = 3.1415926535897932D0 )
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
      COMMON / INOUT  / NINP,NOUT
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PHOT(100,4),NPHOT
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      DIMENSION NPAR(100),XPAR(100)
!------------------------------------------------------------
      DIMENSION Zsic1(4),Zsic2(4),Zsic3(4),Zsic4(4)
      DIMENSION Zosi1(4),Zosi2(4),Zosi3(4),Zosi4(4)
!------------------------------------------------------------
      DIMENSION z1A(4),z1B(4),z1C(4)
      DIMENSION y1A(4),y1B(4),y1C(4)
      DIMENSION u1A(4),u1B(4),u1C(4)
      DIMENSION z2A(4),z2B(4),z2C(4)
      DIMENSION z3A(4),z3B(4),z3C(4)
!------------------------------------------------------------

      IF(MODE.EQ.-1) THEN
*     *******************
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) 'ROBOL2 start initialisation....  '
      WRITE(   6,BXTXT) 'ROBOL2 start initialisation....  '
!

      TH1W=TMINW
      TH2W=TMAXW
      TH1N=TMINN
      TH2N=TMAXN
      IF(KEYTRI.eq.2) THEN
        PAD = (TMAXW-TMINW)/NTHE
        BORWID  = BORNB(CMSENE,TMINW   +PAD, TMAXW   -PAD)
        BORNAR  = BORNB(CMSENE,TMINW +2*PAD, TMAXW -4*PAD)
      ELSE
        BORWID  = BORNB(CMSENE,TH1W,TH2W)
        BORNAR  = BORNB(CMSENE,TH1N,TH2N)
      ENDIF
!.... initialize histos
      NV =   40
      VMA =  1D0
      VMI =  0D0
      AcMA =  0.025D0
      AcMI =      0D0
!....   O(alf1)exp
      KEYGEN = MOD(KEYOPT,10000)/1000
!=============================================================
      JDA = 10000*KEYGEN +2000
      JDU = 10000*KEYGEN +4000
! OPSiW 3-generators
      DO K=1,3
         CALL GBOOK1(JDA+200+K,'B4-(O1+L3) Ener/Asymetri$',NV,VMI,VMA)
         CALL GBOOK1(JDA+210+K,'B4-(O1+L3) Ener/Acolline$',NV,VMI,VMA)
         CALL GBOOK1(JDA+220+K,'B4-(O1+L3) Ener/Asymetri$',NV,VMI,VMA)
!     Acollinearity
         CALL GBOOK1(JDA+230+K,'B4-(O1+L3) Acoll/Asymetr$',NV,AcMI,AcMA)
      ENDDO
! SICAL 3-generators
      DO K=1,3
         CALL GBOOK1(JDA+300+K,'B4-(O1+L3) Ener/Asymetri$',NV,VMI,VMA)
         CALL GBOOK1(JDA+310+K,'B4-(O1+L3) EnerTh/Asymet$',NV,VMI,VMA)
         CALL GBOOK1(JDA+320+K,'B4-(O1+L3) Ener/ClusSize$',NV,VMI,VMA)
      ENDDO
! BARE1, CALO2, SICAL2 
      DO K=1,3
         CALL GBOOK1(JDA+340+K,'B4-(O1+L3) Ener/Asymetri$',NV,VMI,VMA)
         CALL GBOOK1(JDA+350+K,'B4-(O1+L3) Ener/Asymetri$',NV,VMI,VMA)
         CALL GBOOK1(JDA+360+K,'B4-(O1+L3) Ener/Asymetri$',NV,VMI,VMA)
      ENDDO
! SICAL, CALO2, SICAL2 Others
      DO K=1,3
         CALL GBOOK1(JDA+370+K,'O(alf2eVPZ) Ener/Asymetri$',NV,VMI,VMA)
         CALL GBOOK1(JDA+380+K,'O(alf2eVPZ) Ener/Asymetri$',NV,VMI,VMA)
         CALL GBOOK1(JDA+390+K,'O(alf2eVPZ) Ener/Asymetri$',NV,VMI,VMA)
         CALL GBOOK1(JDA+400+K,'O(alf2eZ)   Ener/Asymetri$',NV,VMI,VMA)
         CALL GBOOK1(JDA+410+K,'O(alf2eZ)   Ener/Asymetri$',NV,VMI,VMA)
      ENDDO
! Differences
      IF(KeyGen .EQ. 2 .OR. KeyGen .EQ. 3) THEN
      DO K=1,3
         CALL GBOOK1(JDU+200+K,'B4-B2  Ener/Asymetr OPSiW$',NV,VMI,VMA)
         CALL GBOOK1(JDU+210+K,'LUMLOG                   $',NV,VMI,VMA)
         CALL GBOOK1(JDU+220+K,'LUMLOG                   $',NV,VMI,VMA)
         CALL GBOOK1(JDU+230+K,'LUMLOG     Acoll/Asymetr$',NV,AcMI,AcMA)
         CALL GBOOK1(JDU+260+K,'B4-B2  Ener/Acoll   OPSiW$',NV,VMI,VMA)
         CALL GBOOK1(JDU+270+K,'B4-B2  Acol/Asymet  OPSiW$',NV,VMI,VMA)
         CALL GBOOK1(JDU+290+K,'Miss. O(alf3) LUMLG OPSiW$',NV,VMI,VMA)
!
         CALL GBOOK1(JDU+300+K,'B4-B2  Ener/Asymet SICAL $',NV,VMI,VMA)
         CALL GBOOK1(JDU+310+K,'O(alf2)eB-O(alf2)B SICAL $',NV,VMI,VMA)
         CALL GBOOK1(JDU+320+K,'O(alf2)eB-O(alf2)B SICAL2$',NV,VMI,VMA)
         CALL GBOOK1(JDU+330+K,'O(alf2)eB-O(alf2)B CALO2 $',NV,VMI,VMA)
         CALL GBOOK1(JDU+340+K,'O(alf2)eB-O(alf2)B BARE1 $',NV,VMI,VMA)

         CALL GBOOK1(JDU+350+K,'B4-B2  Ener/Asymet SICAL $',NV,VMI,VMA)
         CALL GBOOK1(JDU+360+K,'B4-B2  Ener/Asymet SICAL2$',NV,VMI,VMA)
         CALL GBOOK1(JDU+370+K,'B4-B2  Ener/Asymet CALO2 $',NV,VMI,VMA)
         CALL GBOOK1(JDU+380+K,'B4-B2  Ener/Asymet BARE1 $',NV,VMI,VMA)

         CALL GBOOK1(JDU+400+K,'B4-B2 Ener/ClusSize SICAL$',NV,VMI,VMA)

         CALL GBOOK1(JDU+410+K,'O(alf2)eB-O(alf2)eA SICAL2$',NV,VMI,VMA)
         CALL GBOOK1(JDU+420+K,'O(alf2)B-O(alf2)A   SICAL2$',NV,VMI,VMA)
         CALL GBOOK1(JDU+430+K,'O(alf2)eB-O(alf1)eB SICAL2$',NV,VMI,VMA)
         CALL GBOOK1(JDU+440+K,'O(alf2)B -O(alf1)B  SICAL2$',NV,VMI,VMA)
         CALL GBOOK1(JDU+450+K,'O(alf3e)B-O(alf2)   SICAL $',NV,VMI,VMA)
         CALL GBOOK1(JDU+460+K,'O(alf3e)B-O(alf2)   BARE1 $',NV,VMI,VMA)
      ENDDO
      ENDIF
      IF(KeyGen .EQ. 3) THEN
! Correction due do LL bug
      DO K=1,3
         CALL GBOOK1(JDU+600+K,'LLbug  Ener/Asymetr OPSiW$',NV,VMI,VMA)
         CALL GBOOK1(JDU+610+K,'LLbug  EneTh/Asymet OPSiW$',NV,VMI,VMA)
         CALL GBOOK1(JDU+700+K,'LLbug  Ener/Asymetr SICAL$',NV,VMI,VMA)
         CALL GBOOK1(JDU+710+K,'LLbug  EneTh/Asymet SICAL$',NV,VMI,VMA)
      ENDDO
      ENDIF
!=============================================================
      NEVGEN=0
      trminf = CMSENE**2* (1D0-COS(tminw))/2D0
      trmaxf = CMSENE**2* (1D0-COS(tmaxw))/2D0
      WRITE(NOUT,BXL1I) KEYTRI,     'Type of Trigger    ','KEYTRI','  '
      WRITE(NOUT,BXL1F) tminw ,     'thet_min     [rad] ','tminw ','  '
      WRITE(NOUT,BXL1F) tmaxw ,     'thet_max     [rad] ','tmaxw ','  '
      WRITE(NOUT,BXL1F) trminf,     't_min     [GeV**2] ','trminf','  '
      WRITE(NOUT,BXL1F) trmaxf,     't_max     [GeV**2] ','trmaxf','  '
      WRITE(NOUT,BXL1I) NPHI,       'No of phi-sectors  ','NPHI  ','  '
      WRITE(NOUT,BXL1I) NTHE,       'No of phi-sectors  ','NTHE  ','  '
      WRITE(NOUT,BXL1F) BORWID,     'Born Wide    [nb]  ','BORWID','  '
      WRITE(NOUT,BXL1F) BORNAR,     'Born Narrow  [nb]  ','BORNAR','  '
      WRITE(NOUT,BXTXT) '  .... Robol1 end initialization '
      WRITE(NOUT,BXCLO)

      ELSEIF(MODE.EQ.0) THEN
*     **********************
      NEVGEN=NEVGEN+1
      WTCRUD = WTCRU1*WTCRU2

      IF(WTCRUD.NE.0D0 ) THEN
!     ***********************
!
!---------------------------------------------------------------
! Event Selections
!---------------------------------------------------------------
! Fidutial theta size like in SICAL
        th1f=0.024d0
        th2f=0.058d0
        ntheta=16
        nphi  =32
        padthe = (th2f-th1f)/ntheta
        padphi = 2d0*pi/nphi

! ***BARE1***, 
! here z1A is for BARE1
! ***CALO1***, calorimetric cones 10mrads,
! here z1B is for CALO1, calorimeter delta-cone
        th1w=th1f
        th2w=th2f
        th1n=th1f   +padthe
        th2n=th2f   -padthe
        dcone = 0.010d0 ! radius of the cone
        dlph  = 0.010d0 ! UNUSED half-size of theta window
        dlth  = 0.010d0 ! UNUSED half-size of phi   window
        CALL trical1(th1n,th2n,th1w,th2w,dcone,dlth,dlph,z1A,z1B,z1C)
!
! ***CALO2***, (theta*phi) rectangles (semi-rings) like in SICAL,
! here y1C is for CALO2 with theta*phi plaquettes
        th1w=th1f   +padthe
        th2w=th2f   -padthe
        th1n=th1f +2*padthe
        th2n=th2f -4*padthe
        dcone = 0.010d0         ! UNUSED radius of the cone
        dlph  = 3d0/2d0 *padphi ! half-size of theta window
        dlth  = 3d0/2d0 *padthe ! half-size of phi   window
        CALL trical1(th1n,th2n,th1w,th2w,dcone,dlth,dlph,y1A,y1B,y1C)

! ***SICAL2***, variant with 3*3 pads as in CALO2, as in WORKSHOP95
! Only z3A is defined (as SICAL2)
        Npad =1
        Nseg =1
        CALL trisic2W(Nphi,Ntheta,Th1f,Th2f,Nseg,Npad,z3A,z3B,z3C)

! ALEPH ****SICAL**** trigger general input:
        THsic1 =.024d0     ! Tmin (rad)  theta_min sical
        THsic2 =.058d0     ! Tmax (rad)  theta_max sical
        Nthe = 16
        Nphi = 32
        CALL TRISIC2(Nphi,Nthe,THsic1,THsic2,Zsic1,Zsic2,Zsic3)

! OPAL ****OSiW**** trigger, general input:
        THosi1 = 0.025204199d0
        THosi2 = 0.057721555d0
        Nthe = 32
        Nphi = 32
        CALL TRIOSiW(Nphi,Nthe,THosi1,THosi2,Zosi1,Zosi2,Zosi3,Zosi4)
!       *************************************************************

!---------------------------------------------------------------
!                  Bhlum2 specific
!---------------------------------------------------------------
      IF(KeyGen.eq.3) THEN
! Absolute
! **********************************************************
! *********** 3-generator class ****************************
! =======  true Event Selections SICAL and OSiW   ==========
!                     OSiW    TRUE exper.
! O(alf2)eB, Energy cut versus Asimetrycity
        CALL Filix(JDA+200,Zosi1,wtcrud*wtset(142))  ! OSiW
! O(alf2)eB, Theorist energy cut versus Asimetricity
        CALL Filix(JDA+210,Zosi3,wtcrud*wtset(142))  ! OSiW
! O(alf2)eB, Energy cut versus Acollinearity
        CALL Fili3(JDA+220,Zosi2,wtcrud*wtset(142))  ! OSiW
! O(alf2)eB, Acollinearity versus Asimetricity
        CALL Filix(JDA+230,Zosi4,wtcrud*wtset(142))  ! OSiW
!                     SICAL   TRUE exper.
! O(alf2)eB, Energy cut versus Asimetrycity
        CALL Filix(JDA+300,Zsic1,wtcrud*wtset(142))  !#SICAL# P.L.
! O(alf2)eB, Theorist energy cut versus Asimetricity
        CALL Filix(JDA+310,Zsic3,wtcrud*wtset(142))  ! SICAL Glasgow
! O(alf2)eB, Energy cut versus Cluster size
        CALL Fili3(JDA+320,Zsic2,wtcrud*wtset(142))  ! SICAL
! =======  Simplified Event Selections as in Wshop95 =======
! O(alf2)eB, Energy cut versus Asimetrycity
        CALL Filiv(jda+340,z1A,wtcrud*wtset(142)) ! BARE1
        CALL Filiv(jda+350,y1C,wtcrud*wtset(142)) ! CALO2
        CALL Filiv(jda+360,z3A,wtcrud*wtset(142)) ! SICAL2
! ********* 3-generator class THE END  *********************
! ----- Other ----- 
! O(alf2)eB,    VP+Z included, Z included
        wtVPZ = wtcrud*( wtset(142)*wtset(2)*wtset(3) +wtset(12))
        CALL Filix(JDA+370,Zsic1,wtVPZ)   ! SICAL
        CALL Filiv(jda+380,  y1C,wtVPZ)   ! CALO2
        CALL Filiv(jda+390,  z3A,wtVPZ)   ! SICAL2
        wtZ = wtcrud*( wtset(142) +wtset(11))
        CALL Filiv(jda+400,  y1C,wtZ)     ! CALO2
        CALL Filiv(jda+410,  z3A,wtZ)     ! SICAL2
! **********************************************************
! ------ Matrix Element Differences  --------------
! **********************************************************
! BHLUMI  4.03 - 2.02  difference [Bhlumi4-Bhlumi2 'Correction']
        WtCorr= WtCrud*(wtset(142)-wtset( 41))
        CALL Filix(JDU+300,Zsic1,WtCorr) !#SICAL# true experim.
        CALL Filix(JDU+200,Zosi1,WtCorr) ! OSiW   true experim.
! Specials on CUT-OFF variations for 
! Energy cut versus Acollinearity
        CALL Fili3(JDU+260,Zosi2,WtCorr) ! OSiW  true experim.
! Acollinearity versus Asimetricity
        CALL Filix(JDU+270,Zosi4,WtCorr) ! OSiW  true experim.
! Energy cut versus Cluster size
        CALL Fili3(JDU+400,Zsic2,WtCorr) ! SICAL true experim.
! END OF [Bhlumi4-Bhlumi2 'Correction']
! ------
!  B Switching-off exponentiation O(alf2)eB-O(alf2)B (Phys. Lett. 95)
        WtUEX =wtcrud*(wtset(172)-wtset(142))
        CALL Filix(JDU+310,Zsic1,WtUEX) !#SICAL#
!  B Switching-off exponentiation   O(alf2)eB-O(alf2)B
        CALL Filiv(JDU+320,  z3A,WtUEX) ! SICAL2
        CALL Filiv(JDU+330,  y1C,WtUEX) ! CALO2
        CALL Filiv(JDU+340,  z1A,WtUEX) ! BARE1
! ------
! BHLUMI  4.03-2.02 difference  [Bhlumi4-Bhlumi2 'Correction']
        WtCorr= WtCrud*(wtset(142)-wtset( 41))
        CALL Filix(JDU+350,Zsic1,WtCorr) !#SICAL#
        CALL Filiv(JDU+360,  z3A,WtCorr) ! SICAL2
        CALL Filiv(JDU+370,  y1C,WtCorr) ! CALO2
        CALL Filiv(JDU+380,  z1A,WtCorr) ! BARE1
! ------
!  (B)-(A)  EXP   Second Order      O(alf2)eB-O(alf2)eA
        CALL Filiv(JDU+410,z3A,wtcrud*(wtset(142)-wtset( 42))) ! SICAL2
!  (B)-(A)  UNEXP Second Order      O(alf2)B-O(alf2)A
        CALL Filiv(JDU+420,z3A,wtcrud*(wtset(172)-wtset( 72))) ! SICAL2
! ------
!  (B) Second-First Order  EXP      O(alf2)eB-O(alf1)eB
        CALL Filiv(JDU+430,z3A,wtcrud*(wtset(142)-wtset(141))) ! SICAL2
!  (B) Second-First Order  UNEXP    O(alf2)B -O(alf1)B
        CALL Filiv(JDU+440,z3A,wtcrud*(wtset( 72)-wtset( 71))) ! SICAL2
! ------
! Influence of LL BUG, correction with respect to version WITH bug
        CALL Filix(JDU+600,Zosi1,wtcrud*(wtset(142)-wtset(242))) ! OSiW
        CALL Filix(JDU+700,Zsic1,wtcrud*(wtset(142)-wtset(242))) ! SICAL
! The same but Theorist's energy cut
        CALL Filix(JDU+610,Zosi3,wtcrud*(wtset(142)-wtset(242))) ! OSiW
        CALL Filix(JDU+710,Zsic3,wtcrud*(wtset(142)-wtset(242))) ! SICAL

      ENDIF
!---------------------------------------------------------------
!                   Oldbis specific
!---------------------------------------------------------------
      IF(KeyGen.eq.1) THEN
! **********************************************************************
! *******  3-generator class difference B4-(O1+L3) *********************
! O(alf1) Energy cut versus Asymetricity
        CALL Filix(JDA+200,Zosi1,WTMOD) ! OSiW
! Theorist energy cut versus Asimetricity
        CALL Filix(JDA+210,Zosi3,WTMOD) ! OSiW
! O(alf1) Energy cut versus Acolineraity
        CALL Fili3(JDA+220,Zosi2,WTMOD) ! OSiW
! O(alf1) Acolineraity versus Asymetricity
        CALL Filix(JDA+230,Zosi4,WTMOD) ! OSiW
!
! O(alf1) Energy cut versus Asymetricity
        CALL Filix(JDA+300,Zsic1,WTMOD) ! SICAL
! Theorist energy cut versus Asimetricity
        CALL Filix(JDA+310,Zsic3,WTMOD) ! SICAL
! O(alf1) Energy cut versus Cluster size
        CALL Fili3(JDA+320,Zsic2,WTMOD) ! SICAL
! =======  Simplified Event Selections as in Wshop95 =======
! O(alf2)eB, Energy cut versus Asimetrycity
        CALL Filiv(jda+340,z1A,WtMod) ! BARE1
        CALL Filiv(jda+350,y1C,WtMod) ! CALO2
        CALL Filiv(jda+360,z3A,WtMod) ! SICAL2
! ******  3-generator class difference THE END    **********************
! **********************************************************************
!
! *******
      ENDIF

!---------------------------------------------------------------
!                     LUMLOG specific
!---------------------------------------------------------------
      IF(KeyGen.eq.2) THEN
! **********************************************************************
! ****************  3-generator class **********************************
!                  O(alf3)exp -O(alf1) 
         WtDifz = WtCrud*(wtset(4)-wtset(12))
! ----------------------------------------------------------
! Sical energy cut versus Asimetricity
        CALL Filix(jdu+200,Zosi1,WtDifz)   ! OSiW
! Theorist energy cut versus Asymetricity
        CALL Filix(jdu+210,Zosi3,WtDifz)   ! OSiW
! Energy cut versus Acolineraity
        CALL Fili3(jdu+220,Zosi2,WtDifz)   ! OSiW
! Acolineraity versus Asimetricity
        CALL Filix(jdu+230,Zosi4,WtDifz)   ! OSiW
! ----------------------------------------------------------
! Sical energy cut versus Asimetricity
        CALL Filix(jdu+300,Zsic1,WtDifz)   ! SICAL
! Theorist energy cut versus Asymetricity
        CALL Filix(jdu+310,Zsic3,WtDifz)   ! SICAL
! Energy cut versus Cluster size
        CALL Fili3(jdu+320,Zsic2,WtDifz)   ! SICAL
! =======  Simplified Event Selections as in Wshop95 =======
! O(alf2)eB, Energy cut versus Asimetrycity
        CALL Filiv(jdu+340,  z1A,WtDifz)   ! BARE1
        CALL Filiv(jdu+350,  y1C,WtDifz)   ! CALO2
        CALL Filiv(jdu+360,  z3A,WtDifz)   ! SICAL2
!=======================================================================
!                     O(alf1) for SABSPV
        WtFirs = WtCrud*wtset(12)
! ----------------------------------------------------------
! Sical energy cut versus Asimetricity
        CALL Filix(jda+200,Zosi1,WtFirs)   ! OSiW
! Theorist energy cut versus Asymetricity
        CALL Filix(jda+210,Zosi3,WtFirs)   ! OSiW
! Energy cut versus Acolineraity
        CALL Fili3(jda+220,Zosi2,WtFirs)   ! OSiW
! Acolineraity versus Asimetricity
        CALL Filix(jda+230,Zosi4,WtFirs)   ! OSiW
! ----------------------------------------------------------
! Sical energy cut versus Asimetricity
        CALL Filix(jda+300,Zsic1,WtFirs)   ! SICAL
! Theorist energy cut versus Asymetricity
        CALL Filix(jda+310,Zsic3,WtFirs)   ! SICAL
! Energy cut versus Cluster size
        CALL Fili3(jda+320,Zsic2,WtFirs)   ! SICAL
! =======  Simplified Event Selections as in Wshop95 =======
! O(alf2)eB, Energy cut versus Asimetrycity
        CALL Filiv(jda+340,  z1A,WtFirs)   ! BARE1
        CALL Filiv(jda+350,  y1C,WtFirs)   ! CALO2
        CALL Filiv(jda+360,  z3A,WtFirs)   ! SICAL2
! **************  THE END 3-generator class    *************************
! **********************************************************************

! --- O(alf2)-O(alf1) UNEXP
        CALL Filix(JDA+400,Zsic3,wtcrud*(wtset(13)-wtset(12)))! SICAL
! Missing third order
        WtMis3 = wtcrud*(wtset(4)-wtset( 42))
        CALL Filiv(JDU+410,  z3A,WtMis3)   ! SICAL2
        CALL Filiv(JDU+420,  y1C,WtMis3)   ! CALO2
        CALL Filiv(JDU+430,  z1A,WtMis3)   ! BARE1
        CALL Filix(JDU+440,Zsic3,WtMis3 )  ! SICAL
        WtMis3u = wtcrud*(wtset(4)-wtset( 13))
        CALL Filix(JDU+450,Zsic3,WtMis3u)  ! SICAL
        CALL Filiv(JDU+460,  z1A,WtMis3u)  ! BARE1
!---------------------------------------------------------------
      ENDIF
!---------------------------------------------------------------
      ENDIF ! (WTCRUD.NE.0D0 )

      ELSEIF(MODE.EQ.1) THEN
!     ***********************
      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '============ ROBOL1 ============='
      CALL BHLUMI(   2,XPAR,NPAR)
!---------------------------------------------------------!
! control printout
      DO K=1,3
         CALL gprint(JDA+200+K)
         CALL gprint(JDA+220+K)
         CALL gprint(JDA+230+K)
!
         CALL gprint(JDA+300+K)
         CALL gprint(JDA+310+K)
         CALL gprint(JDA+320+K)

         CALL gprint(JDA+340+K)
         CALL gprint(JDA+350+K)
         CALL gprint(JDA+360+K)

         CALL gprint(JDA+370+K)
         CALL gprint(JDA+380+K)
         CALL gprint(JDA+390+K)
!
         CALL gprint(JDU+200+K)
         CALL gprint(JDU+260+K)
         CALL gprint(JDU+270+K)
         CALL gprint(JDU+290+K)

         CALL gprint(JDU+300+K)
         CALL gprint(JDU+310+K)
         CALL gprint(JDU+320+K)
         CALL gprint(JDU+330+K)
         CALL gprint(JDU+340+K)

         CALL gprint(JDU+350+K)
         CALL gprint(JDU+360+K)
         CALL gprint(JDU+370+K)
         CALL gprint(JDU+380+K)

         CALL gprint(JDU+400+K)
      ENDDO
      DO K=1,3
         CALL gprint(JDU+410+K)
         CALL gprint(JDU+420+K)
         CALL gprint(JDU+430+K)
         CALL gprint(JDU+440+K)
         CALL gprint(JDU+450+K)
         CALL gprint(JDU+460+K)
      ENDDO
      DO K=1,3
         CALL gprint(JDU+600+K)
         CALL gprint(JDU+610+K)
         CALL gprint(JDU+700+K)
         CALL gprint(JDU+710+K)
      ENDDO

      ENDIF
!     *****
      END


      SUBROUTINE FILIX(JD,Z,WTMOD)
!     ****************************
      IMPLICIT REAL*8(A-H,O-Z)
! Fill WITH symetrisation NW<=>WN
      DIMENSION Z(*)
!
      CALL GF1(JD+1, Z(1) ,WTMOD)
      CALL GF1(JD+2, Z(2) ,WTMOD)
      CALL GF1(JD+3, Z(3) ,WTMOD/2D0)
      CALL GF1(JD+3, Z(4) ,WTMOD/2D0)
      END 

      SUBROUTINE FILI3(JD,Z,WTMOD)
!     ****************************
      IMPLICIT REAL*8(A-H,O-Z)
! Fill WITHOUT symetrisation NW<=>WN
      DIMENSION Z(*)
!
      CALL GF1(JD+1, Z(1) ,WTMOD)
      CALL GF1(JD+2, Z(2) ,WTMOD)
      CALL GF1(JD+3, Z(3) ,WTMOD)
      END 

      SUBROUTINE filiv(jd,z,wtmod)
!     ****************************
      IMPLICIT REAL*8(A-H,O-Z)
! Fills entries WITH symetrisation NW<=>WN
      DIMENSION z(*)
!
      CALL gf1(jd+1, 1d0-z(1) ,wtmod)
      CALL gf1(jd+2, 1d0-z(2) ,wtmod)
      CALL gf1(jd+3, 1d0-z(3) ,wtmod/2d0)
      CALL gf1(jd+3, 1d0-z(4) ,wtmod/2d0)
      END 


