      PROGRAM MAIN
!     ***********************************
! This program is ploting energy dependence of the Z contribution
! The histo file includes "merged" results from several independent
! MC runs at several energies
!------
! To execute: make makfig1-dvi
!             make makfig1-ps
! note that ../prod1/prod1.hst may not exist
! then cd ../; make arch-makfig1-dvi
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
      OPEN( NOUT, file='pubfig.output')
      CALL GOUTPU(NOUT)
      Tesnam    = 'prod3'
      TeXfile   = 'prod3.tex'
      OPEN( NOUT, file='pubfig.output')
      CALL GOUTPU(NOUT)
!--------------------------------------------------------
! Stored histograms and corresponding histograms
!--------------------------------------------------------
      lendan = 0
! Results with changed parameters
      lendan = lendan+1
!-------
      TeXfile   = 'prod1.2M.Z-gamma.tex'
      Dname(lendan)  = '../prod1/prod1.data'
      Hname(lendan)  = '../prod1/prod1.hst'
!--------------------------------------------------------
!----    initialize GPLOT
!--------------------------------------------------------
      CALL GPLINT(0)
      NOUFIG=11
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!--------------------------------------------------------

!==================
      CALL prod1
!==================

! ------------dumping histogram for control -------------------
      NOUTH=20
      OPEN(NOUTH,file='dump.hst')
      CALL GRFILE(NOUTH,DNAME,'N')
      CALL GROUT( 0,ICY,' ')
      CALL GREND(DNAME)
! ------------THE END OF HISTO WRITING -------------------------
      CALL GPLEND
      CLOSE(NOUT)
      END


      SUBROUTINE prod1
!     *****************
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
      WRITE(NOUT,*) '=============  prod1     ====================='
      WRITE(NOUT,*) '=============================================='
!
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

! ==========================================================
! ==========================================================
! ==========================================================
      CALL pfi5
!
! ==========================================================
! ==========================================================
! ==========================================================
      END


      SUBROUTINE pfi5
!     *****************
! fig. 5 in CERN-TH/95-74
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
      DIMENSION ws0(10), ws1(10),wse(10),wsd(10)
      DIMENSION we0(10), we1(10),wee(10),wed(10)
!
      DIMENSION    idl(5)
      CHARACTER*16 capt(6)
      CHARACTER*8  fmt(3)


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

      AMZ  = 91.187D0
      dele = 0.375d0
      emin = AMZ-dele*5.5d0
      emax = AMZ+dele*4.5d0

      call gbook1(21,' SICAL: Born          $',10,emin,emax)
      call gbook1(22,' SICAL: Born+O(al)    $',10,emin,emax)
      call gbook1(23,' SICAL: Born+O(al)exp $',10,emin,emax)
      call gbook1(24,' SICAL: O(al)exp-O(al)$',10,emin,emax)

      DO i=1,10
      KeyGen = 3
      JDA    = 10000*KEYGEN
      KeyGei = KeyGen       +1000000*(i-1)
      JDi    = 10000*KeyGen +1000000*(i-1)
!
      DO itr=1,3
! O(alf2)expB real exp.
        CALL CUMHIS(KeyGei,JDi+100+itr,JDi+400+itr)
        CALL Gdelet(JDi+100+itr)
        CALL CUMHIS(KeyGei,JDi+200+itr,JDi+500+itr)
        CALL Gdelet(JDi+200+itr)
        CALL CUMHIS(KeyGei,JDi+300+itr,JDi+600+itr)
        CALL Gdelet(JDi+300+itr)
      ENDDO

      energy = emin + (i+0.5d0)*dele

      ws0(i) = gi(JDi+403,1)/BorN(2)
      ws1(i) = gi(JDi+503,1)/BorN(2)
      wse(i) = gi(JDi+603,1)/BorN(2)
      wsd(i) = wse(i)-ws1(i)

      we0(i) = gie(JDi+403,1)/BorN(2)
      we1(i) = gie(JDi+503,1)/BorN(2)
      wee(i) = gie(JDi+603,1)/BorN(2)
      wed(i) = wee(i)-we1(i)

      write(6,'(4f20.10)') energy,ws0(i),ws1(i),wse(i)

      IF(i.eq.2) THEN
! For Control: Z contribution single energy
      YMAX =  .003
      YMIN = -.003
c      CALL PL3MC('SICAL: WN, Z-Born Z-O(alf1), Z-O(alf1)exp, 1-Zmin $'
c     $,JDi+403,JDi+503,JDi+603,YMIN,YMAX, BorN(2),BorN(2),BorN(2))
c      CALL PL3MC('SICAL: WW, Z-Born Z-O(alf1), Z-O(alf1)exp, 1-Zmin $'
c     $,JDi+401,JDi+501,JDi+601,YMIN,YMAX, BorW(2),BorW(2),BorW(2))
      ENDIF
      ENDDO
!====================================================================
!                       plotting - SICAL
!====================================================================
      call gpak(21,ws0)
      call gpak(22,ws1)
      call gpak(23,wse)
      call gpak(24,wsd)

      call gpake(21,we0)
      call gpake(22,we1)
      call gpake(23,wee)
      call gpake(24,wed)

      YMAX =  .003
      YMIN = -.003
      CALL Gpltit('Fig. 5 in TH-95-74, Z-Born, O(alf1), O(alf1)exp$')
      CALL GMINIM(21,  ymin)
      CALL GMAXIM(21,  ymax)
      CALL GPLSET('DMOD',4D0)
      CALL GPLOT(21,   ' ','*',0)
      CALL GMINIM(22,  ymin)
      CALL GMAXIM(22,  ymax)
      CALL GPLSET('DMOD',3D0)
      CALL GPLOT(22,   'S','*',0)
      CALL GMINIM(23,  ymin)
      CALL GMAXIM(23,  ymax)
      CALL GPLSET('DMOD',6D0)
      CALL GPLOT(23,   'S','*',0)
! Inset
      CALL Gpltit('Fig. 5, inset, in TH-95-74, O(alf)exp-O(alf)$')
      YMAX =  .0002
      YMIN = -.0004
      CALL GMINIM(24,  ymin)
      CALL GMAXIM(24,  ymax)
      CALL GPLSET('DMOD',3D0)
      CALL GPLOT(24,   ' ','*',0)

!*********************************************************************
!                         the same as TABLE
!*********************************************************************
      CALL Gidopt(21,'SLAN')
      CALL Gidopt(21,'ERRO')
      CALL Gpltit('Fig. 5 in TH-95-74: Z-Born, O(alf1), O(alf1)exp$')
      capt(1)='CMS ene'
      capt(2)='Z-Born'
      capt(3)='O(alf1)'
      capt(4)='O(alf1)exp'
      capt(5)='alf1exp-alf1'
      fmt(1)='F12.3'
      fmt(2)='F10.6'
      fmt(3)='F10.6'
      idl(1)=21
      idl(2)=22
      idl(3)=23
      idl(4)=24
      CALL gpltab(4,idl,capt,fmt,1,1,0)
!=========================================================
      END

