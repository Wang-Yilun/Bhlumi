      SUBROUTINE ReaDat(Dname)
!     ************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*60 Dname
      SAVE
      COMMON / INOUT  / NINP,NOUT
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI
! -- THIS COMMON TRANSFERS TEST NAME FROM DATA SET TO PLOTING
      character*5 tdum

      WRITE(NOUT,*) '   '
      WRITE(NOUT,*) '   '
      WRITE(NOUT,*) '=============================================='
      WRITE(NOUT,*) '==========*********************==============='
      WRITE(NOUT,*) '==========***    PUBFIG     ***==============='
      WRITE(NOUT,*) '==========*********************==============='
      WRITE(NOUT,*) '=============================================='
      WRITE(NOUT,*) '   '
! ======================================
      write(6,*) '|ReaDat||',Dname,  '|||'
      NINPD =4
      OPEN( NINPd, file=Dname) 
! ======================================
      READ( NINPd,'(A5)'   ) tdum
      READ( NINPd,'(8I2)'  ) KAT1,KAT2,KAT3,KAT4,KAT5,KAT6
      READ( NINPd,'(I10)'  ) NEVT,KEYOPT,KEYRAD,KEYTRI
      READ( NINPd,'(F10.0)') CMSENE,TMING,TMAXG,VMAXG,XK0
      READ( NINPd,'(F10.0)') TMINW,TMAXW,TMINN,TMAXN,VMAXE
      READ( NINPd,'(I10)'  ) NPHI,NTHE
! ======================================
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
      WRITE(NOUT,'(2A12/2I12)')
     $  'NPHI','NTHE',
     $   NPHI,  NTHE

      END
