

      SUBROUTINE BHLOG4(MODE,XPAR,NPAR)
C     *********************************  
C=================================================================== 
C===================================================================
C===BBBBBB=====BBB==BBB===BBB=======BBB=======BBBBBB=====BBBBBBB====
C===BBB==BB====BBB==BBB===BBB=======BBB======BBBBBBBB===BBBBBBBBB===
C===BBB==BB====BBB==BBB===BBB=======BBB======BBB==BBB===BBB=========
C===BBBBBB=====BBBBBBBB===BBB=======BBB======BBB==BBB===BBB=========
C===BBBBBBBB===BBBBBBBB===BBB=======BBB======BBB==BBB===BBB==BBBB===
C===BBB==BBB===BBB==BBB===BBBBBBB===BBBBBBB==BBB==BBB===BBB===BBB===
C===BBBBBBBB===BBB==BBB===BBBBBBB===BBBBBBB===BBBBBB=====BBBBBB=====
C===================================================================
C===================================================================
C
****************************************************
c This is simple LL generator with one-line emission
c Its purpose is to test analytical v-integration in O(alf4)prag.
c see BOKER4 and VVRHO routines in SUPERFIG
C     ***********************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      DIMENSION XPAR(*),NPAR(*)
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
      COMMON / INOUT  / NINP,NOUT
      COMMON / TRANSR / TRAN,TRMIN,TRMAX
      COMMON / BHPAR3 / KEYRAD,KEYOPT
      COMMON / BHPAR2 / CMSENE,AMEL  
      COMMON / BHPAR1 / DEL,EPSCM,THMIN,XMIVIS 
C Communicates with VARRAN
      COMMON / RANPAR / KEYRND
      COMMON / BHPCOL / PVI,PVF,QVI,QVF
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      SAVE / INOUT  /,/ TRANSR /,/ BHPAR3 /,/ BHPAR2 /,/ BHPAR1 /
      SAVE / RANPAR /,/ BHPCOL /,/ WGTALL /
C
      SAVE NEVGEN
C...
      IF(MODE.EQ.-1) THEN  
*     *******************
      KEYOPT = NPAR(1)
      KEYRAD = NPAR(2)
      CMSENE = XPAR(1)
      TRMIN  = XPAR(2)
      TRMAX  = XPAR(3)
      EPSCM  = XPAR(4)
      KEYRND = MOD(KEYOPT,10)   
C--
      AMEL  =  0.5111D-3
      TRAN  = (TRMIN+TRMAX)/2

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '  BHLOG4: INPUT PARAMETRES       '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXL1I) KEYOPT,     ' options   switch  ','KEYOPT','N1'
      WRITE(NOUT,BXL1I) KEYRND,     ' VARRAN    switch  ','KEYRND','  '
      WRITE(NOUT,BXL1I) KEYRAD,     ' radiation switch  ','KEYRAD','N2'
      WRITE(NOUT,BXL1F) CMSENE,     ' CMS enegry  [GeV] ','CMSENE','X1'
      WRITE(NOUT,BXL1G) TRMIN ,     ' trasf_min [GeV^2] ','TRMIN ','X2'
      WRITE(NOUT,BXL1G) TRMAX ,     ' trasf_min [GeV^2] ','TRMAX ','X3'
      WRITE(NOUT,BXL1G) TRAN  ,     ' trasfer   [GeV^2] ','TRMAS ','X3'
      WRITE(NOUT,BXL1F) EPSCM ,     ' eps_CM infr. cut  ','EPSCM ','X4'
      WRITE(NOUT,BXCLO)  

      NEVGEN=0

      ELSEIF(MODE.EQ.0) THEN
*     **********************
      NEVGEN=NEVGEN+1 

      KEYGEN = MOD(KEYOPT,10000)/1000

      CALL GENIEK(TRAN,VI,VF,VV,WTBASE,WTPDEL,WTCONV,WTCON2,WTCON3,WTD)
      WTSET(10)=WTBASE
      WTSET(12)=WTCONV
      WTSET(13)=WTCON2
      WTSET(14)=WTCON3
      WTSET(15)=WTPDEL
      WTSET(16)=WTD
      PVI=VI
      PVF=VF
      QVI=0
      QFV=0

      ELSEIF(MODE.EQ.1.OR.MODE.EQ.2) THEN 
*     ***********************************
      NPAR(10)= NEVGEN
      NPAR(20)= NEVGEN
c x-section equas one
      XPAR(10)= 1D0
      XPAR(11)= 0D0    
C for WEIGHTED events
      XPAR(20)= 1D0
      XPAR(21)= 0D0    

      ENDIF
      END

      SUBROUTINE GENIEK(TRAN,VI,VF,VV,
     $                 WTBASE,WTPDEL,WTCONV,WTCON2,WTCON3,WTDIF)
*     **********************************************************
c This is simple LL generator with two-line emission
c Its purpose is to test analytical v-integration in O(alf4)prag.
c see BOKER4 and VVRHO routines in SUPERFIG
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI     = 3.1415926535897932D0, ALFA=1D0/137.03604D0) 
      PARAMETER(ALF1   = ALFA/PI)  
      PARAMETER(AMEL   = 0.5111D-3)
      PARAMETER(CEULER = 0.57721566D0)       
      COMMON / INOUT  / NINP,NOUT     
      dimension drvec(3)
      save
      data icont /0/
      icont = icont+1
C=====================================================
      BILG   = DLOG(TRAN/AMEL**2)
      gami   =   ALF1*(BILG-1D0) 
      GAMFAC = EXP(-CEULER*gami)/DPGAMM(1D0+gami)
      call varran(drvec,2)
      r1 = drvec(1)
      r2 = drvec(2)
c v-variables
c by hand regularization to avoid r1=0.d0
      r1 = max(r1, 1d-10)
      r2 = max(r2, 1d-10)
      vi = r1**(1d0/gami)    
      vf = r2**(1d0/gami)    
      vv = vi+vf
c various weights
      WTBASE = gamfac**2
c
      WTPDEL = 0
      WTCON2 = 0
      DISBAS= 0
      WTCONV= 0
      IF(VV.LT. 0.9999999) then
        WTPDEL = WTBASE * EXP(-2*gami*LOG(1-VI))
        DISBAS=  gamfac*gami*vi**(gami-1)
     $          *gamfac*gami*vf**(gami-1)
        TRANP = TRAN*(1D0-VI)**2/(1D0-VV)
        GAMIP = ALF1*(DLOG(TRANP/AMEL**2) -1D0) 
        GAMFP = EXP(-CEULER*gamip)/DPGAMM(1D0+gamip)
        DISTP =  gamfp**2*gamip**2
     $           *vi**(gamip-1)*vf**(gamip-1)
        WTCONV= WTBASE *DISTP/DISBAS 

        WTCON2= WTBASE *DISTP/DISBAS *EXP(-2*gamip*LOG(1-VI))

*! WTCON3 is intermediate step WTBASE--> WTCONV
*! effective powers in vi**(gamip-1)*vf**(gamip-1)
*! excelent approximation ( below 0.2E-4)
        d  =  ALF1*LOG(1D0-VV)
        DISTZ =  gamfac**2 *gamip**2
     $          *vi**(gami-1-d) *vf**(gami-1+d)
        WTCON3= WTBASE *DISTZ/DISBAS

*------------------------CONV--------------------------
*! the only approximation is as in CON3 effective powers
c  otherwise pure rearrangements
c        d     =    alf1*log(1-vv)
c        DISTZ =  gamfac**2 *(
c pole part of gamip**2 seems to be below 0.3E-4
c anyway, it is integrable by hand if necessary!
c     $   +(vi/vv*(gami+d)**2+ vf/vv*(gami-d)**2 )
c exact residue part of gamip**2
c     $   +4*alf1*(gami-d)*(log(1-vi)-(vi/vv)*log(1-vv))
c     $   +4*alf1**2*(log(1-vi)**2-(vi/vv)*log(1-vv)**2)
c     $        )*vi**(gami-1-d) *vf**(gami-1+d)
c        WTCON4= WTBASE *DISTZ/DISBAS
c! the version directly integrable by hand
        d     =    alf1*log(1-vv)
        DISTZ =  gamfac**2 *(
c pole part of gamip**2 seems to be below 0.3E-4
c anyway, it is integrable by hand if necessary!
     $ +(vi/vv*(gami+d)**2+ vf/vv*(gami-d)**2 )
     $            *vi**(gami-1-d) *vf**(gami-1+d)
c approximate residue part integrable analytically
     $ +4*alf1*(gami-d)   *vv**(2*gami)
     $    *(log(1-vi)-(vi/vv)*log(1-vv))/(vi*vf)
c new terms beyond 2-nd pragmatic
c     $ +4*alf1*gami**2    *vv**(2*gami)
c     $    *(log(1-vi)-(vi/vv)*log(1-vv))/(vi*vf) *log(vi*vf/vv/vv)   
c     $ +4*alf1**2 
c     $    *(log(1-vi)**2-(vi/vv)*log(1-vv)**2)/(vi*vf)
     $        )
        WTCON4= WTBASE *DISTZ/DISBAS
        WTDIF = WTCON4-WTBASE
cc        WTDIF = WTCONV-WTCON4

*------------------------CON2--------------------------
*! truncated WTCON2 such that disruption is below 0.0001
        d  =  ALF1*LOG(1D0-VV)
cc        DISTU =  gamfac**2 *gamip**2
cc     $   *vi**(gami-1-d)*vf**(gami-1+d)
cc     $   *(1 -2*(gami+d )*LOG(1-VI) +2*gami**2*LOG(1-VI)**2  )
*! now pole approximations, disruption is below 0.0001 !!!
cc        DISTU =  
cc     $  +1
cc     $      *gamfac**2 *vi**(gami-1-d)*vf**(gami-1+d)*gamip**2
cc     $  -2*(gami+d )*( vi/vv*LOG(1-vv)*(gami+d)**2 )
cc     $      *gamfac**2 *vi**(gami-1-d)*vf**(gami-1+d)
cc     $  +2*vi/vv*(gami+d)**2*(gami+d)**2 *LOG(1-vv)**2
cc     $      *gamfac**2 *vi**(gami-1-d)*vf**(gami-1+d)
*! just reordered result
cc      gamip =  gami +alf1*log((1-vi)**2/(1-vv))
        DISTU =  gamfac**2 *(
     $    vi**(gami-1-d)*vf**(gami-1+d)*gamip**2
* new terms with respect to WTCONV below
     $  -2*(gami+d)**3 *LOG(1-vv)/vv
     $      *vi**(gami-d)*vf**(gami-1+d)
     $  +2*(gami+d)**4 *LOG(1-vv)**2/vv
     $      *vi**(gami-d)*vf**(gami-1+d)
     $  )
*! approximations below cost -0.0001 (only for analyt. v-integration)
cc        DISTU =  gamfac**2 *(
cc     $    vi**(gami-1-d)*vf**(gami-1+d)*gamip**2
* new terms with respect to WTCONV below
cc     $  -2*(gami+d)**3 *LOG(1-vv)/vv
cc     $      *vf**(gami-1)
cc     $  +2*(gami)**4 *LOG(1-vv)**2/vv
cc     $      *vf**(gami-1)
cc     $  )
        WTCON5= WTBASE *DISTU/DISBAS
cc        WTDIF = WTCON5 -WTCON2
      ENDIF
      END     



      SUBROUTINE BHLOG5(MODE,XPAR,NPAR)
C     *********************************  
c This is simple LL generator for 2-line emission.
c It implements NEW BASELINE distribution and beta0
C     ***********************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      DIMENSION XPAR(*),NPAR(*)
      CHARACTER*80   BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G 
      PARAMETER (PI     = 3.1415926535897932D0, ALFA=1D0/137.03604D0) 
      PARAMETER (ALF1   = ALFA/PI)  
      PARAMETER (CEULER = 0.57721566D0)       
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
      COMMON / TRANSR / TRAN,TRMIN,TRMAX
      COMMON / BHPAR3 / KEYRAD,KEYOPT
      COMMON / BHPAR2 / CMSENE,AMEL  
      COMMON / BHPAR1 / DEL,EPSCM,THMIN,XMIVIS 
C Communicates with VARRAN
      COMMON / RANPAR / KEYRND
      COMMON / BHPCOL / PVI,PVF,QVI,QVF
      COMMON / WGTALL / WTMOD,WTCRU1,WTCRU2,WTSET(300)
      SAVE / INOUT  /,/ TRANSR /,/ BHPAR3 /,/ BHPAR2 /,/ BHPAR1 /
      SAVE / RANPAR /,/ BHPCOL /,/ WGTALL /
C
      SAVE NEVGEN
C...
      IF(MODE.EQ.-1) THEN  
*     *******************
      KEYOPT = NPAR(1)
      KEYRAD = NPAR(2)
      CMSENE = XPAR(1)
      TRMIN  = XPAR(2)
      TRMAX  = XPAR(3)
      EPSCM  = XPAR(4)
      KEYRND = MOD(KEYOPT,10)   
C--
      AMEL  =  0.5111D-3
      TRAN  = (TRMIN+TRMAX)/2
      DEL   = EPSCM

      WRITE(NOUT,BXOPE)
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXTXT) '  BHLOG5: INPUT PARAMETRES       '
      WRITE(NOUT,BXTXT) '*********************************'
      WRITE(NOUT,BXL1I) KEYOPT,     ' options   switch  ','KEYOPT','N1'
      WRITE(NOUT,BXL1I) KEYRND,     ' VARRAN    switch  ','KEYRND','  '
      WRITE(NOUT,BXL1I) KEYRAD,     ' radiation switch  ','KEYRAD','N2'
      WRITE(NOUT,BXL1F) CMSENE,     ' CMS enegry  [GeV] ','CMSENE','X1'
      WRITE(NOUT,BXL1G) TRMIN ,     ' trasf_min [GeV^2] ','TRMIN ','X2'
      WRITE(NOUT,BXL1G) TRMAX ,     ' trasf_min [GeV^2] ','TRMAX ','X3'
      WRITE(NOUT,BXL1G) TRAN  ,     ' trasfer   [GeV^2] ','TRMAS ','X3'
      WRITE(NOUT,BXL1F) EPSCM ,     ' eps_CM infr. cut  ','EPSCM ','X4'
      WRITE(NOUT,BXCLO)  

      NEVGEN=0

      ELSEIF(MODE.EQ.0) THEN
*     **********************
      NEVGEN=NEVGEN+1 

      CALL GENIAK(TRAN,AMEL,DEL,PVi,PVf, wp0,wp1,wp2,wp3)
      CALL GENIAK(TRAN,AMEL,DEL,QVi,QVf, wq0,wq1,wq2,wq3)

c one line case
      WTSET(10) = wp0
      WTSET(11) = Wp1
      WTSET(12) = Wp2
      WTSET(13) = wp3
c two-line case
      WTSET(20) = wp0*wq0
      WTSET(21) = wp1*wq1
      WTSET(22) = wp2*wq2
      WTSET(23) = wp3*wq3

      ELSEIF(MODE.EQ.1.OR.MODE.EQ.2) THEN 
*     ***********************************
      NPAR(10)= NEVGEN
      NPAR(20)= NEVGEN
c x-section equas one
      XPAR(10)= 1D0
      XPAR(11)= 0D0    
C for WEIGHTED events
      XPAR(20)= 1D0
      XPAR(21)= 0D0    
      ENDIF
c     *****
      END

      SUBROUTINE GENIAK(TRAN,AMEL,DEL,Vi,Vf,wt0,wt1,wt2,wt3)
*     ******************************************************
c This is simple LL generator for 2-line emission.
c It implements NEW BASELINE distribution and beta0
c ----------------  distributions ------------------
c wt0 = EXP(DELB1) trivial constant weight
c wt1 = covolution of two NewBaseline (exact)
c wt2 = covolution of two NewBaseline (approx)
c wt3 = convolution with exp(-gam*log(1-vi)+DELB) YFS formfactor
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI     = 3.1415926535897932D0, ALFA=1D0/137.03604D0) 
      PARAMETER (ALF1   = ALFA/PI)  
      PARAMETER (NMAX=25)
      COMMON / INOUT  / NINP,NOUT     
      SAVE   / INOUT  /
      DIMENSION RR(100),VVI(100),VVF(100)

      BILG   = DLOG(TRAN/AMEL**2)
      BETI   = 2D0*ALF1*(BILG-1D0) 
      DELB   = 0.25D0*BETI -ALF1*(0.5D0 +1/6D0*PI**2)
      DELB1  = -ALF1*1/6D0*PI**2
      FYFS   = EXP(DELB)
      FYFS0  = EXP(0.5*DELB1)
      FYFS1  = EXP(DELB1)
      GAMH   = BETI/2
      AVERG  = GAMH*LOG(1D0/DEL)
c initial line
 10   CONTINUE
      CALL POISSG(AVERG,NMAX,NPH1,RR) 
      IF(NPH1.GT.NMAX) GOTO 10  
      Vi   =0
      proi =1
      sumi =0
      DO 80 I=1,NPH1  
        V      = DEL**RR(I)
        VVI(i) = V
        Vi     = Vi +V
        proi   = proi*(GAMH - 0.5D0*ALF1*LOG(1D0-V))/gamh
        sumi   = sumi - 0.5D0*ALF1*LOG(1D0-V)
 80   CONTINUE
      wti1 = FYFS0*proi
      wti2 = FYFS0*(1 +sumi/gamh) 
c final line
 110  CONTINUE
      CALL POISSG(AVERG,NMAX,NPH2,RR) 
      IF(NPH2.GT.NMAX) GOTO 110  
      Vf   =0
      prof =1
      sumf =0
      DO 180 I=1,NPH2  
        V      = DEL**RR(I)
        VVF(i) = V
        Vf     = Vf +V
        prof   = prof*(GAMH - 0.5D0*ALF1*LOG(1D0-V))/gamh
        sumf   = sumf - 0.5D0*ALF1*LOG(1D0-V)
 180  CONTINUE
      wtf1 = FYFS0*prof
      wtf2 = FYFS0*(1 +sumf/gamh) 
c-------------------------
c merge two lines
      VSUM  = Vi + Vf
      wt0   = EXP(DELB1)
      wt1   = wti1*wtf1
      wt2   = wti2*wtf2
      wt3=0
      if(vsum.lt.1d0) then
        GAMHb = GAMH +ALF1*LOG((1-vi)**2/(1-vsum))
        dl    =       ALF1*LOG((1-vi)**2/(1-vsum))
        GAMB  = 2*GAMHb
c recalculate weight with gamma-bar
        proi =1
        sumi =0
        DO 281 I=1,NPH1
        sumi = sumi    +dl -ALF1/2*LOG(1D0-VVI(i))
        proi = proi*(GAMHb -ALF1/2*LOG(1D0-VVI(i)))/gamh
 281    CONTINUE
        prof =1
        sumf =0
        DO 282 I=1,NPH2
        sumf = sumf    +dl -ALF1/2*LOG(1D0-VVF(i))
        prof = prof*(GAMHb -ALF1/2*LOG(1D0-VVF(i)))/gamh
 282    CONTINUE
        wt3 = proi*prof
     $    *exp(-2*alf1*log((1-vi)**2/(1-vsum))*log(1/del))
     $    *exp(DELB -gamb*log(1-vi))
      endif
      END

      SUBROUTINE GENIAL(TRAN,AMEL,DEL,VSUM,wt1,wt2)
*     *********************************************
c This is simple LL generator for one-line emission.
c It implements NEW BASELINE distribution for one line
c WT0 is the old trivial baseline F(gam)*gam*v**(gam-1)
c WT1 is NEW beseline (NBAS) exact
c WT2 is NEW beseline (NBAS) truncated to O(alf3)prag
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI     = 3.1415926535897932D0, ALFA=1D0/137.03604D0) 
      PARAMETER (ALF1   = ALFA/PI)  
      PARAMETER (NMAX=25)
      COMMON / INOUT  / NINP,NOUT     
      SAVE   / INOUT  /
      DIMENSION RR(100)

      GAMH   = ALF1*(LOG(TRAN/AMEL**2)-1)
      AVERG  = GAMH*LOG(1D0/DEL)
 10   CONTINUE
      CALL POISSG(AVERG,NMAX,NPH,RR) 
      IF(NPH.GT.NMAX) GOTO 10  
      VSUM=0D0
      prod1=1
      sum1 =0
      DO 80 I=1,NPH  
        V  = DEL**RR(I)
        VSUM = VSUM +V
        prod1= prod1*(GAMH - 0.5D0*ALF1*LOG(1D0-V))
        sum1 = sum1 - 0.5D0*ALF1*LOG(1D0-V)
 80   CONTINUE
      wt1 =  prod1/gamh**nph
      wt2 =  (1 + sum1/gamh) 
      END

      SUBROUTINE GENIUL(TRAN,DEL,VSUM, wt0,wt1,wt2)
*     *********************************************
c This is simple LL generator with one-line emission.
c It implements NEW BASELINE distribution for one line
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI     = 3.1415926535897932D0, ALFA=1D0/137.03604D0) 
      PARAMETER (ALF1   = ALFA/PI)  
      PARAMETER (CEULER = 0.57721566D0)       
      PARAMETER (AMEL   = 0.5111D-3)
      PARAMETER (NMAX=25)
      COMMON / INOUT  / NINP,NOUT     
      SAVE   / INOUT  /
      DIMENSION RR(100)

ccc      GAMH   = ALF1*LOG(TRAN/AMEL**2)
      GAMH   = 2*ALF1*LOG(TRAN/AMEL**2)
      AVERG  = GAMH*LOG(1D0/DEL)
 10   CONTINUE
      CALL POISSG(AVERG,NMAX,NPH,RR) 
      IF(NPH.GT.NMAX) GOTO 10  
      VSUM=0D0
      prod1=1
      sum1 =0
      DO 80 I=1,NPH  
        V  = DEL**RR(I)
        VSUM = VSUM +V
        prod1= prod1*(GAMH - ALF1*LOG(1D0-V))
        sum1 = sum1 - ALF1*LOG(1D0-V)
 80   CONTINUE
      wt1 =  prod1/gamh**nph
      wt2 =  (1 + sum1/gamh) 

      wt0 =       EXP(-1/6D0*ALF1*PI**2)
      wt1 =  wt1 *EXP(-1/6D0*ALF1*PI**2)
      wt2 =  wt2 *EXP(-1/6D0*ALF1*PI**2)
      END
