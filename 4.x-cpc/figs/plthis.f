      SUBROUTINE PL1AN(TITLE,ID1,YMIN,YMAX)
*     **************************************************
*! Single plot, analytical or MC. 
*!  TITLE      title of the plot
*!  ID1,       plots existing
*!  YMIN,YMAX  uper lower limit in vertical scale of the plot
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE

      CALL GPLTIT(TITLE)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
      CALL GMINIM(ID1,    YMIN)
      CALL GMAXIM(ID1,    YMAX)
      CALL GPRINT(ID1)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(ID1,   ' ',' ',0)
      END

      SUBROUTINE PL2AN(TITLE,ID1,ID2,ID3,YMIN,YMAX,XMAG)
*     **************************************************
*! Two plots and difference 
*!  TITLE      title of the plot
*!  ID1,ID2    two plots existing
*!  ID3        difference (nonexisting)
*!  YMIN,YMAX  uper lower limit in vertical scale of the plot
*!  XMAG       magnification of the difference
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE

      CALL GPLTIT(TITLE)
*! divide by Born
      CALL GOPERA(ID1,'-',ID2,ID3, XMAG, XMAG)  
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
      CALL GMINIM(ID1,    YMIN)
      CALL GMAXIM(ID1,    YMAX)
      CALL GMINIM(ID2,    YMIN)
      CALL GMAXIM(ID2,    YMAX)
      CALL GMINIM(ID3,    YMIN)
      CALL GMAXIM(ID3,    YMAX)
      CALL GPRINT(ID1)
      CALL GPRINT(ID2)
      CALL GPRINT(ID3)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(ID1,   ' ',' ',0)
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(ID2,   'S','*',0)
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(ID3,   'S','*',0)
c purge local histos
      END

      SUBROUTINE HISUNI(ID1,ID2)
C     **************************
C creates unit'ary histogram id2 with binning of id1
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      DIMENSION YY(200),ER(200) 
      CHARACTER*80 TITLE

      CALL GINBO1(ID1,TITLE,NBIV,VMIN,VMAX)
      CALL GBOOK1(ID2,' unit histo',NBIV,VMIN,VMAX)
      DO 30 K=1,NBIV
        YY(K)= 1.D0
        ER(K)= 0.D0
 30   CONTINUE
      CALL GPAK (ID2,YY)
      CALL GPAKE(ID2,ER)
      END


      FUNCTION DSIG0(TRMIN,TRMAX)
C     *******************************
C BORN XSECTION pure t-channel low angle approximation    
C result in nanobarns
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      PARAMETER( PI = 3.1415926535897932D0, ALFINV = 137.03604D0)
      PARAMETER( ALFPI=  1D0/PI/ALFINV ,ALFA=1D0/ALFINV)
      PARAMETER( GNANOB=389.385D-30*1.D33 )

      DSIG0 = 4D0*ALFA**2*PI*(1D0/TRMIN-1D0/TRMAX)*GNANOB  
      END

      FUNCTION BORNB(CMSENE,THA,THB)
C     *******************************
C BORN XSECTION pure t-channel     
C THA,THB are in radians, CMSENE in GEV units
C result in nanobarns
C     ***********************************
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
C     ******************
C Integrand of BORNB
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BORXI=1D0/XI**2*(1D0+(1D0-XI)**2)/2D0
      END

      SUBROUTINE RENHST(KeyGen,CHAK,ID1,ID2)     
C     **************************************
*! This routine normalizes to x-section in NANOBARNS or to unity
*  CHAK = 'NB  '    normal case [nb]
*  CHAK = 'NB10'    log10 x-scale assumed [nb]
*  CHAK = 'UNIT'    normalization to unity
*! id2.ne.id1 required
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
* CMONIT COMMUNICATES WITH GMONIT         
      COMMON / CMONIT/ AVERWT,ERRELA,NEVTOT,NEVACC,NEVNEG,NEVOVE,NEVZER 
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
C
      IF(ID2.EQ.ID1) GOTO 900
      CALL GMONIT(   1,KeyGen,AWTOT,DWTOT,WWMX )
      NEVT    =  NEVTot
      XSCRNB  =  AVERWT
      CALL GINBO1(ID1,TITLE,NBT,TMIN,TMAX)
      FACT=1D0
      IF( CHAK.EQ. 'NB  ') THEN
        FACT = NBT*XSCRNB/(NEVT*(TMAX-TMIN))
      ELSEIF( CHAK.EQ. 'NB10') THEN
        FLN10 = LOG(10.)
        FACT = NBT*XSCRNB/(NEVT*(TMAX-TMIN)*FLN10)
      ELSE
        WRITE(6,*) '+++++ RENHST: wrong chak=',chak
      ENDIF 
C Multiply content
      CALL GOPERA(ID1,'+',ID1,ID2, FACT, 0D0)
      RETURN
 900  WRITE(6,*) '+++++ RENHST: ID1=ID2=',ID1
      END  

      SUBROUTINE CUMHIS(KeyGen,ID1,ID2)
C     *********************************
*! cumulates histogram content starting from underflow
*! and normalizes to the total x-section in NANOBARNS
*! id2.ne.id1 required
*     *********************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
* CMONIT COMMUNICATES WITH GMONIT         
      COMMON / CMONIT/ AVERWT,ERRELA,NEVTOT,NEVACC,NEVNEG,NEVOVE,NEVZER 
      CHARACTER*80 TITLE
      DIMENSION  X(400),ER(100)
      LOGICAL Gexist 

      IF(ID2.EQ.ID1) GOTO 900
      IF(Gexist(ID2)) THEN
        WRITE(6,*) 'CUMHIS: warning histo deleted ID2 = ',ID2
        CALL GDELET(ID2)
      ENDIF
      CALL GMONIT(   1,KeyGen,AWTOT,DWTOT,WWMX )
      NEVT    =  NEVTot
      XSCRNB  =  AVERWT
      IF(NEVT.EQ.0) GOTO 901
      CALL GINBO1(ID1,TITLE,NBT,TMIN,TMAX) 
      SWT  = GI (ID1,0)
      SSWT = GIE(ID1,0)**2
      DO 100 I=1,NBT
      SWT   = SWT + GI (ID1,I)  
      SSWT  = SSWT+ GIE(ID1,I)**2  
C note NEVT in error calc. is for the entire sample related 
C to the crude x-section XCRU including !!! zero weight events !!!!
      XSEC  = 0D0
      ERREL = 0D0
      IF(SWT.NE.0D0 .AND. NEVT.NE.0) THEN
        XSEC  = SWT*(XSCRNB/NEVT)  
        ERREL = SQRT(ABS(SSWT/SWT**2-1D0/FLOAT(NEVT)))
      ENDIF
      X(I)  = XSEC       
      ER(I) = XSEC*ERREL  
 100  CONTINUE
*! store result in id2
      CALL GBOOK1(ID2,TITLE,NBT,TMIN,TMAX)
      CALL GPAK(  ID2,X)
      CALL GPAKE( ID2,ER)   
      RETURN
 900  WRITE(6,*) '+++++ CUMHIS: ID1=ID2=',ID1
      RETURN
 901  WRITE(6,*) '+++++ CUMHIS: EMPTY HISTO ID=',ID1
      END      



c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE PLDMCB
     $   (CHAK,TITLE,ID1,ID2,IDREF,YMIN,YMAX,XMAG,YMAG,BORN1)
*     ***************************************************************
* For seminanalytical part, BHLUMI is assumed
!---------------------------------------------------------------------
!  TITLE      title of the plot
!  ID1        histogram MC (not cumulative)
!  ID2        histogram MC (not cumulative)
!  IDREF      histogram analytical cumulative - background reference! 
!  YMIN,YMAX  uper lower limit in vertical scale of the plot
!  XMAG       magnification for difference ID2-ID1.
!  YMAG       magnification for ID1 and ID2.
!---------------------------------------------------------------------
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      LOGICAL gexist

      BORN=ABS(BORN1)
      CALL GPLTIT(TITLE)
*! turn into cumulative normalized distribution or renormalize
      IDCU1  = ID1+500
      IDCU2  = ID2+500
      if(gexist(IDCU1)) write(6,*) ' ++++PLDMCB: Warning! IDCU1= ',IDCU1
      if(gexist(IDCU2)) write(6,*) ' ++++PLDMCB: Warning! IDCU2= ',IDCU2
      IF(CHAK.eq.'CUMU') THEN
         CALL CUMHIS(3,ID1,IDCU1)
         CALL CUMHIS(3,ID2,IDCU2)
      ELSE
         CALL RENHST(3,CHAK,ID1,IDCU1)
         CALL RENHST(3,CHAK,ID2,IDCU2)
      ENDIF
*! divide by Born
      JUNIT=99999
      IF(BORN1.EQ.0D0) goto 901
      CALL HISUNI(ID1,JUNIT)
      CALL GOPERA(IDCU1,'/',JUNIT,IDCU1, 1D0, BORN) 
      CALL GOPERA(IDCU2,'/',JUNIT,IDCU2, 1D0, BORN) 
*! optionaly subtract one 
      IF(BORN1.LT.0D0) THEN
        CALL GOPERA(IDCU1,'-',JUNIT,IDCU1, 1D0, 1D0) 
        CALL GOPERA(IDCU2,'-',JUNIT,IDCU2, 1D0, 1D0)
      ENDIF 
*! calculate difference with respect to NREF (MC-Ref)
      IF(IDREF.NE.0) THEN
        CALL GOPERA(IDCU1,'-',IDREF,IDCU1, 1D0, 1D0) 
        CALL GOPERA(IDCU2,'-',IDREF,IDCU2, 1D0, 1D0) 
      ENDIF 
*! the difference  ID1-ID2 times XMAG
      IDIFF =  ID1+200
      if(gexist(IDIFF)) write(6,*) ' ++++PLDMCB: Warning! IDIFF= ',IDIFF
      CALL GOPERA(IDCU1,'-',IDCU2,IDIFF, XMAG, XMAG)
*! magnify
      CALL GOPERA(IDCU1,'*',JUNIT,IDCU1, 1D0, YMAG)  
      CALL GOPERA(IDCU2,'*',JUNIT,IDCU2, 1D0, YMAG)  
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GMINIM(IDCU1,    YMIN)
      CALL GMAXIM(IDCU1,    YMAX)
      CALL GIDOPT(IDCU1,   'ERRO')
      CALL GPRINT(IDCU1)
      CALL GMINIM(IDCU2,    YMIN)
      CALL GMAXIM(IDCU2,    YMAX)
      CALL GPRINT(IDCU2)
      CALL GMINIM(IDIFF,    YMIN)
      CALL GMAXIM(IDIFF,    YMAX)
      CALL GIDOPT(IDIFF,   'ERRO')
      CALL GPRINT(IDIFF)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCU1,   ' ',' ',0)
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCU2,   'S','*',0)
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDIFF   ,'S','*',0)
c purge local histos
      CALL GDELET(JUNIT)
      CALL GDELET(IDCU1)
      CALL GDELET(IDCU2)
      CALL GDELET(IDIFF)
      RETURN
 901  WRITE(6,*) ' PLDMCB: wrong born!!!!'
      END


      SUBROUTINE PLDIFB
     $   (CHAK,TITLE,ID,IDAN,IDREF,YMIN,YMAX,XMAG,YMAG,BORN1)
*     ***************************************************************
* For seminanalytical part, BHLUMI is assumed
*! Absolute 
*!  TITLE      title of the plot
*!  ID         histogram MC (not cumulative)
*!  IDAN       histogram analytical cumulative!
*!  IDREF      histogram analytical cumulative - background reference! 
*!  YMIN,YMAX  uper lower limit in vertical scale of the plot
*!  XMAG       magnification for difference MC-Analit.
*!  YMAG       magnification for ID and IDAN.
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      LOGICAL gexist

      BORN=ABS(BORN1)
      CALL GPLTIT(TITLE)
*! turn into cumulative normalized distribution or renormalize
      IDCUM  = ID+400
      if(gexist(IDCUM)) write(6,*) ' +++PLIFB: Warning! IDCUM = ',IDCUM
      IF(CHAK.eq.'CUMU') THEN
         CALL CUMHIS(3,ID,IDCUM)
      ELSE
         CALL RENHST(3,CHAK,ID,IDCUM)
      ENDIF
*! divide by Born
      JUNIT=99999
      IF(BORN1.EQ.0D0) goto 901
      CALL HISUNI(ID,JUNIT)
      CALL GOPERA(IDCUM,'/',JUNIT,IDCUM, 1D0, BORN) 
*! optionaly subtract one 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM,'-',JUNIT,IDCUM, 1D0, 1D0) 
*! calculate difference with respect to NREF (MC-Ref)
      IF(IDREF.NE.0)
     $  CALL GOPERA(IDCUM,'-',IDREF,IDCUM, 1D0, 1D0) 
*! the difference MC-Analit. times XMAG
      IDIFF =  ID+500
      if(gexist(IDIFF)) write(6,*) ' +++PLIFB: Warning! IDIFF = ',IDIFF
      CALL GOPERA(IDCUM,'-',IDAN,IDIFF, XMAG, XMAG)
*! magnify
      IDAN2 = IDAN+1000
      if(gexist(IDAN2)) write(6,*) ' +++PLIFB: Warning! IDAN2 = ',IDIFF
      CALL GOPERA(IDCUM,'*',JUNIT,IDCUM, 1D0, YMAG)  
      CALL GOPERA(IDAN ,'*',JUNIT,IDAN2, 1D0, YMAG)  
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GMINIM(IDCUM,    YMIN)
      CALL GMAXIM(IDCUM,    YMAX)
      CALL GIDOPT(IDCUM,   'ERRO')
      CALL GPRINT(IDCUM)
      CALL GMINIM(IDAN2,    YMIN)
      CALL GMAXIM(IDAN2,    YMAX)
      CALL GPRINT(IDAN2)
      CALL GMINIM(IDIFF,    YMIN)
      CALL GMAXIM(IDIFF,    YMAX)
      CALL GIDOPT(IDIFF,   'ERRO')
      CALL GPRINT(IDIFF)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCUM,   ' ',' ',0)
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDAN2,   'S','*',0)
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDIFF   ,'S','*',0)
c purge local histos
      CALL GDELET(JUNIT)
      CALL GDELET(IDCUM)
      CALL GDELET(IDIFF)
      CALL GDELET(IDAN2)
      RETURN
 901  WRITE(6,*) ' PLIFB: wrong born!!!!'
      END



      SUBROUTINE PLONEB(TITLE,ID,YMIN,YMAX,BORN1)
*     ******************************************
* For seminanalytical part, BHLUMI is assumed
*! One Monte Carlo histo 
*!  TITLE      title of the plot
*!  ID         histogram MC (not cumulative)
*!  YMIN,YMAX  uper lower limit in vertical scale of the plot
*!  Born       content divided by Born
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE

      BORN=ABS(BORN1)
      CALL GPLTIT(TITLE)
*! turn into cumulative normalized distribution
      IDCUM = ID+200
      CALL CUMHIS(3,ID,IDCUM)
*! divide by Born
      JUNIT=99999
      CALL HISUNI(ID,JUNIT)
      CALL GOPERA(IDCUM,'/',JUNIT,IDCUM, 1D0, BORN)  
*! optionaly subtract one 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM,'-',JUNIT,IDCUM, 1D0, 1D0) 
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
      CALL GMINIM(IDCUM,    YMIN)
      CALL GMAXIM(IDCUM,    YMAX)
      CALL GIDOPT(IDCUM,   'ERRO')
      CALL GPRINT(IDCUM)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDCUM,   ' ','*',0)
c purge local histos
      CALL GDELET(IDCUM)
      CALL GDELET(JUNIT)
      END


      SUBROUTINE PLTWOB(TITLE,ID1,ID2,YMIN,YMAX,BORN1)
*     ***********************************************
* For seminanalytical part, BHLUMI is assumed
*! Two Monte Carlo histos
*!  TITLE      title of the plot
*!  ID1,ID2    histograms MC (not cumulative)
*!  YMIN,YMAX  uper lower limit in vertical scale of the plot
*!  Born       content divided by Born
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE

      BORN=ABS(BORN1)
      CALL GPLTIT(TITLE)
*! turn into cumulative normalized distribution
      IDCUM1 = ID1+200
      IDCUM2 = ID2+200
      CALL CUMHIS(3,ID1,IDCUM1)
      CALL CUMHIS(3,ID2,IDCUM2)
*! divide by Born
      JUNIT=99999
      CALL HISUNI(ID1,JUNIT)
      CALL GOPERA(IDCUM1,'/',JUNIT,IDCUM1, 1D0, BORN)  
      CALL GOPERA(IDCUM2,'/',JUNIT,IDCUM2, 1D0, BORN)  
*! optionaly subtract one 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM1,'-',JUNIT,IDCUM1, 1D0, 1D0) 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM2,'-',JUNIT,IDCUM2, 1D0, 1D0) 
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
      CALL GMINIM(IDCUM1,    YMIN)
      CALL GMAXIM(IDCUM1,    YMAX)
      CALL GIDOPT(IDCUM1,   'ERRO')
      CALL GPRINT(IDCUM1)
      CALL GMINIM(IDCUM2,    YMIN)
      CALL GMAXIM(IDCUM2,    YMAX)
      CALL GIDOPT(IDCUM2,   'ERRO')
      CALL GPRINT(IDCUM2)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDCUM1,   ' ','*',0)
      CALL GPLSET('DMOD',4D0)
      CALL GPLOT(IDCUM2,   'S','*',0)
c purge local histos
      CALL GDELET(IDCUM1)
      CALL GDELET(IDCUM2)
      CALL GDELET(JUNIT)
      END

*=======================================================================
*=======================================================================
*=======================================================================
*=======================================================================
*==================== trigger specific part ============================
*=======================================================================
*=======================================================================
*=======================================================================
*=======================================================================


      SUBROUTINE PLDIFF
     $   (TITLE,ID,IDAN,IDREF,YMIN,YMAX,XMAG,YMAG,BORN1)
*     ***************************************************************
*! Absolute 
*!  TITLE      title of the plot
*!  ID         histogram MC already renormalized!!!
*!  IDAN       histogram analytical cumulative!
*!  IDREF      histogram analytical cumulative - background reference! 
*!  YMIN,YMAX  uper lower limit in vertical scale of the plot
*!  XMAG       magnification for difference MC-Analit.
*!  YMAG       magnification for ID and IDAN.
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE
      LOGICAL gexist

      BORN=ABS(BORN1)
      CALL GPLTIT(TITLE)
      IDCUM  = ID+400
      if(gexist(IDCUM)) write(6,*) ' +++PLDIFF: Warning! IDCUM = ',IDCUM
*! divide by Born
      JUNIT=99999
      IF(BORN1.EQ.0D0) goto 901
      CALL HISUNI(ID,JUNIT)
      CALL GOPERA(ID,'/',JUNIT,IDCUM, 1D0, BORN) 
*! optionaly subtract one 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM,'-',JUNIT,IDCUM, 1D0, 1D0) 
*! calculate difference with respect to NREF (MC-Ref)
      IF(IDREF.NE.0)
     $  CALL GOPERA(IDCUM,'-',IDREF,IDCUM, 1D0, 1D0) 
*! the difference MC-Analit. times XMAG
      IDIFF =  ID+500
      if(gexist(IDIFF)) write(6,*) ' +++PLDIFF: Warning! IDIFF = ',IDIFF
      CALL GOPERA(IDCUM,'-',IDAN,IDIFF, XMAG, XMAG)
*! magnify
      IDAN2 = IDAN+1000
      if(gexist(IDAN2)) write(6,*) ' +++PLDIFF: Warning! IDAN2 = ',IDIFF
      CALL GOPERA(IDCUM,'*',JUNIT,IDCUM, 1D0, YMAG)  
      CALL GOPERA(IDAN ,'*',JUNIT,IDAN2, 1D0, YMAG)  
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GMINIM(IDCUM,    YMIN)
      CALL GMAXIM(IDCUM,    YMAX)
      CALL GIDOPT(IDCUM,   'ERRO')
      CALL GPRINT(IDCUM)
      CALL GMINIM(IDAN2,    YMIN)
      CALL GMAXIM(IDAN2,    YMAX)
      CALL GPRINT(IDAN2)
      CALL GMINIM(IDIFF,    YMIN)
      CALL GMAXIM(IDIFF,    YMAX)
      CALL GIDOPT(IDIFF,   'ERRO')
      CALL GPRINT(IDIFF)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCUM,   ' ',' ',0)
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDAN2,   'S','*',0)
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDIFF   ,'S','*',0)
c purge local histos
      CALL GDELET(JUNIT)
      CALL GDELET(IDCUM)
      CALL GDELET(IDIFF)
      CALL GDELET(IDAN2)
      RETURN
 901  WRITE(6,*) ' PLDIFF: wrong born!!!!'
      END


      SUBROUTINE PLD2MC
     $   (TITLE,ID1,ID2,YMIN,YMAX,XMAG,YMAG,BORN1)
*     ***************************************************************
!---------------------------------------------------------------------
!  TITLE      title of the plot
!  ID1        histogram MC (not cumulative)
!  ID2        histogram MC (not cumulative)
!  YMIN,YMAX  uper lower limit in vertical scale of the plot
!  XMAG       magnification for difference ID2-ID1.
!  YMAG       magnification for ID1 and ID2.
!---------------------------------------------------------------------
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE
      LOGICAL gexist

      CALL GPLTIT(TITLE)
*! turn into cumulative normalized distribution or renormalize
      IDCU1  = ID1+2000
      IDCU2  = ID2+2000
      if(gexist(IDCU1)) write(6,*) ' ++++PLDIMC: Warning! IDCU1= ',IDCU1
      if(gexist(IDCU2)) write(6,*) ' ++++PLDIMC: Warning! IDCU2= ',IDCU2
*! divide by Born
      JUNIT=909090
      IF(BORN1.EQ.0D0) goto 901
      CALL HISUNI(ID1,JUNIT)
      CALL GOPERA(ID1,'/',JUNIT,IDCU1, 1D0, abs(BORN1)) 
      CALL GOPERA(ID2,'/',JUNIT,IDCU2, 1D0, abs(BORN1)) 
*! optionaly subtract one 
      IF(BORN1.LT.0D0) THEN
        CALL GOPERA(IDCU1,'-',JUNIT,IDCU1, 1D0, 1D0) 
        CALL GOPERA(IDCU2,'-',JUNIT,IDCU2, 1D0, 1D0)
      ENDIF 
*! the difference  ID1-ID2 times XMAG
      IDIFF =  ID1+200
      if(gexist(IDIFF)) write(6,*) ' ++++PLDIMC: Warning! IDIFF= ',IDIFF
      CALL GOPERA(IDCU1,'-',IDCU2,IDIFF, XMAG, XMAG)
*! magnify
      CALL GOPERA(IDCU1,'*',JUNIT,IDCU1, 1D0, YMAG)  
      CALL GOPERA(IDCU2,'*',JUNIT,IDCU2, 1D0, YMAG)  
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GMINIM(IDCU1,    YMIN)
      CALL GMAXIM(IDCU1,    YMAX)
      CALL GIDOPT(IDCU1,   'ERRO')
      CALL GPRINT(IDCU1)
      CALL GMINIM(IDCU2,    YMIN)
      CALL GMAXIM(IDCU2,    YMAX)
      CALL GPRINT(IDCU2)
      CALL GMINIM(IDIFF,    YMIN)
      CALL GMAXIM(IDIFF,    YMAX)
      CALL GIDOPT(IDIFF,   'ERRO')
      CALL GPRINT(IDIFF)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCU1,   ' ',' ',0)
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCU2,   'S','*',0)
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDIFF   ,'S','*',0)
c purge local histos
      CALL GDELET(JUNIT)
      CALL GDELET(IDCU1)
      CALL GDELET(IDCU2)
      CALL GDELET(IDIFF)
      RETURN
 901  WRITE(6,*) ' PLDIMC: wrong born!!!!'
      END

      SUBROUTINE PLC2MC
     $   (TITLE,ID1,ID2,IDIFF,YMIN,YMAX,XMAG,YMAG,BORN1)
*     ***************************************************************
!---------------------------------------------------------------------
* The same as PLD2MC but the difference is kept
!  TITLE      title of the plot
!  ID1        histogram MC (not cumulative)
!  ID2        histogram MC (not cumulative)
!  YMIN,YMAX  uper lower limit in vertical scale of the plot
!  XMAG       magnification for difference ID2-ID1.
!  YMAG       magnification for ID1 and ID2.
!---------------------------------------------------------------------
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE
      LOGICAL gexist

      CALL GPLTIT(TITLE)
*! turn into cumulative normalized distribution or renormalize
      IDCU1  = ID1+2000
      IDCU2  = ID2+2000
      if(gexist(IDCU1)) write(6,*) ' ++++PLDIMC: Warning! IDCU1= ',IDCU1
      if(gexist(IDCU2)) write(6,*) ' ++++PLDIMC: Warning! IDCU2= ',IDCU2
*! divide by Born
      JUNIT=909090
      IF(BORN1.EQ.0D0) goto 901
      CALL HISUNI(ID1,JUNIT)
      CALL GOPERA(ID1,'/',JUNIT,IDCU1, 1D0, abs(BORN1)) 
      CALL GOPERA(ID2,'/',JUNIT,IDCU2, 1D0, abs(BORN1)) 
*! optionaly subtract one 
      IF(BORN1.LT.0D0) THEN
        CALL GOPERA(IDCU1,'-',JUNIT,IDCU1, 1D0, 1D0) 
        CALL GOPERA(IDCU2,'-',JUNIT,IDCU2, 1D0, 1D0)
      ENDIF 
*! the difference  ID1-ID2 times XMAG
      if(gexist(IDIFF)) write(6,*) '++++PLC2IMC: Warning! IDIFF= ',IDIFF
      CALL GOPERA(IDCU1,'-',IDCU2,IDIFF, XMAG, XMAG)
*! magnify
      CALL GOPERA(IDCU1,'*',JUNIT,IDCU1, 1D0, YMAG)  
      CALL GOPERA(IDCU2,'*',JUNIT,IDCU2, 1D0, YMAG)  
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GMINIM(IDCU1,    YMIN)
      CALL GMAXIM(IDCU1,    YMAX)
      CALL GIDOPT(IDCU1,   'ERRO')
      CALL GPRINT(IDCU1)
      CALL GMINIM(IDCU2,    YMIN)
      CALL GMAXIM(IDCU2,    YMAX)
      CALL GPRINT(IDCU2)
      CALL GMINIM(IDIFF,    YMIN)
      CALL GMAXIM(IDIFF,    YMAX)
      CALL GIDOPT(IDIFF,   'ERRO')
      CALL GPRINT(IDIFF)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCU1,   ' ',' ',0)
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCU2,   'S','*',0)
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDIFF   ,'S','*',0)
c purge local histos
      CALL GDELET(JUNIT)
      CALL GDELET(IDCU1)
      CALL GDELET(IDCU2)
      RETURN
 901  WRITE(6,*) ' PLC2IMC: wrong born!!!!'
      END


      SUBROUTINE PL1MC(TITLE,ID,YMIN,YMAX,BORN1)
*     ******************************************
*! One Monte Carlo histo 
*!  TITLE      title of the plot
*!  ID         histogram MC cumulative, normalized
*!  YMIN,YMAX  uper lower limit in vertical scale of the plot
*!  Born       content divided by Born
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE

      BORN=ABS(BORN1)
      CALL GPLTIT(TITLE)
      IDCUM = ID+2000
*! divide by Born
      JUNIT=909090
      CALL HISUNI(ID,JUNIT)
      CALL GOPERA(ID,'/',JUNIT,IDCUM, 1D0, abs(BORN1))  
*! optionaly subtract one 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM,'-',JUNIT,IDCUM, 1D0, 1D0) 
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
      CALL GMINIM(IDCUM,    YMIN)
      CALL GMAXIM(IDCUM,    YMAX)
      CALL GIDOPT(IDCUM,   'ERRO')
      CALL GPRINT(IDCUM)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDCUM,   ' ','*',0)
c purge local histos
      CALL GDELET(JUNIT)
      CALL GDELET(IDCUM)
      END

      SUBROUTINE PL2MC(TITLE,ID1,ID2,YMIN,YMAX,BORN1,BORN2)
*     ****************************************************
*! Two Monte Carlo histos
*!  TITLE      title of the plot
*!  ID1,ID2    histograms MC, cumulative, normalized
*!  YMIN,YMAX  uper lower limit in vertical scale of the plot
*!  Born       content divided by Born
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE

      BORN=ABS(BORN1)
      CALL GPLTIT(TITLE)
      IDCUM1 = ID1+2000
      IDCUM2 = ID2+2000
*! divide by Born
      JUNIT=909090
      CALL HISUNI(ID1,JUNIT)
      CALL GOPERA(ID1,'/',JUNIT,IDCUM1, 1D0, abs(BORN1))  
      CALL GOPERA(ID2,'/',JUNIT,IDCUM2, 1D0, abs(BORN2))  
*! optionaly subtract one 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM1,'-',JUNIT,IDCUM1, 1D0, 1D0) 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM2,'-',JUNIT,IDCUM2, 1D0, 1D0) 
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
      CALL GMINIM(IDCUM1,    YMIN)
      CALL GMAXIM(IDCUM1,    YMAX)
      CALL GIDOPT(IDCUM1,   'ERRO')
      CALL GPRINT(IDCUM1)
      CALL GMINIM(IDCUM2,    YMIN)
      CALL GMAXIM(IDCUM2,    YMAX)
      CALL GIDOPT(IDCUM2,   'ERRO')
      CALL GPRINT(IDCUM2)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDCUM1,   ' ','*',0)
      CALL GPLSET('DMOD',4D0)
      CALL GPLOT(IDCUM2,   'S','*',0)
c purge local histos
      CALL GDELET(JUNIT)
      CALL GDELET(IDCUM1)
      CALL GDELET(IDCUM2)
      END



      SUBROUTINE PL3MC(TITLE,ID1,ID2,ID3,YMIN,YMAX,BORN1,BORN2,BORN3)
*     ***************************************************************
! Plots Three Monte Carlo histos
!  TITLE          title of the plot
!  ID1,ID2,ID3    histograms MC, cumulative, normalized
!  YMIN,YMAX      uper lower limit in vertical scale of the plot
!  Born           content divided by abs(Born)
!                 If Born<0, abs(Born) subtracted
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      CHARACTER*80 TITLE

      CALL GPLTIT(TITLE)
*! turn into cumulative normalized distribution
      IDCUM1 = ID1+2000
      IDCUM2 = ID2+2000
      IDCUM3 = ID3+2000
*! divide by Born
      JUNIT=909090
      CALL HISUNI(ID1,JUNIT)
      CALL GOPERA(ID1,'/',JUNIT,IDCUM1, 1D0, abs(BORN1))  
      CALL GOPERA(ID2,'/',JUNIT,IDCUM2, 1D0, abs(BORN2))  
      CALL GOPERA(ID3,'/',JUNIT,IDCUM3, 1D0, abs(BORN3))  
*! optionaly subtract one 
      IF(BORN1.LT.0D0)
     $  CALL GOPERA(IDCUM1,'-',JUNIT,IDCUM1, 1D0, 1D0) 
      IF(BORN2.LT.0D0)
     $  CALL GOPERA(IDCUM2,'-',JUNIT,IDCUM2, 1D0, 1D0) 
      IF(BORN3.LT.0D0)
     $  CALL GOPERA(IDCUM3,'-',JUNIT,IDCUM3, 1D0, 1D0) 
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
      CALL GMINIM(IDCUM1,    YMIN)
      CALL GMAXIM(IDCUM1,    YMAX)
      CALL GIDOPT(IDCUM1,   'ERRO')
      CALL GPRINT(IDCUM1)
      CALL GMINIM(IDCUM2,    YMIN)
      CALL GMAXIM(IDCUM2,    YMAX)
      CALL GIDOPT(IDCUM2,   'ERRO')
      CALL GPRINT(IDCUM2)
      CALL GMINIM(IDCUM3,    YMIN)
      CALL GMAXIM(IDCUM3,    YMAX)
      CALL GIDOPT(IDCUM3,   'ERRO')
      CALL GPRINT(IDCUM3)
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GPLSET('DMOD',1D0)
      CALL GPLOT(IDCUM1,   ' ','*',0)
      CALL GPLSET('DMOD',2D0)
      CALL GPLOT(IDCUM2,   'S','*',0)
      CALL GPLSET('DMOD',3D0)
      CALL GPLOT(IDCUM3,   'S','*',0)
c purge local histos
      CALL GDELET(JUNIT)
      CALL GDELET(IDCUM1)
      CALL GDELET(IDCUM2)
      CALL GDELET(IDCUM3)
      END



      SUBROUTINE CUMHYB(ID1,ID2,IDCUM)
C     ********************************
! Renormalizes and cumulates two histos ID1,ID2
! ID1 produced by OLDBIS and ID2 produces by LUMLOG
! Total result is stored in IDCUM
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  

      KeyGen = 1
      CALL CUMHIS(KeyGen,ID1,ID1+1000)
      KeyGen = 2
      CALL CUMHIS(KeyGen,ID2,ID2+1000)
      CALL GOPERA(ID1+1000,'+',ID2+1000,IDCUM, 1D0, 1D0)
      CALL GDELET(ID1+1000)
      CALL GDELET(ID2+1000)
      CALL GIDOPT(IDCUM,'ERRO')
      CALL GPRINT(IDCUM)
      END



*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

