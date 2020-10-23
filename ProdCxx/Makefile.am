# ****** This is input data file for Automake ****** 
# Automake creates Makefile using all information defined below
INCLUDES =  -I$(top_srcdir)/6xx  -I$(top_srcdir)/ProdCxx -I$(ROOTINCDIR) $(all_includes)
METASOURCES = AUTO
SUBDIRS = Prod2 88GeV 89GeV 91GeV 92GeV 93GeV 94GeV
#///////////////////////////////////////////////////////////////
lib_LTLIBRARIES = libProd.la
CLEANFILES = $(dict1C) $(dict1H)
#///////////////////////////////////////////////////////////////
#   ROOT Dictionary 1 file definition
dict1C = ProdDict.cxx
dict1H = ProdDict.h
#  Header files:
include_HEADERS1 = TRobolProd.h TSemaf.h TRig.h
noinst_HEADERS1  = LinkDef1.h
#  ROOT Directory generation:
ROOTCINTHDRS1= $(include_HEADERS1) $(noinst_HEADERS1)
$(dict1C) : $(ROOTCINTHDRS1)
	$(ROOTCINT) -f $(dict1C) -c $(INCLUDES) $(ROOTCINTHDRS1)
#---------------------------------------------------------------
#  Local library 1
#  Additional flags:
nodist_libProd_la_SOURCES = $(dict1C)
libProd_la_SOURCES =  TRobolProd.cxx TSemaf.cxx TRig.cxx
libProd_la_LIBADD =  $(top_builddir)/6xx/libBHL6.la \
    $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS)
libProd_la_LDFLAGS = -L$(ROOTLIBDIR)
#---------------------------------------------------------------

######## BHLUMI4 f77 libraries #############
$(top_builddir)/400/lib400.a:
	(cd $(top_builddir)/400/; make)
$(top_builddir)/gLIB/libGlib.a:
	(cd $(top_builddir)/gLIB/; make)
BHL4_Libs= \
	$(top_builddir)/400/lib400.a  $(top_builddir)/gLIB/libGlib.a 
############################################

#///////////////////////////////////////////////////////////////
# Local executables with ROOT persistency
bin_PROGRAMS = MainProd
AM_FFLAGS = @MY_F77_FLAGS@
#-----------------------------------------------------------------------
LOADLIST = \
	$(top_builddir)/6xx/libBHL6.la \
	$(top_builddir)/ProdCxx/libProd.la \
	$(BHL4_Libs) \
	$(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran -lFoam $(FLIBS)
#-----------------------------------------------------------------------
MainProd_SOURCES = MainProd.cxx 
MainProd_LDADD   = $(LOADLIST)
MainProd_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#///////////////////////////////////////////////////////////////
#-----------------------------------------------------------------------
Prod2-start: all
	(cd ./Prod2; make start)
Prod2-stop:
	(cd ./Prod2; make stop)
#-----------------------------------------------------------------------
88GeV-start: all
	(cd ./88GeV; make start)
88GeV-stop:
	(cd ./88GeV; make stop)
#-----------------------------------------------------------------------
89GeV-start: all
	(cd ./89GeV; make start)
89GeV-stop:
	(cd ./89GeV; make stop)
#-----------------------------------------------------------------------
90GeV-start: all
	(cd ./90GeV; make start)
90GeV-stop:
	(cd ./90GeV; make stop)
#-----------------------------------------------------------------------
91GeV-start: all
	(cd ./91GeV; make start)
91GeV-stop:
	(cd ./91GeV; make stop)
#-----------------------------------------------------------------------
92GeV-start: all
	(cd ./92GeV; make start)
92GeV-stop:
	(cd ./92GeV; make stop)
#-----------------------------------------------------------------------
93GeV-start: all
	(cd ./93GeV; make start)
93GeV-stop:
	(cd ./93GeV; make stop)
#-----------------------------------------------------------------------
94GeV-start: all
	(cd ./94GeV; make start)
94GeV-stop:
	(cd ./94GeV; make stop)
#-----------------------------------------------------------------------
start-3: 89GeV-start 91GeV-start 93GeV-start
stop-3:  89GeV-stop  91GeV-stop  93GeV-stop
#-----------------------------------------------------------------------
start-4: 88GeV-start 90GeV-start 92GeV-start 94GeV-start
stop-4:  88GeV-stop  90GeV-stop  92GeV-stop  94GeV-stop
#-----------------------------------------------------------------------
SLfarm-all:
	(cd 88GeV; make SLfarm3; make SLsubmitall;)
	(cd 89GeV; make SLfarm3; make SLsubmitall;)
	(cd 90GeV; make SLfarm3; make SLsubmitall;)
	(cd 91GeV; make SLfarm3; make SLsubmitall;)
	(cd 92GeV; make SLfarm3; make SLsubmitall;)
	(cd 93GeV; make SLfarm3; make SLsubmitall;)
	(cd 94GeV; make SLfarm3; make SLsubmitall;)
#-----------------------------------------------------------------------
SLfarm-combine:
	(cd 88GeV; make combine;)
	(cd 89GeV; make combine;)
	(cd 90GeV; make combine;)
	(cd 91GeV; make combine;)
	(cd 92GeV; make combine;)
	(cd 93GeV; make combine;)
	(cd 94GeV; make combine;)
SLfarm-import:
	(cd 88GeV; scp jadach@192.245.169.227:~/work/GITSPACE/BHLUMI/bhlumi6/ProdCxx/88GeV/histo.root .)
	(cd 89GeV; scp jadach@192.245.169.227:~/work/GITSPACE/BHLUMI/bhlumi6/ProdCxx/89GeV/histo.root .)
	(cd 90GeV; scp jadach@192.245.169.227:~/work/GITSPACE/BHLUMI/bhlumi6/ProdCxx/90GeV/histo.root .)
	(cd 91GeV; scp jadach@192.245.169.227:~/work/GITSPACE/BHLUMI/bhlumi6/ProdCxx/91GeV/histo.root .)
	(cd 92GeV; scp jadach@192.245.169.227:~/work/GITSPACE/BHLUMI/bhlumi6/ProdCxx/92GeV/histo.root .)
	(cd 93GeV; scp jadach@192.245.169.227:~/work/GITSPACE/BHLUMI/bhlumi6/ProdCxx/93GeV/histo.root .)
	(cd 94GeV; scp jadach@192.245.169.227:~/work/GITSPACE/BHLUMI/bhlumi6/ProdCxx/94GeV/histo.root .)
SLfarm-stop:
	(cd 88GeV; make farm-stop)
	(cd 89GeV; make farm-stop)
	(cd 90GeV; make farm-stop)
	(cd 91GeV; make farm-stop)
	(cd 92GeV; make farm-stop)
	(cd 93GeV; make farm-stop)
	(cd 94GeV; make farm-stop)
#-----------------------------------------------------------------------
Plot1-run:
	(cd ../ProdCxxPlots; make Plot1; ./Plot1)
#-----------------------------------------------------------------------
