#================================================================
INCLUDES =  -I$(top_srcdir)/6xx -I$(ROOTINCDIR) $(all_includes)
METASOURCES = AUTO
AM_FFLAGS = @MY_F77_FLAGS@
#---------------------------------------------------------------
######### ROOT specific part  ##############
dictC = BHL6Dict.cxx
dictH = BHL6Dict.h
CLEANFILES = $(dictC) $(dictH)
# ROOT Directory generation:
noinst_HEADERS  = LinkDef.h
include_HEADERS = TMCgenBHL6.h
ROOTCINTHDRS= $(include_HEADERS) $(noinst_HEADERS)
$(dictC) : $(ROOTCINTHDRS)
	$(ROOTCINT)  -f $(dictC) -c $(INCLUDES) $(ROOTCINTHDRS)
#---------------------------------------------------------------
# part interpreted by automake, edited by automake manager
lib_LTLIBRARIES = libBHL6.la
nodist_libBHL6_la_SOURCES =$(dictC)
libBHL6_la_SOURCES = TMCgenBHL6.cxx

######## BHLUMI4 f77 libraries #############
$(top_builddir)/400/lib400.a:
	(cd $(top_builddir)/400/; make)
$(top_builddir)/gLIB/libGlib.a:
	(cd $(top_builddir)/gLIB/; make)
BHL4_Libs= \
	$(top_builddir)/400/lib400.a  $(top_builddir)/gLIB/libGlib.a 
############################################

##########################################################
######## Executables for testing ROOT persistency ########
bin_PROGRAMS = TestMini1 TestMini2 TestMini3
LOADLIST = \
	$(top_builddir)/6xx/libBHL6.la \
	$(BHL4_Libs) \
	$(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran -lFoam $(FLIBS)
#/////////////////////////////////////////////////////////
TestMini1_SOURCES = TestMini1.cxx 
TestMini1_LDADD   = $(LOADLIST)
TestMini1_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#/////////////////////////////////////////////////////////
TestMini2_SOURCES = TestMini2.cxx 
TestMini2_LDADD   = $(LOADLIST)
TestMini2_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#/////////////////////////////////////////////////////////
TestMini3_SOURCES = TestMini3.cxx 
TestMini3_LDADD   = $(LOADLIST)
TestMini3_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#/////////////////////////////////////////////////////////
test1-start: TestMini1
	(cd ../ make)
	./TestMini1
test23-start: TestMini2 TestMini3
	(cd ../ make)
	(./TestMini2; ./TestMini3 )
#/////////////////////////////////////////////////////////
