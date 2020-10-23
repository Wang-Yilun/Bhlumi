INCLUDES = -I$(top_srcdir)/6xx  -I$(top_srcdir)/ProdCxx $(all_includes) -I$(ROOTINCDIR)
METASOURCES = AUTO

bin_PROGRAMS = Plot1 Plot2 Plot3

######## BHLUMI4 f77 libraries #############
$(top_builddir)/400/lib400.a:
	(cd $(top_builddir)/400/; make)
$(top_builddir)/gLIB/libGlib.a:
	(cd $(top_builddir)/gLIB/; make)
BHL4_Libs= \
	$(top_builddir)/400/lib400.a  $(top_builddir)/gLIB/libGlib.a 
############################################
LOADLIST = \
	$(top_builddir)/6xx/libBHL6.la \
	$(BHL4_Libs) \
	$(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran $(FLIBS)
#-----------------------------------------------------------------------

Plot1_SOURCES = HSTplot.cxx HisNorm.cxx Plot1.cxx
Plot1_LDADD   = $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran $(LOADLIST)
Plot1_LDFLAGS = -L$(ROOTLIBDIR) -no-install -g

Plot2_SOURCES = HSTplot.cxx HisNorm.cxx Plot2.cxx
Plot2_LDADD   = $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran $(LOADLIST)
Plot2_LDFLAGS = -L$(ROOTLIBDIR) -no-install -g

Plot3_SOURCES = HSTplot.cxx HisNorm.cxx Plot3.cxx
Plot3_LDADD   = $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran $(LOADLIST)
Plot3_LDFLAGS = -L$(ROOTLIBDIR) -no-install -g

######################### private #############################
SUFFIXES:  .ps .dvi .tex .txp .html .bib .blg .bbl .aux .fig .pdf .eps .jpg
#==============================================================
#==============================================================
.txp.eps:
	(latex $*.txp; dvips -E $*.dvi -o $*.eps;)
.txp.pdf:
	( latex $*.txp; dvips -E $*.dvi -o $*.eps ; ps2pdf -dEPSCrop $*.eps ;)
#	($(FIG2PDF)  $*.fig $*.pdf )


#################
plot1: Plot1
	./Plot1
plot1-run: Plot1
	./Plot1
#################
plot3: Plot3
	./Plot3
################# old 1996 benchmark ##########
####### make OldBench-pdf
Plot1.phony: Plot1
	(./Plot1; touch Plot1.phony;)
BM1 = OldBench
TabOldBench.txp:   	Plot1.phony
TabVP1.txp:        	Plot1.phony
TabVP2.txp:        	Plot1.phony
TabVP3.txp:        	Plot1.phony
#
TabDiff.txp:        Plot2.phony
TabDiff2.txp:       Plot2.phony
TabDiff3.txp:       Plot2.phony
######
$(BM1).pdf: $(BM1).tex TabOldBench.pdf TabVP1.pdf TabVP2.pdf TabVP3.pdf
	pdflatex $(BM1).tex
$(BM1)-pdf: $(BM1).pdf
	evince  $(BM1).pdf
########
TabVP1-pdf: Plot1.phony TabVP1.pdf
	(xpdf TabVP1.pdf )
########
TabVP3-pdf:  Plot1.phony TabVP3.pdf
	(xpdf TabVP3.pdf )
###############################################
############### LEP VP study ##################
# make Tab91GeV-pdf
Tab88GeV.txp:		Plot2.phony
Tab89GeV.txp:		Plot2.phony
Tab90GeV.txp:		Plot2.phony
Tab91GeV.txp:		Plot2.phony
Tab92GeV.txp:		Plot2.phony
Tab93GeV.txp:		Plot2.phony
Tab94GeV.txp:		Plot2.phony
Plot2.phony: Plot2
	(./Plot2; touch Plot2.phony;)
Tab91GeV-pdf: Tab91GeV.pdf
	xpdf Tab91GeV.pdf
Tab93GeV-pdf: Tab93GeV.pdf
	xpdf Tab93GeV.pdf
Tab89GeV-pdf: Tab89GeV.pdf
	xpdf Tab89GeV.pdf
TabDiff-pdf: TabDiff.pdf
	xpdf TabDiff.pdf
####### make NewVPs-pdf
PDFlist2= Tab88GeV.pdf Tab89GeV.pdf Tab90GeV.pdf Tab91GeV.pdf Tab92GeV.pdf Tab93GeV.pdf Tab94GeV.pdf TabDiff.pdf TabDiff2.pdf TabDiff3.pdf
BM2 = NewVPs
$(BM2).tex:	           $(PDFlist2)
$(BM2).pdf: $(BM2).tex $(PDFlist2)
	pdflatex $(BM2).tex
$(BM2)-pdf: $(BM2).pdf
	evince  $(BM2).pdf
############### End ofLEP VP study ############
###############################################
############### Z study #######################
# make Zstudy-pdf
TabZ_88GeV.txp:    Plot3.phony
TabZ_89GeV.txp:    Plot3.phony
TabZ_90GeV.txp:    Plot3.phony
TabZ_91GeV.txp:    Plot3.phony
TabZ_92GeV.txp:    Plot3.phony
TabZ_93GeV.txp:    Plot3.phony
TabZ_94GeV.txp:    Plot3.phony
TabZet1.txp:       Plot3.phony
TabZet2.txp:       Plot3.phony
TabZdel.txp:       Plot3.phony

Plot3.phony: Plot3
	(./Plot3; touch Plot3.phony;)
####### make Zstudy-pdf
PDFlist3= TabZ_88GeV.pdf TabZ_89GeV.pdf TabZ_90GeV.pdf TabZ_91GeV.pdf TabZ_92GeV.pdf TabZ_93GeV.pdf TabZ_94GeV.pdf \
          TabZet1.pdf TabZet2.pdf TabZdel.pdf
BM3 = Zstudy
$(BM3).tex:	            $(PDFlist3)
$(BM3).pdf:  $(BM3).tex $(PDFlist3)
	pdflatex $(BM3).tex
$(BM3)-pdf:  $(BM3).pdf
	evince   $(BM3).pdf

###############################################################
klean:
	(rm -f *.o *.exe *~ *.dvi *.log *.aux dump.hst core lmap *.phony)
	(rm -f *.txp chi afb output-* *.hst)
###############################################################
