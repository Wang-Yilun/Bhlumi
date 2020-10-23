INCLUDES = -I$(top_srcdir)/4xCPC/figs

METASOURCES = AUTO
bin_PROGRAMS = PlotProd2

PlotProd2_SOURCES = PlotProd2.f plthis.f readat.f
PlotProd2_LDADD = $(top_builddir)/400/lib400.a $(top_builddir)/gLIB/libGlib.a

AM_FFLAGS = @MY_F77_FLAGS@

XDVI   = xdvi      -bg white -fg black  -s 4 
XPS    = gv -bg white -fg black  -magstep -1 
#XPS    = gv


#----------------------------------------------------------------------
Fp2 = PlotProd2
$(Fp2).tex: 	$(Fp2) ../prod2/bhl.hst
	(time   ./$(Fp2) )
$(Fp2)-dvi: $(Fp2).tex#
	(latex	$(Fp2).tex; 	$(XDVI) $(Fp2).dvi)
$(Fp2)-ps: $(Fp2).tex
	(latex	$(Fp2).tex; dvips $(Fp2).dvi -o; $(XPS) $(Fp2).ps)
#----------------------------------------------------------------------
