INCLUDES = -I$(top_srcdir)/4xCPC
SUBDIRS =  prod2 figs
METASOURCES = AUTO
bin_PROGRAMS = demo demo2

demo_SOURCES = demo.f
demo_LDADD = $(top_builddir)/400/lib400.a $(top_builddir)/gLIB/libGlib.a

demo2_SOURCES = demo2.f  trical.f trigas.f silicon.f
demo2_LDADD = $(top_builddir)/400/lib400.a $(top_builddir)/gLIB/libGlib.a

AM_FFLAGS = @MY_F77_FLAGS@

#-----------------------------------------------------------------------
#                   Demo simple test
#-----------------------------------------------------------------------
demo-start: demo
	(cd demo-work;   time ../demo )
	(cd demo-work;   latex demo.tex;  dvips demo.dvi -o; gv demo.ps )

#-----------------------------------------------------------------------
#                   Demo2 benchmarks
#-----------------------------------------------------------------------
prod2-start:  demo2
	(cd prod2;  make start)
prod2-stop:
	(cd prod2;  make stop)
prod2-cont:   demo2
	(cd prod2;  make cont)
prod2-ps:
	(cd figs;   make PlotProd2-ps)
#---------------------------

