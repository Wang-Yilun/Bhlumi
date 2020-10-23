INCLUDES = -I$(top_srcdir)/400
METASOURCES = AUTO

lib_LIBRARIES = lib400.a
lib400_a_SOURCES = bhllog.f bhlum4.f bhlumi.f lumlog.f m2agzi.f modl2a.f \
	modl2b.f oldbis.f repi.f  model.f TRIG.f hadr5x19.f vp_DHMZ.f vp_knt_v3_0_1.f
	
lib400_a_LIBADD = $(top_builddir)/gLIB/libGlib.a

AM_FFLAGS = @MY_F77_FLAGS@

