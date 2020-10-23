INCLUDES = -I$(top_srcdir)/gLIB
METASOURCES = AUTO

lib_LIBRARIES = libGlib.a
libGlib_a_SOURCES = glibk.f yfslib.f

bin_PROGRAMS = combine
combine_SOURCES = combine.f
combine_LDADD = $(top_builddir)/gLIB/libGlib.a

AM_FFLAGS = @MY_F77_FLAGS@

