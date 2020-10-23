AC_INIT(configure.ac)
AM_INIT_AUTOMAKE(KKMCee, 4.24)
AC_PREFIX_DEFAULT(${PWD}/)
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_F77([gfortran f95])


AM_CONFIG_HEADER(config.h)

# Do compilation tests using the compiler, preprocessor, and file extensions
AC_LANG_FORTRAN77

AC_F77_LIBRARY_LDFLAGS

MY_F77_FLAGS="-O -g -fno-second-underscore -fno-automatic -ffixed-line-length-132 -fbackslash"
AC_SUBST(MY_F77_FLAGS)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Uncomment only for making distro to be used without automake
###AM_MAINTAINER_MODE
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ROOT_PATH(3.02/00, , AC_MSG_ERROR( ROOT has not been found.))
#HAS_CERNLIB( , AC_MSG_ERROR( CERNLIB has not been found))
AC_LANG(C)
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM([[#include<stdio.h> ]],
			[[const char hw[] = "Hello, World\n";]],
                       [[fputs (hw, stdout);]])],
	[AC_MSG_RESULT(TESTPROG OK)], [AC_MSG_RESULT(TestProg Failed)]
	);

##### macros defining different queue and class 
##### for NQS batch scripts
QUEUE_CMD
CLASS_CMD 

AC_OUTPUT(Makefile \
400/Makefile \
6xx/Makefile \
gLIB/Makefile \
4xCPC/Makefile \
4xCPC/prod2/Makefile \
4xCPC/figs/Makefile \
ProdCxx/Makefile \
ProdCxx/Prod2/Makefile \
ProdCxx/88GeV/Makefile \
ProdCxx/89GeV/Makefile \
ProdCxx/90GeV/Makefile \
ProdCxx/91GeV/Makefile \
ProdCxx/92GeV/Makefile \
ProdCxx/93GeV/Makefile \
ProdCxx/94GeV/Makefile \
ProdCxxPlots/Makefile \
)
