AUTOMAKE_OPTIONS = foreign 1.6
SUBDIRS =   gLIB 400 6xx 4xCPC ProdCxx ProdCxxPlots
### Subdir below is for limited distribution like v4.22m 
### ( make distclean;  autoreconf -i --force ; ./configure) plus uncommenting in Makefile.am
###SUBDIRS = ???

# Automake 1.9.1 manual suggtests putting all local-macrros in m4 subdir
# and using this:
ACLOCAL_AMFLAGS = -I m4

#Trick: we want to distribute file called "makefile" in /work1.
#We want to force the use of "Makefile" (instead of "makefile", which
#has a precedence) in the whole installation process.
AM_MAKEFLAGS = -f Makefile

DATE = 2018-11
VERS = v5.01
