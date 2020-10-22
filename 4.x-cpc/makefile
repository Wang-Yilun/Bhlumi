###########################################################################
# Collections of commands:
#
#   Simple demo
# 	make demo		basic demonstration run 100k events
#
#   Production runs demo2
#	make prod2-start	start   BHLUM4 run
#	make prod2-stop 	stop    BHLUM4 run
#	make prod2-cont 	restart BHLUM4 run
#	make prod2-dvi		view results (xdvi)
#
#	make llog2-start	start   LUMLOG run
#	make llog2-stop 	stop    LUMLOG run
#	make llog2-cont 	restart LUMLOG run
#	make llog2-dvi 		view results (xdvi)
#
#	make obis2-start	start   OLDBIS run
#	make obis2-stop 	stop    OLDBIS run
#	make obis2-cont 	restart OLDBIS run
#	make obis2-dvi 		view results (xdvi)
#
#   ======== Plots based on ACTUAL histograms ========
#
#    	make BeMaTabs-dvi	CPC table  benchmark
#    	make BeMaFigs-dvi	CPC figure benchmark
#
#	make makfig2-dvi	th/95-38 plots, Bhlumi-(Oldbis+Lumlog)
#	make makfig4-dvi	th/96-156 plots, Missing O(alf3)LL
#	make Fig-MisO3-ps       again th/96-156, Missing O(alf3)LL
#
#   ======== Plots based on ARCHIVE histograms ========
#
#    	make arch-BeMaTabs-dvi	CPC table  benchmark
#    	make arch-BeMaFigs-dvi	CPC figure benchmark
#
#	make arch-makfig1-dvi	Z-contribution, th/95-74 plot
#	make arch-makfig2-dvi	th/95-38 plots, Bhlumi-(Oldbis+Lumlog)
#	make arch-makfig4-dvi	th/96-156 plots, Missing O(alf3)LL
#	make arch-Fig-MisO3-ps	again th/96-156, Missing O(alf3)LL
#
#	make arch-makfig3-dvi		exercise on dummy parameters in MC
#    	make arch-BeMaTabs-100k-dvi	CPC figure benchmark, 100k events
#
#	For more see figs/makefile
###########################################################################
FFLAGS  =  -fno-automatic -ffixed-line-length-132 -fbackslash -C -O
LDFLAGS = 
F77 = gfortran
LD  = $(F77)
#-----------------------------------------------------------------------
.f.o:
	$(F77) $(FFLAGS) -c $<
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
EXTER_LIB1  = ../400/400.a
$(EXTER_LIB1): ../400/*.f
	cd ../400; make
#-----------------------------------------------------------------------
EXTER_LIB2  = ../lib/glib.a
$(EXTER_LIB2): ../lib/*.f
	cd ../lib; make
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
COMMAND = demo.exe
OBJECTS = demo.o
# Link procedure
$(COMMAND): $(OBJECTS)                     $(EXTER_LIB1) $(EXTER_LIB2)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS)  $(EXTER_LIB1) $(EXTER_LIB2)
#-----------------------------------------------------------------------
COMMAND2 = demo2.exe
OBJECTS2 = demo2.o trigas.o silicon.o trical.o 
$(COMMAND2): $(OBJECTS2)                   $(EXTER_LIB1) $(EXTER_LIB2)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS2) $(EXTER_LIB1) $(EXTER_LIB2)
#-----------------------------------------------------------------------
#                   Demo entries
#-----------------------------------------------------------------------
demo: $(COMMAND)
	(cd demo;   time ../demo.exe )
	(cd demo;   latex demo.tex;  xdvi demo; dvips demo.dvi -o )
demo-diff-rs:
	(cd demo; diff demo.output.100k.keypia.eq.2.keywgt.eq.0.rs demo.output)
	(cd demo; diff demo.tex.100k.keypia.eq.2.keywgt.eq.0.rs    demo.tex)
BeMaTabs-dvi:
	(cd figs; make BeMaTabs-dvi)
BeMaFigs-dvi:
	(cd figs; make BeMaFigs-dvi)
makfig2-dvi:
	(cd figs; make makfig2-dvi)
makfig4-dvi:
	(cd figs; make makfig4-dvi)
Fig-MisO3-ps:
	(cd figs; make Fig-MisO3-ps)
#-----------------------------------------------------------------------
arch-makfig1-dvi:
	(cp prod1/prod1.hst.2M.10ene  prod1/prod1.hst)
	(cd figs; make makfig1-dvi)
copy2arch:
	(cp prod2/prod2.data.2037M  prod2/prod2.data)
	(cp obis2/obis2.data.2118M  obis2/obis2.data)
	(cp llog2/llog2.data.2114M  llog2/llog2.data)
	(cp prod2/prod2.hst.2037M   prod2/bhl.hst)
	(cp obis2/obis2.hst.2118M   obis2/bhl.hst)
	(cp llog2/llog2.hst.2114M   llog2/bhl.hst)
arch-makfig2-dvi:	copy2arch
	(cd figs; make makfig2-dvi)
arch-BeMaTabs-dvi:	copy2arch
	(cd figs; make BeMaTabs-dvi)
arch-BeMaFigs-dvi:	copy2arch
	(cd figs; make BeMaFigs-dvi)
arch-BeMaTabs-100k-dvi:
	(cp prod2/prod2.hst.100k   prod2/bhl.hst)
	(cp obis2/obis2.hst.100k   obis2/bhl.hst)
	(cp llog2/llog2.hst.100k   llog2/bhl.hst)
	(cd figs; make BeMaTabs-dvi)
arch-makfig4-dvi:
	(cp llog2/llog2.data.537M   llog2/llog2.data)
	(cp llog2/llog2.hst.537M    llog2/bhl.hst)
	(cd figs; make makfig4-dvi)
arch-Fig-MisO3-ps: arch-makfig4-dvi
	(cd figs; make Fig-MisO3-ps)
arch-makfig3-dvi:
	(cp prod3/prod3.data.eps.e-5            prod3/prod3.data)
	(cp prod3/prod3.hst.2083M.eps.e-5       prod3/bhl.hst)
	(cd figs; make makfig3-dvi)
arch-bmakfig3b-dvi:
	(cp prod3/prod3.data.keywgt.eq.0        prod3/prod3.data)
	(cp prod3/prod3.hst.2090M.keywgt.eq.0   prod3/bhl.hst)
	(cd figs; make makfig3-dvi)
arch-makfig3c-dvi:
	(cp prod3/prod3.data.keywgt.eq.1        prod3/prod3.data)
	(cp prod3/prod3.hst.1973M.keywgt.eq.1   prod3/bhl.hst)
	(cd figs; make makfig3-dvi)
arch-makfig3d-dvi:
	(cp prod3/prod3.data.keywgt.eq.2        prod3/prod3.data)
	(cp prod3/prod3.hst.2028M.keywgt.eq.2   prod3/bhl.hst)
	(cd figs; make makfig3-dvi)
#-----------------------------------------------------------------------
#                   MC Production
#-----------------------------------------------------------------------
prod1-start:  $(COMMAND2)
	(cd prod1;  make start)
prod1-stop:
	(cd prod1;  make stop)
prod1-cont:  $(COMMAND2)
	(cd prod1;  make cont)
#---------------------------
prod2-start:  $(COMMAND2)
	(cd prod2;  make start)
prod2-stop:
	(cd prod2;  make stop)
prod2-cont:  $(COMMAND2)
	(cd prod2;  make cont)
prod2-dvi:
	(cd figs; make prod2-dvi)
#---------------------------
llog2-start:  $(COMMAND2)
	(cd llog2;  make start)
llog2-stop:
	(cd llog2;  make stop)
llog2-cont:  $(COMMAND2)
	(cd llog2;  make cont)
llog2-dvi:
	(cd figs; make llog2-dvi)
#---------------------------
obis2-start:   $(COMMAND2)
	(cd obis2;  make start)
obis2-stop:
	(cd obis2;  make stop)
obis2-cont:  $(COMMAND2)
	(cd obis2;  make cont)
obis2-dvi:
	(cd figs; make obis2-dvi)
#---------------------------
prod3-start:  $(COMMAND2)
	(cd prod3;  make start)
prod3-stop:
	(cd prod3;  make stop)
prod3-cont:  $(COMMAND2)
	(cd prod3;  make cont)
#-----------------------------------------------------------------------
#                  update production makefiles
#-----------------------------------------------------------------------
makprod:
#-----------------------------------
	(sed  -e 's#<main>#demo2#' -e 's#<dset>#prod2#'  -e 's#<class>#long#'\
		./makefile.template > ./prod2/makefile)
#-----------------------------------
	(sed  -e 's#<main>#demo2#' -e 's#<dset>#obis2#'  -e 's#<class>#long#'\
		./makefile.template > ./obis2/makefile)
#-----------------------------------
	(sed  -e 's#<main>#demo2#' -e 's#<dset>#llog2#'  -e 's#<class>#long#'\
		./makefile.template > ./llog2/makefile)
#-----------------------------------
	(sed  -e 's#<main>#demo2#' -e 's#<dset>#prod3#'  -e 's#<class>#long#'\
		./makefile.template > ./prod3/makefile)
#-----------------------------------------------------------------------
#        Freeze/unfreeze CVS directories
#-----------------------------------------------------------------------
cvs-freeze:
	(cd ../; find ./ -name CVS -exec chmod -R -w {} \; )
cvs-unfreeze:
	(cd ../; find ./ -name CVS -exec chmod -R +w {} \; )
#-----------------------------------------------------------------------
#                   Cleaning
#-----------------------------------------------------------------------
clean:
	(rm -f  *.o  *.a  *~  *.exe  core  *.gz)
	(cd figs; make clean)
Clean:  clean
	(cd ../lib; make clean)
	(cd ../400; make clean)
	(cd ../dok; make clean)
	(cd ./prod2;  make clean)
	(cd ./obis2;  make clean)
	(cd ./llog2;  make clean)
	(cd ./prod3;  make clean)
CLEAN:  Clean
	(cd ../dok;   make CLEAN)
	(cd ./prod2;  make CLEAN)
	(cd ./obis2;  make CLEAN)
	(cd ./llog2;  make CLEAN)
	(cd ./prod3;  make CLEAN)
#-----------------------------------------------------------------------
gtar: Clean
	(cp ./README ../ )
	(cd ../../;   gtar cvzf   bhlp-4x-prod.tar.gz \
			bhlp-4x/4.x-cpc \
			bhlp-4x/400 \
			bhlp-4x/lib )
	(cd ../../;   chmod -w    bhlp-4x-prod.tar.gz )
#-----------------------------------------------------------------------
#                          Export version
#-----------------------------------------------------------------------
# All steps are:   cvs rtag   -R <tag> bhlp-4.x-cpc
#                  cvs export -r <tag> bhlp-4.x-cpc
#                  cd bhlp-4.x-cpc/4.x-cpc
#                  make export
#-----------------------------------------------------------------------
export: CLEAN
	(cp ./README ../ )
	(cd ../../;   gtar cvzf   bhlp-4.x-export.tar.gz  bhlp-4.x-cpc)
	(cd ../../;   chmod -w    bhlp-4.x-export.tar.gz )
#-----------------------------------------------------------------------
EXPORT: 
	(cd ../../;  cvs rtag    -R v404m  bhlp-4.x-cpc )
	(cd ../../;  cvs export  -r v404m  bhlp-4.x-cpc )
	(cd ../../bhlp-4.x-cpc/4.x-cpc; make export )
#----------------------------------------------------------------------------
#### tags: v404a (12sept), v404b (13sept) v404c (16sept) v404d (17sept)
#### tags: v404e (18sept), v404f (20sept) v404g (20sept) v404h (20sept)
#### tags: v404i (20sept), v404j (21sept) v404k (21sept) v404l (23sept) 
#### tags: v404m (25sept) 