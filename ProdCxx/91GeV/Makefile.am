#=======================================================================
# Invokes all kind of scripts which setup (Start.C) and executes MAIN
#    make start; make stop
# Seting up input files/scripts and runing on the NZ42 PC-farm
#    make sfarm6; make ssubmit6; make q-sem; make q-nev
#    make farm-stop; make combine
# Farming under SLURM batch system on CC1
#    make SLfarm6; make SLsubmitall; ...
#=======================================================================
MAIN=MainProd
DSET=91GeV
## Local qubCMC=qsub, CLASS= qunlimited
subCMD=qsub
CLASS=qunlimitted
##############
check_all:
	(cd ../; make)
start_n: check_all
	rm -f ./pro.output ./bhl_f77.output./histo.root ./mcgen.root
	$(ROOTEXEC)  -b -q -l ./Start.C
start:  start_n
	time ../${MAIN} &
stop:
	$(ROOTEXEC)  -b -q -l ../farming/Stop.C
cont:
	time ../${MAIN} &
debug:  check_all
	ddd --directory="../" ../../bin/MainPr
#########################################################################
####      NEW UNIVERSAL VERSION tested on the local farm
#########################################################################
####   version for SLURM on NZ42 cluster ####
qfarm3:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmQ.C("$(DSET)",3)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
qfarm40:	farm-clean start_n
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmQ.C("$(DSET)",40)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
qsubmit6:
	($(ROOTEXEC)  -b -q -l '../farming/SubmFarmQ.C("$(subCMD)","$(CLASS)",6)' )
qqsubmitall:
	($(ROOTEXEC)  -b -q -l '../farming/SubmFarmQ.C("$(subCMD)","$(CLASS)",200)' )
#########################################################################
sfarm6:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmS.C("$(DSET)",6)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
sfarm24:	farm-clean start_n
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmS.C("$(DSET)",24)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
sfarm40:	farm-clean start_n
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmS.C("$(DSET)",40)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
sfarm80:	farm-clean start_n
	(echo ">>>Set up MC gener:"; $(ROOTEXEC) -b -q -l ./Start.C)
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmS.C("$(DSET)",80)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
ssubmitall:
	($(ROOTEXEC)  -b -q -l '../farming/SubmFarmS.C("sbatch",200)' )
#########################################################################
####   version for CC1 cluster ####
SLfarm3:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmSL.C("$(DSET)",3)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLfarm24:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmSL.C("$(DSET)",24)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLfarm82:	farm-clean start_n
	(echo ">>>Set up farm dir:"; $(ROOTEXEC) -b -q -l '../farming/SetFarmSL.C("$(DSET)",82)' ;)
	(ln -s ../../${MAIN} ./farm/${DSET}.exe; echo ">>>>>>>>... DONE!!!!!!!")
SLsubmitall:
	($(ROOTEXEC)  -b -q -l '../farming/SubmFarmS.C("sbatch",200)' )
########################################################################
farm-gtar:
	(gtar -cvzf ../../../${DSET}-farm.tar.gz farm)
#-----------------------------------------------------------------------
farm-stop:
	($(ROOTEXEC)  -b -q -l ../farming/StopFarm.C)
#--------- cobmine (add) results from several subdirectories -----------
combine:
#	($(ROOTSYS)/bin/hadd -f ./histo.root ./farm/*/histo.root)
	(hadd -f ./histo.root ./farm/*/histo.root)
	(cp -p ./farm/1/mcgen.root ./)
#-----------------------------------------------------------------------
#---- query
q-sem:
	(ls -altr ./farm/*/mcgen.root)
q-nev:
	($(ROOTEXEC)  -b -q -l ../farming/q-nev.C)
#-----------------------------------------------------------------------
farm-clean:
	(rm  -R -f farm *.bak )
#========================================================================
#========================================================================
