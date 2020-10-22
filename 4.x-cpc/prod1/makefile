#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#             Special  makefile for prod1
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#             Not from ../makefile.template!!! 
#             See for instance 'merge' in stream-diting of combine.list
#
# It is for running 10 parallel jobs with 10 different energies.
# Z contribution is tested.
#
#   make start		Run demo2.f at E=MZ
#   make stop		Stop the job
#   make cont		Continue the job
#
#
#   make farmSP		Create 10 subdirectiories each for diffrent energy
#   make llsubmit	Submitt 10 jobs
#   make combine10      Merge results into single file prod1.hst
#
#---------------------------------------------------------------------------
DSET=prod1
MAIN=demo2
#-----------------------------------------------------------------------
# Interactive:  run start stop and continue
#-----------------------------------------------------------------------
start:
	(cp ../iniseed/iniseed.0 ./iniseed     ;\
	 cp ../semaphore.start   ./semaphore   ;\
	 time ../${MAIN}.exe  < ./$(DSET).data )&
#-----------------------------------------------------------------------
stop:
	(cp ../semaphore.stop ./semaphore )&
#-----------------------------------------------------------------------
cont:
	(cp semaphore  semaphore.bak           ;\
	 cp bhl.hst    bhl.hst.bak             ;\
	 cp bhl.output bhl.output.bak          ;\
         time ../${MAIN}.exe  < ./$(DSET).data 	)&
#---------------------------------------------------------------------------
#            initialize, submit and stop
#---------------------------------------------------------------------------
farmSP:
	(csh ./farm-init   ${DSET} "1 2 3 4 5 6 7 8 9 10" SP850)
#
llsubmit:
	(csh ../farm-llsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10")
#
stop:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10")
#
#---------------------------------------------------------------------------
#---- cobmine results
combine:
	(cd ../../lib; make )
combine10:  combine
	(sed -e 's#?.hst#${DSET}.hst#'\
	-e 's#?#10#'\
	-e 's#add#merge#'\
	../combine.list | ../../lib/combine.exe;\
	more ${DSET}.hst)
#---------------------------------------------------------------------------
#---- query
q-semaphore:
	(ls -altr farm/*/semaphore)
q-nev:
	(cat farm/*/semaphore    | grep NEV)
q-cont:
	(cat farm/*/semaphore    | grep CONT)
q-ene:
	(cat farm/*/${DSET}.data | grep CMSENE)
#---------------------------------------------------------------------------
clean:                                                     
	(rm -f *.o *.a *~ *.exe core *.err *.out)
CLEAN:  clean
	(rm  -R farm )
#---------------------------------------------------------------------------
