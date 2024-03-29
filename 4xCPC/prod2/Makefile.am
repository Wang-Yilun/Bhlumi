#-----------------------------------------------------------------------
#       Template for makefile in production directory prod2
#-----------------------------------------------------------------------
MAIN=demo2
DSET=prod2
CLASS=long
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Interactive:  run start stop and continue
#-----------------------------------------------------------------------
start:
	(cp ../iniseed/iniseed.0 ./iniseed     ;\
	 cp ../semaphore.start   ./semaphore   ;\
	 time ../${MAIN}  < ./$(DSET).data )&
#-----------------------------------------------------------------------
stop:
	(cp ../semaphore.stop ./semaphore )&
#-----------------------------------------------------------------------
cont:
	(cp semaphore  semaphore.bak           ;\
	 cp bhl.hst    bhl.hst.bak             ;\
	 cp bhl.output bhl.output.bak          ;\
         time ../${MAIN}  < ./$(DSET).data 	)&
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#   Batch on farm:  initialize, submit and stop
#-----------------------------------------------------------------------
farm12-SP:
	(csh ../farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MAIN} ${DSET} SP ${CLASS} )
farm24-SP:
	(csh ../farm-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	      13 14 15 16 17 18 19 20 21 22 23 24 " \
	${MAIN} ${DSET} SP ${CLASS} )
farm6-CSF:
	(csh ../farm-init "1 2 3 4 5 6" \
	${MAIN} ${DSET} CSF ${CLASS} )
farm12-CSF:
	(csh ../farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MAIN} ${DSET} CSF ${CLASS} )
farm6-BASTA:
	(csh ../farm-init "1 2 3 4 5 6" \
	${MAIN} ${DSET} BASTA ${CLASS} )
farm12-BASTA:
	(csh ../farm-init "1 2 3 4 5 6 7 8 9 10 11 12" \
	${MAIN} ${DSET} BASTA ${CLASS} )
farm18-BASTA:
	(csh ../farm-init  "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18" \
	${MAIN} ${DSET} BASTA ${CLASS} )
#-----------------------------------------------------------------------
farm-gtar:
	(gtar -cvzf ../../../${DSET}-farm.tar.gz farm)
#-----------------------------------------------------------------------
submit6:
	(csh ../farm-submit ${DSET} "1 2 3 4 5 6 ")
submit12:
	(csh ../farm-submit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12")
submit18:
	(csh ../farm-submit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 ")
llsubmit6:
	(csh ../farm-llsubmit ${DSET} "1 2 3 4 5 6")
llsubmit12:
	(csh ../farm-llsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12")
llsubmit18:
	(csh ../farm-llsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18")
llsubmit24:
	(csh ../farm-llsubmit ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#-----------------------------------------------------------------------
stop6:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 ")
stop12:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12")
stop18:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 ")
stop24:
	(csh ../farm-stop   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#
cphst24:
	(csh ../farm-cphst   ${DSET} "1 2 3 4 5 6 7 8 9 10 11 12 \
	                13 14 15 16 17 18 19 20 21 22 23 24 ")
#-----------------------------------------------------------------------
#--------- cobmine (add) results from several subdirectories -----------
combine:
	(cd ../../lib; make combine.exe)
combine6:  combine
	(sed -e 's#?.hst#bhl.hst#' -e 's#?#6#' ../combine.list \
	| ../../lib/combine.exe; more bhl.hst)
combine12:  combine
	(sed -e 's#?.hst#bhl.hst#' -e 's#?#12#' ../combine.list \
	| ../../lib/combine.exe; more bhl.hst)
combine18:  combine
	(sed -e 's#?.hst#bhl.hst#' -e 's#?#18#' ../combine.list \
	| ../../lib/combine.exe; more bhl.hst)
combine24:  combine
	(sed -e 's#?.hst#bhl.hst#' -e 's#?#24#' ../combine.list \
	| ../../lib/combine.exe; more bhl.hst)
#-----------------------------------------------------------------------
#---- query
q-semaphore:
	(ls -altr farm/*/semaphore)
q-nev:
	(cat farm/*/semaphore | grep NEV)
q-cont:
	(cat farm/*/semaphore | grep CONT)
#-----------------------------------------------------------------------
clean:                                                     
	(rm -f *.o *.a *~ *.exe core *.err *.out)
CLEAN:  clean
	(rm  -R -f farm )
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Miscelaneous
#-----------------------------------------------------------------------
script-SP:
	(csh ../farm-script-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	 13 14 15 16 17 18 19 20 21 22 23 24 " ${MAIN} ${DSET} SP ${CLASS} )
script-SP850:
	(csh ../farm-script-init "1 2 3 4 5 6 7 8 9 10 11 12 \
	 13 14 15 16 17 18 19 20 21 22 23 24 " ${MAIN} ${DSET} SP850 ${CLASS} )
#-----------------------------------------------------------------------