SHELL   = /bin/bash
CMP     = ifort # gfortran 
LIBS    = -mkl #-llapack -lblas #${LINK_LAPACK} #-mkl #-llapack -lblas
FRAMEWORK = # -framework Accelerate
FLAGS   = -O3 #-check all #-no-wrap-margin #-fbounds-check
DEBUG   = #-g -O0 #-ftrapuv #-check all -traceback #-O3
OBJrun  = test1.o constants.o

run.x: test1.o
	${CMP} -o run.x ${OBJrun} ${LIBS} ${FRAMEWORK} ${OMP} ${DEBUG} ${FLAGS}

test1.o: test1.f90 constants.o
	${CMP} -c ${DEBUG}  $<

constants.o: constants.f90 
	${CMP} -c ${DEBUG}  $<

clean:
	rm -f *.mod *.o *~ *.exe *.x *.dat

core: run.x
	ulimit -c unlimited; ./run.x	
