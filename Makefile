SHELL   = /bin/bash
CMP     = ifort # gfortran 
LIBS    = -mkl #-llapack -lblas #${LINK_LAPACK} #-mkl #-llapack -lblas
FRAMEWORK = # -framework Accelerate
FLAGS   = -O1 -r8 #-check all #-no-wrap-margin #-fbounds-check
DEBUG   = #-g -O0 #-ftrapuv #-check all -traceback #-O3
HAN     = -fp-model precise
OBJrun  = dynamics_main.o constants.o


run.x: dynamics_main.o
	${CMP} -o run.x ${OBJrun} ${LIBS} ${FRAMEWORK} ${OMP} ${DEBUG} ${FLAGS} ${HAN}

dynamics_main.o: dynamics_main.f90 constants.o
	${CMP} -c ${DEBUG}   ${FLAGS} ${HAN} $<

constants.o: constants.f90 
	${CMP} -c ${DEBUG}  ${FLAGS} ${HAN} $<

clean:
	rm -f *.mod *.o *~ *.exe *.x *.dat

core: run.x
	ulimit -c unlimited; ./run.x	
