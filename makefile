include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


CFLAGS += -pedantic -std=c99

: Posson1D.o
	-${CLINKER} -o Posson1D Posson1D.o  ${PETSC_LIB}
	${RM} Posson1D.o

distclean:
	@rm -f *~ Posson1D *tmp
