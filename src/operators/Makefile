.SUFFIXES: .F .o

OBJS = module_vector_reconstruction.o

all: operators

operators: $(OBJS)
	ar -ru libops.a $(OBJS)

module_vector_reconstruction.o:

clean:
	$(RM) *.o *.mod *.f90 libops.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework
