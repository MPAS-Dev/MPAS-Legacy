.SUFFIXES: .F .o

OBJS = module_RBF_interpolation.o module_vector_reconstruction.o module_spline_interpolation.o

all: operators

operators: $(OBJS)
	ar -ru libops.a $(OBJS)

module_vector_reconstruction.o: module_RBF_interpolation.o
module_RBF_interpolation.o:
module_spline_interpolation:

clean:
	$(RM) *.o *.mod *.f90 libops.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework
