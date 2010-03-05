.SUFFIXES: .F .o

OBJS = module_test_cases.o \
       module_vector_reconstruction.o \
       module_time_integration.o \
       mpas_interface.o

all: core_sw

core_sw: $(OBJS)
	ar -ru libdycore.a $(OBJS)

module_test_cases.o: 

module_vector_reconstruction.o:

module_time_integration.o: module_vector_reconstruction.o

mpas_interface.o: module_test_cases.o module_time_integration.o module_vector_reconstruction.o

clean:
	$(RM) *.o *.mod libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework
	$(RM) $*.f90
