.SUFFIXES: .F .o

OBJS = module_test_cases.o \
       module_time_integration.o \
       module_vector_reconstruction.o \
       module_global_diagnostics.o \
       mpas_interface.o

all: core_hyd

core_hyd: $(OBJS)
	ar -ru libdycore.a $(OBJS)

module_test_cases.o:

module_time_integration.o: module_vector_reconstruction.o

module_vector_reconstruction.o: 

module_global_diagnostics.o: 

mpas_interface.o: module_global_diagnostics.o module_vector_reconstruction.o module_test_cases.o module_time_integration.o

clean:
	$(RM) *.o *.mod libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework
	$(RM) $*.f90
