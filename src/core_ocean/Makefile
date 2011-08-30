.SUFFIXES: .F .o

OBJS = module_mpas_core.o \
       module_test_cases.o \
       module_advection.o \
       module_time_integration.o \
       module_global_diagnostics.o

all: core_hyd

core_hyd: $(OBJS)
	ar -ru libdycore.a $(OBJS)

module_test_cases.o:

module_advection.o:

module_time_integration.o: 

module_global_diagnostics.o: 

module_mpas_core.o: module_advection.o module_global_diagnostics.o module_test_cases.o module_time_integration.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
