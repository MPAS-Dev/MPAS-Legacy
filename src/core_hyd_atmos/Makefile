.SUFFIXES: .F .o

OBJS = mpas_atmh_mpas_core.o \
       mpas_atmh_test_cases.o \
       mpas_atmh_time_integration.o \
       mpas_atmh_advection.o

all: core_hyd

core_hyd: $(OBJS)
	ar -ru libdycore.a $(OBJS)

mpas_atmh_test_cases.o: 

mpas_atmh_time_integration.o: 

mpas_atmh_advection.o: 

mpas_atmh_mpas_core.o: mpas_atmh_advection.o mpas_atmh_test_cases.o mpas_atmh_time_integration.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
