.SUFFIXES: .F .o

OBJS = module_mpas_subdriver.o \
       mpas.o

all: $(OBJS)

module_mpas_subdriver.o: 

mpas.o: module_mpas_subdriver.o

clean:
	$(RM) *.o *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../core_$(CORE) -I../external/esmf_time_f90
