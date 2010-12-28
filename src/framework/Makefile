.SUFFIXES: .F .o

ifdef ZOLTAN_HOME
   ZOLTANOBJ = module_zoltan_interface.o
endif

OBJS = module_mpas_framework.o \
       module_timer.o \
       module_configure.o \
       module_constants.o \
       module_grid_types.o \
       module_hash.o \
       module_sort.o \
       module_block_decomp.o \
       module_dmpar.o \
       module_io_input.o \
       module_io_output.o \
       $(ZOLTANOBJ) \
       streams.o

all: framework

framework: $(OBJS)
	ar -ru libframework.a $(OBJS)

module_mpas_framework.o: module_dmpar.o module_io_input.o module_io_output.o module_grid_types.o module_configure.o module_timer.o

module_configure.o: module_dmpar.o

module_grid_types.o: module_dmpar.o

module_dmpar.o: module_sort.o streams.o

module_block_decomp.o: module_grid_types.o module_hash.o module_configure.o

module_io_input.o: module_grid_types.o module_dmpar.o module_block_decomp.o module_sort.o module_configure.o $(ZOLTANOBJ)

module_io_output.o: module_grid_types.o module_dmpar.o module_sort.o module_configure.o

clean:
	$(RM) *.o *.mod *.f90 libframework.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES)

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -c $<
