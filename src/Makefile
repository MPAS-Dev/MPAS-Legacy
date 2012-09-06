.SUFFIXES: .F .c .o

all: mpas

mpas: reg_includes externals frame ops dycore drver
	$(LINKER) $(LDFLAGS) -o $(CORE)_model.exe driver/*.o -L. -ldycore -lops -lframework $(LIBS) -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time

reg_includes: 
	( cd registry; $(MAKE) CC="$(SCC)" )
	( cd inc; $(CPP) ../core_$(CORE)/Registry | ../registry/parse > Registry.processed)

externals: reg_includes
	( cd external; $(MAKE) FC="$(FC)" SFC="$(SFC)" CC="$(CC)" SCC="$(SCC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)" CPP="$(CPP)" NETCDF="$(NETCDF)" CORE="$(CORE)" )

frame: reg_includes externals
	( cd framework; $(MAKE) all ) 
	ln -sf framework/libframework.a libframework.a

ops: reg_includes externals frame
	( cd operators; $(MAKE) all ) 
	ln -sf operators/libops.a libops.a

dycore: reg_includes externals frame ops
	( cd core_$(CORE); $(MAKE) all ) 
	ln -sf core_$(CORE)/libdycore.a libdycore.a

drver:  reg_includes externals frame ops dycore
	( cd driver; $(MAKE) all ) 

clean:
	$(RM) $(CORE)_model.exe libframework.a libops.a libdycore.a
	( cd registry; $(MAKE) clean )
	( cd external; $(MAKE) clean )
	( cd framework; $(MAKE) clean )
	( cd operators; $(MAKE) clean )
	( cd inc; rm -f *.inc Registry.processed )
	if [ -d core_$(CORE) ] ; then \
	   ( cd core_$(CORE); $(MAKE) clean ) \
	fi;
	( cd driver; $(MAKE) clean )
