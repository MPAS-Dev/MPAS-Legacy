.SUFFIXES: .F .c .o

all: mpas

mpas: reg_includes externals frame ops dycore drver
	$(FC) $(LDFLAGS) -o $(CORE)_model.exe driver/*.o -L. -ldycore -lops -lframework $(LIBS) -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time

reg_includes: 
	( cd registry; make CC="$(SCC)" )
	( cd inc; ../registry/parse ../core_$(CORE)/Registry )

externals:
	( cd external; make FC="$(FC)" SFC="$(SFC)" CC="$(CC)" SCC="$(SCC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)" CPP="$(CPP)" NETCDF="$(NETCDF)" CORE="$(CORE)" )

frame: 
	( cd framework; make all ) 
	ln -sf framework/libframework.a libframework.a

ops: 
	( cd operators; make all ) 
	ln -sf operators/libops.a libops.a

dycore: 
	( cd core_$(CORE); make all ) 
	ln -sf core_$(CORE)/libdycore.a libdycore.a

drver: 
	( cd driver; make all ) 

clean:
	$(RM) $(CORE)_model.exe libframework.a libops.a libdycore.a
	( cd registry; make clean )
	( cd external; make clean )
	( cd framework; make clean )
	( cd operators; make clean )
	( cd inc; rm -f *.inc )
	if [ -d core_$(CORE) ] ; then \
	   ( cd core_$(CORE); make clean ) \
	fi;
	( cd driver; make clean )
