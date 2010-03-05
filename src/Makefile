.SUFFIXES: .F .c .o

all: mpas

mpas: reg_includes externals frame dycore drver
	$(FC) $(LDFLAGS) -o $(CORE)_model.exe driver/*.o -L. -ldycore -lframework $(LIBS)

reg_includes: 
	( cd registry; make CC="$(SCC)" )
	( cd inc; ../registry/parse ../core_$(CORE)/Registry )

externals:
	( cd external; make FC="$(FC)" SFC="$(SFC)" CC="$(CC)" SCC="$(SCC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)" CPP="$(CPP)" RANLIB="$(RANLIB)" NETCDF="$(NETCDF)" CORE="$(CORE)" )

frame: 
	( cd framework; make all ) 
	ln -sf framework/libframework.a libframework.a

dycore: 
	( cd core_$(CORE); make all ) 
	ln -sf core_$(CORE)/libdycore.a libdycore.a

drver: 
	( cd driver; make all ) 

clean:
	$(RM) $(CORE)_model.exe libframework.a libdycore.a
	( cd registry; make clean )
	( cd external; make clean )
	( cd framework; make clean )
	( cd inc; rm -f *.inc )
	if [ -d core_$(CORE) ] ; then \
	   ( cd core_$(CORE); make clean ) \
	fi;
	( cd driver; make clean )
