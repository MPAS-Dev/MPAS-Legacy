# NETCDF should point to the directory holding the netcdf lib and bin subdirs

CC = gcc
INCLUDES = -I$(NETCDF)/include
LIBS = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf -L/usr/lib -lstdc++

transdual: nc2vtk_dual.cpp
	$(CC) $(INCLUDES) -o $@ $<  $(LIBS)

clean:
	rm transdual
