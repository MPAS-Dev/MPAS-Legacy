# NETCDF should point to the directory holding the netcdf lib and bin subdirs

CC = gcc
INCLUDES = -I$(NETCDF)/include
LIBS = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf -L/usr/lib -lstdc++
#MYDEBUG = -g
MYDEBUG =

proj3dsphere: proj3dsphere.o convert.o
	$(CC) $(MYDEBUG) $(INCLUDES) -o $@ $< convert.o $(LIBS)

proj3dsphere.o: proj3dsphere.cpp
	$(CC) $(MYDEBUG) $(INCLUDES) -c $<

convert.o: convert.c
	$(CC) $(MYDEBUG) -c $<

clean:
	rm *.o proj3dsphere
