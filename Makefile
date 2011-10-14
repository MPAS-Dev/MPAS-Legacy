CORE=ocean
#MODEL_FORMULATION = -DNCAR_FORMULATION
MODEL_FORMULATION = -DLANL_FORMULATION

FILE_OFFSET = -DOFFSET64BIT

#########################
# Section for Zoltan TPL
#########################
ifdef ZOLTAN_HOME
   ZOLTAN_DEFINE = -DHAVE_ZOLTAN
endif
#########################


dummy:
	@( echo "try one of:"; \
	echo "   make xlf"; \
	echo "   make pgi"; \
	echo "   make ifort"; \
	echo "   make gfortran"; \
	)

xlf:
	( make all \
	"FC = mpxlf90" \
	"CC = mpcc" \
	"SFC = xlf90" \
	"SCC = xlc" \
	"FFLAGS = -qrealsize=8 -g -C " \
	"CFLAGS = -g" \
	"LDFLAGS = -g -C" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_MPI $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )
 
ftn:
	( make all \
	"FC = ftn" \
	"CC = cc" \
	"SFC = ftn" \
	"SCC = cc" \
	"FFLAGS = -i4 -r8 -gopt -O2 -Mvect=nosse -Kieee -convert big_endian" \
	"CFLAGS = -fast" \
	"LDFLAGS = " \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

pgi:
	( make all \
	"FC = mpif90" \
	"CC = mpicc" \
	"SFC = pgf90" \
	"SCC = pgcc" \
	"FFLAGS = -r8 -O3 -byteswapio -Mfree" \
	"CFLAGS = -O3" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

pgi-llnl:
	( make all \
	"FC = mpipgf90" \
	"CC = pgcc" \
	"SFC = pgf90" \
	"SCC = pgcc" \
	"FFLAGS = -i4 -r8 -g -O2 -byteswapio" \
	"CFLAGS = -fast" \
	"LDFLAGS = " \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

pgi-serial:
	( make all \
	"FC = pgf90" \
	"CC = pgcc" \
	"SFC = pgf90" \
	"SCC = pgcc" \
	"FFLAGS = -r8 -O0 -g -Mbounds -Mchkptr -byteswapio -Mfree" \
	"CFLAGS = -O0 -g" \
	"LDFLAGS = -O0 -g -Mbounds -Mchkptr" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

ifort-serial:
	( make all \
	"FC = ifort" \
	"CC = gcc" \
	"SFC = ifort" \
	"SCC = gcc" \
	"FFLAGS = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS = -O3 -m64" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

ifort-papi:
	( make all \
	"FC = mpif90" \
	"CC = gcc" \
	"SFC = ifort" \
	"SCC = gcc" \
	"FFLAGS = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS = -O3 -m64" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_PAPI -D_MPI -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" \
	"PAPILIBS = -L$(PAPI)/lib -lpapi" )

ifort-papi-serial:
	( make all \
	"FC = ifort" \
	"CC = gcc" \
	"SFC = ifort" \
	"SCC = gcc" \
	"FFLAGS = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS = -O3 -m64" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_PAPI -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" \
	"PAPILIBS = -L$(PAPI)/lib -lpapi" )

ifort:
	( make all \
	"FC = mpif90" \
	"CC = gcc" \
	"SFC = ifort" \
	"SCC = gcc" \
	"FFLAGS = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS = -O3 -m64" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

gfortran:
	( make all \
	"FC = mpif90" \
	"CC = mpicc" \
	"SFC = gfortran" \
	"SCC = gcc" \
	"FFLAGS = -O3 -m64 -ffree-line-length-none -fdefault-real-8 -fconvert=big-endian -ffree-form" \
	"CFLAGS = -O3 -m64" \
	"LDFLAGS = -O3 -m64" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

gfortran-serial:
	( make all \
	"FC = gfortran" \
	"CC = gcc" \
	"SFC = gfortran" \
	"SCC = gcc" \
	"FFLAGS = -O3 -m64 -ffree-line-length-none -fdefault-real-8 -fconvert=big-endian -ffree-form" \
	"CFLAGS = -O3 -m64" \
	"LDFLAGS = -O3 -m64" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) $(EXPAND_LEVELS) -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

g95:
	( make all \
	"FC = mpif90" \
	"CC = mpicc" \
	"SFC = g95" \
	"SCC = gcc" \
	"FFLAGS = -O3 -ffree-line-length-huge -r8 -fendian=big" \
	"CFLAGS = -O3" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -D_MPI -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

g95-serial:
	( make all \
	"FC = g95" \
	"CC = gcc" \
	"SFC = g95" \
	"SCC = gcc" \
	"FFLAGS = -O3 -ffree-line-length-huge -r8 -fendian=big" \
	"CFLAGS = -O3" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )


CPPINCLUDES = -I../inc -I$(NETCDF)/include -I$(PAPI)/include
FCINCLUDES = -I../inc -I$(NETCDF)/include -I$(PAPI)/include
LIBS = -L$(NETCDF)/lib -lnetcdf $(PAPILIBS)

RM = rm -f
CPP = cpp -C -P -traditional
RANLIB = ranlib

#########################
# Section for Zoltan TPL
#########################
ifdef ZOLTAN_HOME
   ifdef ZOLTAN_INC_PATH
      FCINCLUDES += -I$(ZOLTAN_INC_PATH)
   else
      FCINCLUDES += -I$(ZOLTAN_HOME)/include
   endif

   ifdef ZOLTAN_LIB_PATH
      LIBS += -L$(ZOLTAN_LIB_PATH) -lzoltan
   else
      LIBS += -L$(ZOLTAN_HOME)/lib -lzoltan
   endif
endif
#########################


ifdef CORE

all: mpas_main

mpas_main: 
	cd src; make FC="$(FC)" \
                     CC="$(CC)" \
                     CFLAGS="$(CFLAGS)" \
                     FFLAGS="$(FFLAGS)" \
                     LDFLAGS="$(LDFLAGS)" \
                     RM="$(RM)" \
                     CPP="$(CPP)" \
                     CPPFLAGS="$(CPPFLAGS)" \
                     LIBS="$(LIBS)" \
                     CPPINCLUDES="$(CPPINCLUDES)" \
                     FCINCLUDES="$(FCINCLUDES)" \
                     CORE="$(CORE)"
	if [ ! -e $(CORE)_model.exe ]; then ln -s src/$(CORE)_model.exe .; fi

clean:
	cd src; make clean RM="$(RM)" CORE="$(CORE)"
	$(RM) $(CORE)_model.exe

else

all: errmsg
clean: errmsg
errmsg:
	@echo "************ ERROR ************"
	@echo "No CORE specified. Quitting."
	@echo "************ ERROR ************"

endif
