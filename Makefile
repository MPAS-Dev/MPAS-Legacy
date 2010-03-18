#MODEL_FORMULATION = -DNCAR_FORMULATION
MODEL_FORMULATION = -DLANL_FORMULATION

EXPAND_LEVELS = -DEXPAND_LEVELS=26
#FILE_OFFSET = -DOFFSET64BIT

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
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) $(EXPAND_LEVELS) -D_MPI $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )
 
pgi:
	( make all \
	"FC = mpif90" \
	"CC = mpicc" \
	"SFC = pgf90" \
	"SCC = pgcc" \
	"FFLAGS = -r8 -O3" \
	"CFLAGS = -O3" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) $(EXPAND_LEVELS) -D_MPI -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

pgi-serial:
	( make all \
	"FC = pgf90" \
	"CC = pgcc" \
	"SFC = pgf90" \
	"SCC = pgcc" \
	"FFLAGS = -r8 -O0 -g -Mbounds -Mchkptr" \
	"CFLAGS = -O0 -g" \
	"LDFLAGS = -O0 -g -Mbounds -Mchkptr" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) $(EXPAND_LEVELS) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

ifort:
	( make all \
	"FC = mpif90" \
	"CC = gcc" \
	"SFC = ifort" \
	"SCC = gcc" \
	"FFLAGS = -real-size 64 -O3" \
	"CFLAGS = -O3 -m64" \
	"LDFLAGS = -O3" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) $(EXPAND_LEVELS) -D_MPI -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

gfortran:
	( make all \
	"FC = mpif90" \
	"CC = mpicc" \
	"SFC = gfortran" \
	"SCC = gcc" \
	"FFLAGS = -O3 -m64 -ffree-line-length-none" \
	"CFLAGS = -O3 -m64" \
	"LDFLAGS = -O3 -m64" \
	"CORE = $(CORE)" \
	"CPPFLAGS = -DRKIND=8 $(MODEL_FORMULATION) $(EXPAND_LEVELS) -D_MPI -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )



CPPINCLUDES = -I../inc -I$(NETCDF)/include
FCINCLUDES = -I../inc -I$(NETCDF)/include
LIBS = -L$(NETCDF)/lib -lnetcdf -L$(BLAS)/lib -L$(LAPACK)/lib -lblas -llapack

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