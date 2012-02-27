#MODEL_FORMULATION = -DNCAR_FORMULATION
MODEL_FORMULATION = -DLANL_FORMULATION

# This flag must be off for nersc hopper:
FILE_OFFSET = -DOFFSET64BIT

#########################
# Section for Zoltan TPL
#########################
ifdef ZOLTAN_HOME
   ZOLTAN_DEFINE = -DHAVE_ZOLTAN
endif
#########################


dummy:
	( make error )

xlf:
	( make all \
	"FC_PARALLEL = mpxlf90" \
	"CC_PARALLEL = mpcc" \
	"FC_SERIAL = xlf90" \
	"CC_SERIAL = xlc" \
	"FFLAGS_OPT = -O3 -qrealsize=8" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )
 
ftn:
	( make all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -i4 -r8 -gopt -O2 -Mvect=nosse -Kieee -convert big_endian" \
	"CFLAGS_OPT = -fast" \
	"LDFLAGS_OPT = " \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

pgi:
	( make all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = pgf90" \
	"CC_SERIAL = pgcc" \
	"FFLAGS_OPT = -r8 -O3 -byteswapio -Mfree" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -r8 -O0 -g -Mbounds -Mchkptr -byteswapio -Mfree" \
	"CFLAGS_DEBUG = -O0 -g" \
	"LDFLAGS_DEBUG = -O0 -g -Mbounds -Mchkptr" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

pgi-nersc:
	( make all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -r8 -O3 -byteswapio -Mfree" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

pgi-llnl:
	( make all \
	"FC_PARALLEL = mpipgf90" \
	"CC_PARALLEL = pgcc" \
	"FC_SERIAL = pgf90" \
	"CC_SERIAL = pgcc" \
	"FFLAGS_OPT = -i4 -r8 -g -O2 -byteswapio" \
	"CFLAGS_OPT = -fast" \
	"LDFLAGS_OPT = " \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

ifort:
	( make all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = gcc" \
	"FC_SERIAL = ifort" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -real-size 64 -O3 -convert big_endian -FR" \
	"CFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3" \
	"FFLAGS_DEBUG = -real-size 64 -g -convert big_endian -FR -CU -CB -check all" \
	"CFLAGS_DEBUG = -g -m64" \
	"LDFLAGS_DEBUG = -g" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

gfortran:
	( make all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = gfortran" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -O3 -m64 -ffree-line-length-none -fdefault-real-8 -fconvert=big-endian -ffree-form" \
	"CFLAGS_OPT = -O3 -m64" \
	"LDFLAGS_OPT = -O3 -m64" \
	"FFLAGS_DEBUG = -g -m64 -ffree-line-length-none -fdefault-real-8 -fconvert=big-endian -ffree-form -fbounds-check" \
	"CFLAGS_DEBUG = -g -m64" \
	"LDFLAGS_DEBUG = -g -m64" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE -m64 $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

g95:
	( make all \
	"FC_PARALLEL = mpif90" \
	"CC_PARALLEL = mpicc" \
	"FC_SERIAL = g95" \
	"CC_SERIAL = gcc" \
	"FFLAGS_OPT = -O3 -ffree-line-length-huge -r8 -fendian=big" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

pathscale-nersc:
	( make all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -r8 -O3 -freeform -extend-source" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

cray-nersc:
	( make all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -default64 -O3 -f free" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

intel-nersc:
	( make all \
	"FC_PARALLEL = ftn" \
	"CC_PARALLEL = cc" \
	"FC_SERIAL = ftn" \
	"CC_SERIAL = cc" \
	"FFLAGS_OPT = -real-size 64 -O3 -FR" \
	"CFLAGS_OPT = -O3" \
	"LDFLAGS_OPT = -O3" \
	"CORE = $(CORE)" \
	"DEBUG = $(DEBUG)" \
	"SERIAL = $(SERIAL)" \
	"USE_PAPI = $(USE_PAPI)" \
	"CPPFLAGS = $(MODEL_FORMULATION) -DUNDERSCORE $(FILE_OFFSET) $(ZOLTAN_DEFINE)" )

CPPINCLUDES = -I../inc -I$(NETCDF)/include
FCINCLUDES = -I../inc -I$(NETCDF)/include
LIBS = -L$(NETCDF)/lib -lnetcdf

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

ifeq "$(DEBUG)" "true"

ifndef FFLAGS_DEBUG
	FFLAGS=$(FFLAGS_OPT)
	CFLAGS=$(CFLAGS_OPT)
	LDFLAGS=$(LDFLAGS_OPT)
	DEBUG_MESSAGE="Debug flags are not defined for this compile group. Defaulting to Optimized flags"
else # FFLAGS_DEBUG IF
	FFLAGS=$(FFLAGS_DEBUG)
	CFLAGS=$(CFLAGS_DEBUG)
	LDFLAGS=$(LDFLAGS_DEBUG)
	DEBUG_MESSAGE="Debugging is on."
endif # FFLAGS_DEBUG IF

else # DEBUG IF
	FFLAGS=$(FFLAGS_OPT)
	CFLAGS=$(CFLAGS_OPT)
	LDFLAGS=$(LDFLAGS_OPT)
	DEBUG_MESSAGE="Debugging is off."
endif # DEBUG IF

ifeq "$(SERIAL)" "true"
	FC=$(FC_SERIAL)
	CC=$(CC_SERIAL)
	SFC=$(FC_SERIAL)
	SCC=$(CC_SERIAL)
	SERIAL_MESSAGE="Serial version is on."
else # SERIAL IF
	FC=$(FC_PARALLEL)
	CC=$(CC_PARALLEL)
	SFC=$(FC_SERIAL)
	SCC=$(CC_SERIAL)
	CPPFLAGS += -D_MPI
	SERIAL_MESSAGE="Parallel version is on."
endif # SERIAL IF

ifeq "$(USE_PAPI)" "true"
	CPPINCLUDES += -I$(PAPI)/include -D_PAPI
	FCINCLUDES += -I$(PAPI)/include
	LIBS += -L$(PAPI)/lib -lpapi
	PAPI_MESSAGE="Papi libraries are on."
else # USE_PAPI IF
	PAPI_MESSAGE="Papi libraries are off."
endif # USE_PAPI IF

all: mpas_main

mpas_main: 
	cd src; make FC="$(FC)" \
                 CC="$(CC)" \
                 SFC="$(SFC)" \
                 SCC="$(SCC)" \
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
	@echo ""
	@echo $(DEBUG_MESSAGE)
	@echo $(SERIAL_MESSAGE)
	@echo $(PAPI_MESSAGE)
clean:
	cd src; make clean RM="$(RM)" CORE="$(CORE)"
	$(RM) $(CORE)_model.exe
error: errmsg

else # CORE IF

all: error
clean: errmsg
error: errmsg
	@echo "************ ERROR ************"
	@echo "No CORE specified. Quitting."
	@echo "************ ERROR ************"
	@echo ""

endif # CORE IF

errmsg:
	@echo ""
	@echo "Usage: make target CORE=[core] [options]"
	@echo ""
	@echo "Example targets:"
	@echo "    ifort"
	@echo "    gfortran"
	@echo "    xlf"
	@echo "    pgi"
	@echo ""
	@echo "Availabe Cores:"
	@cd src; ls -d core_* | grep ".*" | sed "s/core_/    /g"
	@echo ""
	@echo "Available Options:"
	@echo "    SERIAL=true - builds serial version. Default is parallel version."
	@echo "    DEBUG=true  - builds debug version. Default is optimized version."
	@echo "    USE_PAPI=true   - builds version using PAPI for timers and hardware counters. Default is off."
	@echo ""
	@echo "Ensure that NETCDF (and PAPI if USE_PAPI=true) are environment variables"
	@echo "that point to the absolute paths for the libraries."
	@echo ""

