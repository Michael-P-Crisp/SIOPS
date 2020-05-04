F90 = gfortran
COMPILER_OPTIONS = -O3  -ffree-line-length-none #  -O0 -g -fbacktrace -fcheck=all  #-ffpe-trap=zero,overflow,underflow -g -fbacktrace -fcheck=all
LINKER_OPTIONS =  #-O0 -g -fbacktrace -fcheck=all -ffree-line-length-none #-ffpe-trap=zero,overflow,underflow -g -fbacktrace -ffast-math

execname = soil

SOURCES_F90 := $(shell find . -name '*.F90')
OBJECTS = $(SOURCES_F90:%.F90=%.o)

MODULES = sim3de.mod piecewise_CMD.mod LAS1D.mod variables.mod checkfiles.mod READINM.mod edivide.mod sim3de.mod fem_geom.mod fem_main.mod fem_stuff.mod fem_prep.mod SETUP_SI.mod si_stats.mod output_results.mod  \
 pwl_interp_2d_scattered.mod extrafuncs.mod getdiff.mod SIM2SD.mod qsort_c_module.mod int2d.mod soilgen.mod WRITESOILS.mod PROCESS_CK.mod SI.mod reducem.mod getperfmulti.mod ESETT.mod PROCESS_SI_multi.mod PROCESS_SI.mod ga.mod SETUP_EA.mod weights.mod

${execname}: ${MODULES} ${OBJECTS}
	${F90} ${OBJECTS} -o ${execname} ${LINKER_OPTIONS}

%.o: %.F90
	${F90} -c ${COMPILER_OPTIONS} $< -o $@

%.mod: %.o
	${F90} -c ${COMPILER_OPTIONS} $< -o $@

all: soil

clean:
	rm -fr *.mod ${execname} *.o

