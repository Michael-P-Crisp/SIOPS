F90 = gfortran
COMPILER_OPTIONS = -O3 -ffree-line-length-none #-ffpe-trap=zero,overflow,underflow -g -Og -fbacktrace -fcheck=all
LINKER_OPTIONS =  

execname = soil

#SOURCES_F90 := $(shell find . -name '*.F90')

SOURCES_F90 = akima.F90 akima_ind.F90 akima_nocheck.F90 checkfiles.F90 corn2d.F90 corn3d.F90 daxpy.F90 dchol2.F90 DCVAA2.F90 dcvaa3.F90 DCVAB2.F90 dcvab3.F90 DCVIT2.F90 dcvit3.F90 DCVMT2.F90 dcvmt3.F90 ddot.F90 despile.F90 detcheck.F90 DLAFR2.F90 dlafs3.F90 DLAVX2.F90 dlavx3.F90 DLSEP2.F90 dlsep3.F90 DLSFR2.F90 dlsfr3.F90 DLSPX2.F90 dlspx3.F90 dot3d.F90 dot3db.F90 dsifa.F90 dsisl.F90 dswap.F90 edge3d.F90 edivide.F90 ESETT.F90 extrafuncs.F90 fem_2d.F90 fem_3d.F90 fem_geom.F90 fem_main.F90 fem_prep.F90 fem_stress.F90 fem_stuff.F90 ga.F90 gausv.F90 getboundary.F90 getcosts.F90 getdiff.F90 getheights.F90 getperfmulti.F90 idamax.F90 int2d.F90 intr2d.F90 intr3d.F90 iseed.F90 LAS1D.F90 LAS2G.F90 LAS2I.F90 las3g.F90 las3i.F90 Lnblnk.F90 MAIN.F90 output_results.F90 piecewise_CMD.F90 plan3d.F90 prepmulti.f90 PROCESS_CK.F90 PROCESS_SI.F90 PROCESS_SI_multi.F90 pwl_interp_2d_scattered.F90 qsort_c_module.F90 r8lib.F90 randu.F90 READINM.F90 reducem.F90 SETUP_EA.F90 SETUP_SI.F90 SI.F90 side2d.F90 side3d.F90 SIm.F90 SIM2SD.F90 SIM2SD_init.F90 sim3de.F90 sim3de_init.F90 sim_det_rand.F90 si_stats.F90 soilgen.F90 soilstat.F90 thin1d.F90 variables.F90 vnorm.F90 weights.F90 WRITESOILS.F90
OBJECTS = $(SOURCES_F90:%.F90=%.o)

MODULES = sim3de.mod piecewise_CMD.mod LAS1D.mod variables.mod checkfiles.mod READINM.mod edivide.mod fem_geom.mod fem_main.mod fem_stuff.mod fem_prep.mod SETUP_SI.mod si_stats.mod output_results.mod  \
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

