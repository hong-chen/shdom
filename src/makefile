# Makefile for SHDOM and related programs in distribution

useOMP = no
# pick one of the following two
#╭────────────────────────────────────────────────╮#
FC = gfortran
#FC = mpif90
#╰────────────────────────────────────────────────╯#

FORTRAN_VERSION = $(shell $(FC) -dumpversion)
VERSION_LIMIT = 10
ifeq ($(shell test $(FORTRAN_VERSION) -lt $(VERSION_LIMIT); echo $$?), 0)
	FFLAGS = -O3 -Wargument-mismatch -ffree-line-length-none
else
	FFLAGS = -O3 -fallow-argument-mismatch -ffree-line-length-none
endif

CURRENT_DIR = $(shell pwd)
BIN_DIR     = $(abspath $(CURRENT_DIR)/../build/bin)
TARGETS = bin shdom shdom_ propgen make_mie_table make_tmatrix_table plotscattab put

ifeq ($(useOMP), yes)
	# openMPI
	FFLAGS += -fopenmp
	SHDOMOBJS = shdSolv_omp.o shdComp_omp.o shdIO_inp.o shdIO_out.o shdSurf_brdf.o fftpack.o globals.o
else
	SHDOMOBJS = shdSolv.o     shdComp.o     shdIO_inp.o shdIO_out.o shdSurf_brdf.o fftpack.o globals.o
endif

ifeq ($(FC), mpif90)
	UseMPI = yes
else
	UseMPI = no
endif

ifeq ($(UseMPI),no)
	SHDOMMPIOBJ = shdom_nompi.o
else
	SHDOMMPIOBJ = shdom_mpi.o
endif

SHDOMOBJS += $(SHDOMMPIOBJ)

F90 = $(FC)
F90FLAGS = $(FFLAGS)
LIBS =

all: $(TARGETS)

bin:
	mkdir -p $(BIN_DIR)

# this for namelist input
shdom: shdom.o $(SHDOMOBJS)
	$(F90) $(F90FLAGS) shdom.o $(SHDOMOBJS) $(LIBS) -o $(BIN_DIR)/shdom

shdom.o: shdom.f90 $(SHDOMOBJS)
	$(F90) -c $(F90FLAGS) shdom.f90

# this for user input
shdom_: shdom_.o $(SHDOMOBJS)
	$(F90) $(F90FLAGS) shdom_.o $(SHDOMOBJS) $(LIBS) -o $(BIN_DIR)/shdom_

shdom_.o: shdom_.f90 $(SHDOMOBJS)
	$(F90) -c $(F90FLAGS) shdom_.f90

shdSolv.o: shdSolv.f90 $(SHDOMMPIOBJ)
	$(F90) -c $(F90FLAGS) shdSolv.f90

shdSolv_omp.o: dev/shdSolv_omp.f90 $(SHDOMMPIOBJ)
	$(F90) -c $(F90FLAGS) dev/shdSolv_omp.f90

shdComp.o: shdComp.f90 globals.o $(SHDOMMPIOBJ)
	$(F90) -c $(F90FLAGS) shdComp.f90

shdComp_omp.o: dev/shdComp_omp.f90 globals.o $(SHDOMMPIOBJ)
	$(F90) -c $(F90FLAGS) dev/shdComp_omp.f90

shdIO_inp.o: shdIO_inp.f90 globals.o
	$(F90) -c $(F90FLAGS) shdIO_inp.f90

shdIO_out.o: shdIO_out.f90 globals.o
	$(F90) -c $(F90FLAGS) shdIO_out.f90

shdSurf_brdf.o: shdSurf_brdf.f90
	$(F90) -c $(F90FLAGS) shdSurf_brdf.f90

shdom_nompi.o: shdom_nompi.f90
	$(F90) -c $(F90FLAGS) shdom_nompi.f90

shdom_mpi.o: shdom_mpi.f90
	$(F90) -c $(F90FLAGS) shdom_mpi.f90

fftpack.o: fftpack.f90
	$(F90) -c $(F90FLAGS) fftpack.f90

globals.o: globals.f90
	$(F90) -c $(F90FLAGS) globals.f90

propgen: tools/shdTool_propgen.f90
	$(F90) $(F90FLAGS) tools/shdTool_propgen.f90 -o $(BIN_DIR)/propgen

make_tmatrix_table: tools/shdTool_make_tmx_tab.f90
	$(F90) $(F90FLAGS) tools/shdTool_make_tmx_tab.f90 tools/shdScat_ind_wat.f90 tools/shdScat_tmx_wig.f90 tools/_lapack.f -o $(BIN_DIR)/make_tmatrix_table

make_mie_table: tools/shdTool_make_mie_tab.f90
	$(F90) $(F90FLAGS) tools/shdTool_make_mie_tab.f90 tools/shdScat_ind_wat.f90 tools/shdScat_mie_wig.f90 -o $(BIN_DIR)/make_mie_table

plotscattab: tools/shdTool_plot_sca_tab.f90
	$(F90) $(F90FLAGS) tools/shdTool_plot_sca_tab.f90 -o $(BIN_DIR)/plotscattab

put: tools/shdTool_put.f90
	$(F90) $(F90FLAGS) tools/shdTool_put.f90 -o $(BIN_DIR)/put

.f.o :; $(FC) -c $(FFLAGS) $*.f

clean:
	@ rm -f *.o *.mod $(TARGETS) *~
