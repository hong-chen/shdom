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

# dummy option
ifeq ($(useOMP), yes)
	# openMPI
	FFLAGS += -fopenmp
endif
SHDOMOBJS = globals.o shdUtil.o shdSfc.o shdPha.o shdInit.o shdGrid_ind.o shdGrid_cal.o \
			shdDir.o shdSrc.o shdSHDO.o shdSolv_int.o shdSolv_rte.o shdComp.o shdIO_inp.o shdIO_out.o fftpack.o

ifeq ($(FC), mpif90)
	UseMPI = yes
else
	UseMPI = no
endif

ifeq ($(UseMPI),no)
	SHDOMMPIOBJ = shdExt_mpi_.o
else
	SHDOMMPIOBJ = shdExt_mpi.o
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

shdInit.o: shdInit.f90 globals.o shdGrid_ind.o shdUtil.o shdSfc.o
	$(F90) -c $(F90FLAGS) shdInit.f90

shdGrid_ind.o: shdGrid_ind.f90 globals.o
	$(F90) -c $(F90FLAGS) shdGrid_ind.f90

shdGrid_cal.o: shdGrid_cal.f90 globals.o shdUtil.o shdSrc.o shdDir.o shdGrid_ind.o $(SHDOMMPIOBJ)
	$(F90) -c $(F90FLAGS) shdGrid_cal.f90

shdDir.o: shdDir.f90 globals.o
	$(F90) -c $(F90FLAGS) shdDir.f90

shdSrc.o: shdSrc.f90 globals.o shdPha.o shdUtil.o shdGrid_ind.o
	$(F90) -c $(F90FLAGS) shdSrc.f90

shdSHDO.o: shdSHDO.f90 globals.o shdPha.o shdSfc.o
	$(F90) -c $(F90FLAGS) shdSHDO.f90

shdSolv_int.o: shdSolv_int.f90 $(SHDOMMPIOBJ) shdSfc.o shdInit.o shdGrid_ind.o shdSHDO.o
	$(F90) -c $(F90FLAGS) shdSolv_int.f90

shdSolv_rte.o: shdSolv_rte.f90 $(SHDOMMPIOBJ) shdSfc.o shdInit.o shdGrid_ind.o shdSHDO.o shdDir.o shdPha.o shdSrc.o shdGrid_cal.o shdSolv_int.o
	$(F90) -c $(F90FLAGS) shdSolv_rte.f90

shdComp.o: shdComp.f90 globals.o shdPha.o shdUtil.o shdSrc.o shdInit.o
	$(F90) -c $(F90FLAGS) shdComp.f90

shdIO_inp.o: shdIO_inp.f90 globals.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdIO_inp.f90

shdIO_out.o: shdIO_out.f90 globals.o shdGrid_ind.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdIO_out.f90

shdSfc.o: shdSfc.f90 globals.o
	$(F90) -c $(F90FLAGS) shdSfc.f90

shdPha.o: shdPha.f90 globals.o
	$(F90) -c $(F90FLAGS) shdPha.f90

shdExt_mpi_.o: shdExt_mpi_.f90
	$(F90) -c $(F90FLAGS) shdExt_mpi_.f90

shdExt_mpi.o: shdExt_mpi.f90
	$(F90) -c $(F90FLAGS) shdExt_mpi.f90

shdUtil.o: shdUtil.f90 globals.o
	$(F90) -c $(F90FLAGS) shdUtil.f90

globals.o: globals.f90
	$(F90) -c $(F90FLAGS) globals.f90

fftpack.o: fftpack.f90
	$(F90) -c $(F90FLAGS) fftpack.f90

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
