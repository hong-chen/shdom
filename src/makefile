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
SHDOMOBJS = globals.o shdUtil.o shdSfc.o shdPha.o shdInit.o\
			shdGrid.o shdSrc.o shdSHDO.o shdIntg.o shdSolv.o\
			shdComp.o shdInp.o shdOut.o fftpack.o

ifeq ($(FC), mpif90)
	UseMPI = yes
else
	UseMPI = no
endif

ifeq ($(UseMPI),no)
	SHDOMMPIOBJ = shdMPI_.o
else
	SHDOMMPIOBJ = shdMPI.o
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

shdSolv.o: shdSolv.f90 $(SHDOMMPIOBJ) shdSfc.o shdInit.o shdGrid.o shdSHDO.o shdPha.o shdSrc.o shdIntg.o shdComp.o
	$(F90) -c $(F90FLAGS) shdSolv.f90

shdIntg.o: shdIntg.f90 $(SHDOMMPIOBJ) shdSfc.o shdInit.o shdGrid.o shdSHDO.o
	$(F90) -c $(F90FLAGS) shdIntg.f90

shdComp.o: shdComp.f90 globals.o shdPha.o shdSrc.o shdInit.o shdGrid.o
	$(F90) -c $(F90FLAGS) shdComp.f90

shdInp.o: shdInp.f90 globals.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdInp.f90

shdOut.o: shdOut.f90 globals.o shdGrid.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdOut.f90

shdMPI_.o: shdMPI_.f90
	$(F90) -c $(F90FLAGS) shdMPI_.f90

shdMPI.o: shdMPI.f90
	$(F90) -c $(F90FLAGS) shdMPI.f90

shdSrc.o: shdSrc.f90 globals.o shdUtil.o shdPha.o shdGrid.o shdInit.o
	$(F90) -c $(F90FLAGS) shdSrc.f90

shdInit.o: shdInit.f90 globals.o shdUtil.o shdSfc.o
	$(F90) -c $(F90FLAGS) shdInit.f90

shdSHDO.o: shdSHDO.f90 globals.o shdPha.o shdSfc.o
	$(F90) -c $(F90FLAGS) shdSHDO.f90

shdGrid.o: shdGrid.f90 globals.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdGrid.f90

shdSfc.o: shdSfc.f90 globals.o
	$(F90) -c $(F90FLAGS) shdSfc.f90

shdPha.o: shdPha.f90 globals.o
	$(F90) -c $(F90FLAGS) shdPha.f90

shdUtil.o: shdUtil.f90 globals.o
	$(F90) -c $(F90FLAGS) shdUtil.f90

globals.o: globals.f90
	$(F90) -c $(F90FLAGS) globals.f90

fftpack.o: fftpack.f90
	$(F90) -c $(F90FLAGS) fftpack.f90

# =============================================================

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
