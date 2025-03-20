# Makefile for SHDOM and related programs in distribution

# pick one of the following two
#╭────────────────────────────────────────────────╮#
FC = gfortran
#FC = mpif90
#╰────────────────────────────────────────────────╯#

FORTRAN_VERSION = $(shell $(FC) -dumpversion)
VERSION_LIMIT = 10
ifeq ($(shell test $(FORTRAN_VERSION) -lt $(VERSION_LIMIT); echo $$?), 0)
	FFLAGS = -O3 -Wargument-mismatch
else
	FFLAGS = -O3 -fallow-argument-mismatch
endif

CURRENT_DIR = $(shell pwd)
BIN_DIR     = $(abspath $(CURRENT_DIR)/../bin)

TARGETS = bin shdom shdom_ propgen make_mie_table make_tmatrix_table plotscattab put

SHDOMOBJS = shdSolv_1.o shdSolv_2.o shdUtil_io.o shdSurf_brdf.o fftpack.o globals.o

ifeq ($(FC), mpif90)
	UseMPI = yes
else
	UseMPI = no
endif

ifeq ($(UseMPI),no)
	SHDOMOBJS += shdom_nompi.o
else
	SHDOMOBJS += shdom_mpi.o
endif

F90 = $(FC)
F90FLAGS = $(FFLAGS)
LIBS =

UseNetcdf = no

NetcdfHome = /usr/local/netcdf-3.6.3
NetcdfInl  = -I$(NetcdfHome)/include
NetcdfLib  = -L$(NetcdfHome)/lib -lnetcdf

ifeq ($(UseNetcdf),no)
	SHDOMOBJS += shdom_nonc.o
else
	SHDOMOBJS += shdom_nc.o
	F90FLAGS += $(NetcdfInl)
	LIBS += $(NetcdfLib)
endif

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

shdSolv_1.o: shdSolv_1.f90
	$(F90) -c $(F90FLAGS) shdSolv_1.f90

shdSolv_2.o: shdSolv_2.f90
	$(F90) -c $(F90FLAGS) shdSolv_2.f90

shdUtil_io.o: shdUtil_io.f90 globals.o
	$(F90) -c $(F90FLAGS) shdUtil_io.f90

shdSurf_brdf.o: shdSurf_brdf.f90
	$(F90) -c $(F90FLAGS) shdSurf_brdf.f90

shdom_nompi.o: shdom_nompi.f90
	$(F90) -c $(F90FLAGS) shdom_nompi.f90

shdom_mpi.o: shdom_mpi.f90
	$(F90) -c $(F90FLAGS) shdom_mpi.f90

shdom_nonc.o: shdom_nonc.f90
	$(F90) -c $(F90FLAGS) shdom_nonc.f90

shdom_nc.o: shdom_nc.f90
	$(F90) -c $(F90FLAGS) shdom_nc.f90

fftpack.o: fftpack.f90
	$(F90) -c $(F90FLAGS) fftpack.f90

globals.o: globals.f90
	$(F90) -c $(F90FLAGS) globals.f90

txt2nc: shdTool_txt2nc
	$(F90) $(F90FLAGS) shdTool_txt2nc.o shdom_nc.o $(LIBS) -o $(BIN_DIR)/txt2nc

shdTool_txt2nc.o: tools/shdTool_txt2nc.f90 shdom_nc.o
	$(F90) -c $(F90FLAGS) tools/shdTool_txt2nc.f90

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
