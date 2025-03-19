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

TARGETS = bin shdom _shdom propgen make_mie_table make_tmatrix_table plotscattab put

SHDOMOBJS = shdSolv_1.o shdSolv_2.o shdUtil_io.o shdSurf_brdf.o shdSurf_ocean.o _fftpack.o

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
	SHDOMOBJS += shdom_nonetcdf.o
else
	SHDOMOBJS += shdom_netcdf.o
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
_shdom: _shdom.o $(SHDOMOBJS)
	$(F90) $(F90FLAGS) _shdom.o $(SHDOMOBJS) $(LIBS) -o $(BIN_DIR)/_shdom

_shdom.o: _shdom.f90 $(SHDOMOBJS)
	$(F90) -c $(F90FLAGS) _shdom.f90

shdSolv_1.o: shdSolv_1.f90
	$(F90) -c $(F90FLAGS) shdSolv_1.f90

shdSolv_2.o: shdSolv_2.f90
	$(F90) -c $(F90FLAGS) shdSolv_2.f90

shdUtil_io.o: shdUtil_io.f90
	$(F90) -c $(F90FLAGS) shdUtil_io.f90

shdSurf_brdf.o: shdSurf_brdf.f90
	$(F90) -c $(F90FLAGS) shdSurf_brdf.f90

shdSurf_ocean.o: shdSurf_ocean.f90
	$(F90) -c $(F90FLAGS) shdSurf_ocean.f90

shdom_mpi.o: shdom_mpi.f90
	$(F90) -c $(F90FLAGS) shdom_mpi.f90

shdom_nompi.o: shdom_nompi.f90
	$(F90) -c $(F90FLAGS) shdom_nompi.f90

shdom_nonetcdf.o: shdom_nonetcdf.f90
	$(F90) -c $(F90FLAGS) shdom_nonetcdf.f90

shdom_netcdf.o: shdom_netcdf.f90
	$(F90) -c $(F90FLAGS) shdom_netcdf.f90

txt2nc: shdTool_txt2nc
	$(F90) $(F90FLAGS) shdTool_txt2nc.o shdom_netcdf.o $(LIBS) -o $(BIN_DIR)/txt2nc

shdTool_txt2nc.o: shdTool_txt2nc.f90 shdom_netcdf.o
	$(F90) -c $(F90FLAGS) shdTool_txt2nc.f90

propgen: shdTool_propgen.f90
	$(F90) $(F90FLAGS) shdTool_propgen.f90 -o $(BIN_DIR)/propgen

make_tmatrix_table: shdTool_make_tmx_tab.f90
	$(F90) $(F90FLAGS) shdTool_make_tmx_tab.f90 shdScat_ind_wat.f90 shdScat_tmx_wig.f90 _lapack.f -o $(BIN_DIR)/make_tmatrix_table

make_mie_table: shdTool_make_mie_tab.f90
	$(F90) $(F90FLAGS) shdTool_make_mie_tab.f90 shdScat_ind_wat.f90 shdScat_mie_wig.f90 -o $(BIN_DIR)/make_mie_table

plotscattab: shdTool_plot_sca_tab.f90
	$(F90) $(F90FLAGS) shdTool_plot_sca_tab.f90 -o $(BIN_DIR)/plotscattab

put: shdTool_put.f90
	$(F90) $(F90FLAGS) shdTool_put.f90 -o $(BIN_DIR)/put

.f.o :; $(FC) -c $(FFLAGS) $*.f

clean:
	@ rm -f *.o *.mod $(TARGETS) *~
