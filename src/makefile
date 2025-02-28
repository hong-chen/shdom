# by Hong Chen
# Notes for Mac:
# export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"
# FFLAGS = -O3 -fallow-argument-mismatch
# PATH="/Users/hchen/Work/mygit/shdom/bin:$PATH"
#
#
# Makefile for SHDOM and related programs in distribution
#
# Set the UseMPI and UseNetcdf variables to make the desired shdom version.

FC     = gfortran                        # your Fortran 90 compiler here (same as used by mpif90)
CC     = gcc
CFLAGS =
FFLAGS = -O3 -fallow-argument-mismatch   # your flags here
#FC     = pgf90
#FFLAGS = -fast
MPIFC  = mpif90    # MPI Fortran 90 compiling command for shdom90
MPIFLAGS = -O3     #   MPI compiling flags
LIBS   =

CURRENT_DIR = $(shell pwd)
BIN_DIR     = $(abspath $(CURRENT_DIR)/../bin)
BUILD_DIR   = $(abspath $(CURRENT_DIR)/../build)

TARGETS = shdom propgen make_mie_table make_tmatrix_table plotscattab put

SHDOMOBJS = shdomsub1.o shdomsub2.o shdomsub3.o fftpack.o shdSfc_ocean.o

UseMPI = no
UseNetcdf = no

NetcdfHome = /usr/local/netcdf-3.6.3
NetcdfInl  = -I$(NetcdfHome)/include
NetcdfLib  = -L$(NetcdfHome)/lib -lnetcdf


ifeq ($(UseMPI),no)
  SHDOMOBJS += shdom_nompi.o
  F90 = $(FC)
  F90FLAGS = $(FFLAGS)
else
  SHDOMOBJS += shdom_mpi.o
  F90 = $(MPIFC)
  F90FLAGS = $(MPIFLAGS)
endif

ifeq ($(UseNetcdf),no)
  SHDOMOBJS += shdom_nonetcdf.o
else
  SHDOMOBJS += shdom_netcdf.o
  F90FLAGS += $(NetcdfInl)
  LIBS   += $(NetcdfLib)
endif

all: $(TARGETS)

shdom: shdom.o $(SHDOMOBJS)
	$(F90) $(F90FLAGS) shdom.o $(SHDOMOBJS) $(LIBS) -o $(BIN_DIR)/shdom

shdom.o: shdom.f90 $(SHDOMOBJS)
	$(F90) -c $(F90FLAGS) shdom.f90

shdomsub1.o: shdomsub1.f
	$(F90) -c $(F90FLAGS) shdomsub1.f

shdomsub2.o: shdomsub2.f
	$(F90) -c $(F90FLAGS) shdomsub2.f

shdomsub3.o: shdomsub3.f
	$(F90) -c $(F90FLAGS) shdomsub3.f

shdom_mpi.o: shdom_mpi.f90
	$(F90) -c $(F90FLAGS) shdom_mpi.f90

shdom_nompi.o: shdom_nompi.f
	$(F90) -c $(F90FLAGS) shdom_nompi.f

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

make_tmatrix_table: shdTool_make_tmatrix_table.f90
	$(F90) $(F90FLAGS) shdTool_make_tmatrix_table.f90 shdSca_indexwatice.f shdSca_tmatrixwig.f lpd.f -o $(BIN_DIR)/make_tmatrix_table

make_mie_table: shdTool_make_mie_table.f90
	$(F90) $(F90FLAGS) shdTool_make_mie_table.f90 shdSca_indexwatice.f shdSca_miewig.f -o $(BIN_DIR)/make_mie_table

plotscattab: shdTool_plotscattab.f90
	$(F90) $(F90FLAGS) shdTool_plotscattab.f90 -o $(BIN_DIR)/plotscattab

put: shdTool_put.c
	$(CC) $(CFLAGS) shdTool_put.c -o $(BIN_DIR)/put

.f.o :; $(FC) -c $(FFLAGS) $*.f

clean:
	@ rm -f *.o *.mod $(TARGETS) *~
