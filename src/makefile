# Makefile for SHDOM and related programs in distribution
#
# Set the UseMPI and UseNetcdf variables to make the desired shdom version.

FC     = gfortran                        # your Fortran 90 compiler here (same as used by mpif90)
FFLAGS = -O3 -fallow-argument-mismatch   # your flags here
#FC     = pgf90
#FFLAGS = -fast
MPIFC  = mpif90    # MPI Fortran 90 compiling command for shdom90
MPIFLAGS = -O3     #   MPI compiling flags
LIBS   =

CURRENT_DIR = $(shell pwd)
BIN_DIR     = $(abspath $(CURRENT_DIR)/../bin)
BUILD_DIR   = $(abspath $(CURRENT_DIR)/../build)

SHDOMOBJS   = shdomsub1.o shdomsub2.o shdomsub3.o fftpack.o ocean_brdf.o

all:
	@ echo $(BUILD_DIR)
