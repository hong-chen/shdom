# Makefile for SHDOM and related programs in distribution
#
# Set the UseMPI and UseNetcdf variables to make the desired shdom version.

FC     = gfortran  # your Fortran 90 compiler here (same as used by mpif90)
FFLAGS = -O3       #   your flags here
#FC     = pgf90
#FFLAGS = -fast
MPIFC  = mpif90    # MPI Fortran 90 compiling command for shdom90
MPIFLAGS = -O3     #   MPI compiling flags
LIBS   =

SHDOMOBJS   = shdomsub1.o shdomsub2.o shdomsub3.o fftpack.o ocean_brdf.o

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

all: shdom

shdom: shdom.o $(SHDOMOBJS)
	$(F90) $(F90FLAGS) -o shdom shdom.o $(SHDOMOBJS) $(LIBS)

shdom.o: shdom_main.f90 $(SHDOMOBJS)
	$(F90) -c $(F90FLAGS) shdom_main.f90

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

.f.o :; $(FC) -c $(FFLAGS) $*.f
