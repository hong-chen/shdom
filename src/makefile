#╭────────────────────────────────────────────────────────────────────────────╮!
# @LICENSE_HEADER_START@
#
#   This file is part of SHDOMNG.
#
#   --
#   SHDOMNG: Spherical Harmonic Discrete Ordinate Method - Next Generation
#   --
#
#   Copyright (C) 2025-
#   EaR³T Community: https://discord.gg/ntqsguwaWv
#
#   SHDOMNG is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   SHDOMNG is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with SHDOMNG. If not, see <http://www.gnu.org/licenses/>.
#
# @LICENSE_HEADER_END@
#
# Acknowledgement:
#   SHDOM was originally developed by K. Franklin Evans and is available
#   at https://coloradolinux.com/shdom.
#
# History:
#   Created by:  K. Franklin Evans
#                evans@nit.colorado.edu
#
#   Modified by: Hong Chen
#                hong.chen-1@colorado.edu
#╰────────────────────────────────────────────────────────────────────────────╯!

# Makefile for SHDOM and related programs in distribution

# pick one of the following two
#╭────────────────────────────────────────────────╮#
FC = gfortran
#FC = mpifort
#╰────────────────────────────────────────────────╯#

FORTRAN_VERSION = $(shell $(FC) -dumpversion)
VERSION_LIMIT = 10
ifeq ($(shell test $(FORTRAN_VERSION) -lt $(VERSION_LIMIT); echo $$?), 0)
	FFLAGS = -O3 -Ofast -Wargument-mismatch -ffree-line-length-none
else
	FFLAGS = -O3 -Ofast -fallow-argument-mismatch -ffree-line-length-none
endif

# dummy option, placeholder for OMP implementation in the future
useOMP = no
ifeq ($(useOMP), yes)
	# openMP
	FFLAGS += -fopenmp
endif

CURRENT_DIR = $(shell pwd)
BIN_DIR     = $(abspath $(CURRENT_DIR)/../build/bin)
TARGETS = bin shdom propgen make_mie_table make_tmatrix_table plotscattab put

SHDOMOBJS = globals.o shdUtil.o shdSfc.o shdScat.o shdInit.o\
			shdGrid.o shdSrc.o shdShdo.o shdIntg.o shdSolv.o\
			shdRad.o shdInp.o shdOut.o shdMpi.o fftpack.o

ifeq ($(FC), mpif90)
	FFLAGS += -cpp -DHAS_MPI=TRUE
else
	FFLAGS += -cpp -DHAS_MPI=FALSE
endif

F90 = $(FC)
F90FLAGS = $(FFLAGS)
LIBS =

all: $(TARGETS)

bin:
	mkdir -p $(BIN_DIR)

# this is for namelist input (exe: shdom)
shdom: shdom.o $(SHDOMOBJS)
	$(F90) $(F90FLAGS) shdom.o $(SHDOMOBJS) $(LIBS) -o $(BIN_DIR)/shdom

shdom.o: shdom.f90 $(SHDOMOBJS)
	$(F90) -c $(F90FLAGS) shdom.f90

shdMpi.o: shdMpi.f90 shdInit.o shdUtil.o shdSrc.o shdRad.o globals.o
	$(F90) -c $(F90FLAGS) shdMpi.f90

shdOut.o: shdOut.f90 globals.o shdGrid.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdOut.f90

shdInp.o: shdInp.f90 globals.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdInp.f90

shdSolv.o: shdSolv.f90 shdMpi.o shdSfc.o shdInit.o shdGrid.o shdShdo.o shdScat.o shdSrc.o shdIntg.o shdRad.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdSolv.f90

shdIntg.o: shdIntg.f90 shdMpi.o shdSfc.o shdInit.o shdGrid.o shdUtil.o shdShdo.o
	$(F90) -c $(F90FLAGS) shdIntg.f90

shdGrid.o: shdGrid.f90 shdMpi.o globals.o shdUtil.o shdRad.o shdSrc.o
	$(F90) -c $(F90FLAGS) shdGrid.f90

shdSrc.o: shdSrc.f90 globals.o shdUtil.o shdScat.o shdShdo.o shdInit.o
	$(F90) -c $(F90FLAGS) shdSrc.f90

shdShdo.o: shdShdo.f90 globals.o shdScat.o shdSfc.o shdUtil.o fftpack.o
	$(F90) -c $(F90FLAGS) shdShdo.f90

shdRad.o: shdRad.f90 globals.o shdScat.o shdSrc.o shdInit.o
	$(F90) -c $(F90FLAGS) shdRad.f90

shdInit.o: shdInit.f90 globals.o shdUtil.o shdSfc.o
	$(F90) -c $(F90FLAGS) shdInit.f90

shdScat.o: shdScat.f90 globals.o
	$(F90) -c $(F90FLAGS) shdScat.f90

shdSfc.o: shdSfc.f90 globals.o shdUtil.o
	$(F90) -c $(F90FLAGS) shdSfc.f90

shdUtil.o: shdUtil.f90 globals.o
	$(F90) -c $(F90FLAGS) shdUtil.f90

globals.o: globals.f90
	$(F90) -c $(F90FLAGS) globals.f90

fftpack.o: fftpack.f90
	$(F90) -c $(F90FLAGS) fftpack.f90

# this block is for tools that pre-process/post-process radiative properties
# for SHDOM - plan to replace them by EaR³T in the future
# =============================================================
propgen: tools/shdTool_propgen.f90
	$(F90) $(F90FLAGS) tools/shdTool_propgen.f90 -o $(BIN_DIR)/propgen

make_tmatrix_table: tools/shdTool_make_tmx_tab.f90
	$(F90) $(F90FLAGS) tools/shdTool_make_tmx_tab.f90 tools/shdPha_ind_wat.f90 tools/shdPha_tmx_wig.f90 tools/_lapack.f -o $(BIN_DIR)/make_tmatrix_table

make_mie_table: tools/shdTool_make_mie_tab.f90
	$(F90) $(F90FLAGS) tools/shdTool_make_mie_tab.f90 tools/shdPha_ind_wat.f90 tools/shdPha_mie_wig.f90 -o $(BIN_DIR)/make_mie_table

plotscattab: tools/shdTool_plot_sca_tab.f90
	$(F90) $(F90FLAGS) tools/shdTool_plot_sca_tab.f90 -o $(BIN_DIR)/plotscattab

put: tools/shdTool_put.f90
	$(F90) $(F90FLAGS) tools/shdTool_put.f90 -o $(BIN_DIR)/put
# =============================================================

clean:
	@ rm -f *.o *.mod $(TARGETS) *~
