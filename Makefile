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
#   This file was referenced from https://github.com/cd10kfsu/F90GIO
#╰────────────────────────────────────────────────────────────────────────────╯!

CURRENT_DIR = $(shell pwd)
BUILD_DIR = $(abspath $(CURRENT_DIR)/build)

all:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR); \
		cmake $(CURRENT_DIR); \
		make

#omp:
#    mkdir -p $(BUILD_DIR)
#    cd $(BUILD_DIR); \
#        cmake $(CURRENT_DIR) \
#            -DUSE_OMP=on; \
#        make

mpi:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR); \
		cmake $(CURRENT_DIR) \
			-DUSE_MPI=on; \
		make

#nc:
#    mkdir -p $(BUILD_DIR)
#    cd $(BUILD_DIR); \
#        cmake $(CURRENT_DIR) \
#            -DUSE_NETCDF=on \
#            -DNC_INC="/opt/homebrew/Cellar/netcdf-fortran/4.6.1_1/include" \
#            -DNC_LIBS="-L/opt/homebrew/Cellar/netcdf-fortran/4.6.1_1/lib -lnetcdff"; \
#        make VERBOSE=1

clean:
	rm -rf $(BUILD_DIR)
