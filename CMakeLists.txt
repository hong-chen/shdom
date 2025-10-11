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
#   Created by:  Cheng Da
#                cda@umd.edu
#
#   Modified by: Hong Chen
#                hong.chen-1@colorado.edu
#╰────────────────────────────────────────────────────────────────────────────╯!

cmake_minimum_required(VERSION 3.10)

project(SHDOMNG LANGUAGES Fortran)

# Fortran Compiler and Flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-O3 -ffree-line-length-none")
    add_definitions(-cpp)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10.0)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wargument-mismatch")
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()
endif()

# openMP Support
option(USE_OMP "build shdom with openmp" OFF)
if(USE_OMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# MPI Support
option(USE_MPI "build shdom with mpi" OFF)
if(USE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

# Check the compatibility of the parallelization options
if(USE_OMP AND USE_MPI)
    message(FATAL_ERROR "USE_OMP and USE_MPI cannot be enabled at the same time. \nSet one of them OFF")
endif()

# check NetCDF support
option(USE_NETCDF "build shdom with netcdf" OFF)
if (USE_NETCDF)
    if (DEFINED NC_INC AND DEFINED NC_LIBS )
        message(STATUS "NC_INC= ${NC_INC}")
        message(STATUS "NC_LIBS= ${NC_LIBS}")
    else()
        message(FATAL_ERROR "USE_NETCDF=on, but NC_INC, NC_LIBS not set")
    endif()
endif()

add_subdirectory(src)
