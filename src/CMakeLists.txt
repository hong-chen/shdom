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
#   Cheng Da, cda@umd.edu
#   Hong Chen, hong.chen-1@colorado.edu
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
#╰────────────────────────────────────────────────────────────────────────────╯!

#
# src codes
#
set(SHDOM_LIB_SRCS
    globals.f90
    shdUtil.f90
    shdSfc.f90
    shdPha.f90
    shdGrid.f90
    shdInit.f90
    shdSrc.f90
    shdSHDO.f90
    shdIntg.f90
    shdSolv.f90
    shdRad.f90
    shdInp.f90
    shdOut.f90
    shdMPI.f90
    fftpack.f90)

# dummy option
#if (USE_OMP)
#    set(SHDOM_LIB_SRCS ${SHDOM_LIB_SRCS})
#else()
#    set(SHDOM_LIB_SRCS ${SHDOM_LIB_SRCS})
#endif()

if (USE_MPI)
    set(HAS_MPI TRUE)
    set(FPPS_FLAGS ${FPPS_FLAGS} HAS_MPI)
else()
    set(HAS_MPI FALSE)
endif()

#if (USE_NETCDF)
#    set(SHDOM_LIB_SRCS ${SHDOM_LIB_SRCS} tools/shdExt_nc.f90)
#else()
#    set(SHDOM_LIB_SRCS ${SHDOM_LIB_SRCS} tools/shdExt_nc_.f90)
#endif()

#
# set compiling & linking flags
#
set(SHDOM_LD_LIBS "")
set(SHDOM_INC "")
if (USE_OMP)
    set(SHDOM_LD_LIBS ${SHDOM_LD_LIBS} OpenMP::OpenMP_Fortran)
endif()
if (USE_MPI)
    set(SHDOM_LD_LIBS ${SHDOM_LD_LIBS} MPI::MPI_Fortran)
endif()
if (USE_NETCDF)
    set(SHDOM_LD_LIBS ${SHDOM_LD_LIBS} ${NC_LIBS})
    set(SHDOM_INC ${SHDOM_INC} ${NC_INC})
endif()

#
# build libraries
#
add_library(shdomNG ${SHDOM_LIB_SRCS})
target_compile_definitions(shdomNG PUBLIC ${FPPS_FLAGS})
target_include_directories(shdomNG PUBLIC "$<BUILD_INTERFACE:${SHDOM_INC}>")
target_link_libraries(shdomNG PUBLIC ${SHDOM_LD_LIBS})

#
# build executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(exes shdom)
foreach (x ${exes})
    add_executable(${x} ${x}.f90)
    target_include_directories(${x} PUBLIC "$<BUILD_INTERFACE:${SHDOM_INC}>")
    target_link_libraries(${x} PRIVATE shdomNG ${SHDOM_LD_LIBS})
    install(TARGETS ${x})
endforeach()


add_subdirectory(tools)
