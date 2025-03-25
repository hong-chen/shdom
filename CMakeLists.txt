# This file is referenced from https://github.com/cd10kfsu/F90GIO, https://github.com/UMD-AOSC/Ocean-LETKF

cmake_minimum_required(VERSION 3.10)

project(SHDOM LANGUAGES Fortran)

# Fortran Compiler and Flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-O3 -ffree-line-length-none")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10.0)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wargument-mismatch")
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()
endif()

# MPI Support
option(USE_MPI "build shdom with mpi" OFF)
if(USE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

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
