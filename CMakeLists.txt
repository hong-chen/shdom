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
