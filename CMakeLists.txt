# This file is referenced from https://github.com/cd10kfsu/F90GIO, https://github.com/UMD-AOSC/Ocean-LETKF

cmake_minimum_required(VERSION 3.10)

project(SHDOM LANGUAGES Fortran)

option(BUILD_SHDOM_MPI "build shdom with mpi" OFF)
