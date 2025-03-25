#
# src codes
#
set(SHDOM_LIB_SRCS
    shdSolv_1.f90
    shdSolv_2.f90
    shdUtil_io.f90
    shdSurf_brdf.f90
    fftpack.f90)
if (USE_MPI)
    set(SHDOM_LIB_SRCS ${SHDOM_LIB_SRCS} shdom_mpi.f90)
else()
    set(SHDOM_LIB_SRCS ${SHDOM_LIB_SRCS} shdom_nompi.f90)
endif()

if (USE_NETCDF)
    set(SHDOM_LIB_SRCS ${SHDOM_LIB_SRCS} shdom_nc.f90)
else()
    set(SHDOM_LIB_SRCS ${SHDOM_LIB_SRCS} shdom_nonc.f90)
endif()

#
# set compiling & linking flags
#
set(SHDOM_LD_LIBS "")
set(SHDOM_INC "")
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
add_library(shdom_f90 ${SHDOM_LIB_SRCS})
target_include_directories(shdom_f90 PUBLIC "$<BUILD_INTERFACE:${SHDOM_INC}>")
target_link_libraries(shdom_f90 PUBLIC ${SHDOM_LD_LIBS})

#
# build executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(exes shdom shdom_)
foreach (x ${exes})
    add_executable(${x} ${x}.f90)
    target_include_directories(${x} PUBLIC "$<BUILD_INTERFACE:${SHDOM_INC}>")
    target_link_libraries(${x} PRIVATE shdom_f90 ${SHDOM_LD_LIBS})
    install(TARGETS ${x})
endforeach()


add_subdirectory(tools)

