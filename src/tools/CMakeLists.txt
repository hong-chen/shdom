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
#   CMake related contribution was from Cheng Da (cda@umd.edu).
#╰────────────────────────────────────────────────────────────────────────────╯!

#
# src codes
#
set(MIE_LIB_SRCS shdPha_ind_wat.f90 shdPha_mie_wig.f90)

set(TMX_LIB_SRCS shdPha_ind_wat.f90 shdPha_tmx_wig.f90 _lapack.f)

#
# libraries
#
add_library(mie ${MIE_LIB_SRCS})
target_link_libraries(mie PUBLIC ${SHDOMNG_LD_LIBS})

add_library(tmx ${TMX_LIB_SRCS})
target_link_libraries(tmx PUBLIC ${SHDOMNG_LD_LIBS})

#
# executables
#
add_executable(make_tmatrix_table shdTool_make_tmx_tab.f90)
target_link_libraries(make_tmatrix_table PRIVATE tmx ${SHDOMNG_LD_LIBS})

add_executable(make_mie_table shdTool_make_mie_tab.f90)
target_link_libraries(make_mie_table PRIVATE mie ${SHDOMNG_LD_LIBS})

add_executable(propgen shdTool_propgen.f90)

add_executable(plotscattab shdTool_plot_sca_tab.f90)

add_executable(put shdTool_put.f90)

set (exes make_tmatrix_table make_mie_table propgen plotscattab put)

if (USE_NETCDF)
    add_executable(txt2nc shdTool_txt2nc.f90 shdExt_nc.f90)
    target_include_directories(txt2nc PUBLIC "$<BUILD_INTERFACE:${SHDOMNG_INC}>")
    target_link_libraries(txt2nc PRIVATE ${SHDOMNG_LD_LIBS})

    set(exes ${exes} txt2nc)
endif()

foreach(x ${exes})
    install(TARGETS ${x})
endforeach()
