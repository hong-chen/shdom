#
# src codes
#
set(MIE_LIB_SRCS shdScat_ind_wat.f90 shdScat_mie_wig.f90)

set(TMTX_LIB_SRCS shdScat_ind_wat.f90 shdScat_tmx_wig.f90 _lapack.f)

#
# libraries
#
add_library(mie ${MIE_LIB_SRCS})
target_link_libraries(mie PUBLIC ${SHDOM_LD_LIBS})

add_library(tmtx ${TMTX_LIB_SRCS})
target_link_libraries(tmtx PUBLIC ${SHDOM_LD_LIBS})

#
# executables
#
add_executable(make_tmatrix_table shdTool_make_tmx_tab.f90)
target_link_libraries(make_tmatrix_table PRIVATE tmtx ${SHDOM_LD_LIBS})

add_executable(make_mie_table shdTool_make_mie_tab.f90)
target_link_libraries(make_mie_table PRIVATE mie ${SHDOM_LD_LIBS})

add_executable(propgen shdTool_propgen.f90)

add_executable(plotscattab shdTool_plot_sca_tab.f90)

add_executable(put shdTool_put.f90)

set (exes make_tmatrix_table make_mie_table propgen plotscattab put)

if (USE_NETCDF)
    add_executable(txt2nc shdTool_txt2nc.f90 shdom_nc.f90)
    target_include_directories(txt2nc PUBLIC "$<BUILD_INTERFACE:${SHDOM_INC}>")
    target_link_libraries(txt2nc PRIVATE ${SHDOM_LD_LIBS})

    set(exes ${exes} txt2nc)
endif()

foreach(x ${exes})
    install(TARGETS ${x})
endforeach()
