cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
enable_language(C)
enable_language(Fortran)

# This is where the CCSE lib was installed
set(CCSE_DIR $ENV{BOXLIB_HOME})

# Set some defaults
set(BL_SPACEDIM 3 CACHE INT "Dimension of BoxLib build")
set(ENABLE_MPI 0 CACHE INT "Enable build with MPI")
set(ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
set(BL_PRECISION "DOUBLE" CACHE INT "Precision of BoxLib build")
set(BL_USE_PARTICLES 0 CACHE INT "Include Particles classes in BoxLib build")
set(ENABLE_PROFILING 0 CACHE INT "Include profiling information in BoxLib build")
set(ENABLE_BACKTRACE 1 CACHE INT "Include backtrace information in BoxLib build")


set(CMAKE_MODULE_PATH ${CCSE_DIR}/cmake)
find_package(CCSE REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIR})
include_directories(model)
include_directories(model/include)
link_directories(${CCSE_LIBRARY_DIR})

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)

project(MFIXEXA)

set(CXX_sources
  model/main.cpp
  model/mfix_level.cpp  )
set(CXX_includes
  model/mfix_level.H
  model/include/mfix_F.H)
set(F90_sources
  model/adjust_a_g.f90
  model/adjust_dt_mod.f90
  model/adjust_leq.f90
  model/bc_mod.f90
  model/calc_cell.f90
  model/calc_coeff.f90
  model/calc_d.f90
  model/calc_mflux.f90
  model/calc_mu_g.f90
  model/calc_outflow.f90
  model/calc_trd_g.f90
  model/check_bqend.f90
  model/check_convergence.f90
  model/check_data/check_axis.f90
  model/check_data/check_bc_dem.f90
  model/check_data/check_bc_geometry.f90
  model/check_data/check_bc_inflow.f90
  model/check_data/check_bc_outflow.f90
  model/check_data/check_bc_walls.f90
  model/check_data/check_boundary_conditions.f90
  model/check_data/check_dmp_prereqs.f90
  model/check_data/check_gas_phase.f90
  model/check_data/check_geometry.f90
  model/check_data/check_geometry_prereqs.f90
  model/check_data/check_ic_common_discrete.f90
  model/check_data/check_initial_conditions.f90
  model/check_data/check_numerics.f90
  model/check_data/check_output_control.f90
  model/check_data/check_point_sources.f90
  model/check_data/check_run_control.f90
  model/check_data/check_solids_common_all.f90
  model/check_data/check_solids_common_discrete.f90
  model/check_data/check_solids_dem.f90
  model/check_data/check_solids_model_prereqs.f90
  model/check_data/check_solids_phases.f90
  model/check_data_20.f90
  model/check_data_30.f90
  model/check_plane.f90
  model/constant_mod.f90
  model/conv_dif_u_g.f90
  model/conv_dif_v_g.f90
  model/conv_dif_w_g.f90
  model/conv_pp_g.f90
  model/conv_rop.f90
  model/corner.f90
  model/correct_0.f90
  model/deprecated.f90
  model/des/calc_collision_wall_mod.f90
  model/des/calc_drag_des.f90
  model/des/calc_force_dem.f90
  model/des/calc_grad_des.f90
  model/des/calc_pg_grad.f90
  model/des/cffctowall.f90
  model/des/cfnewvalues.f90
  model/des/cfrelvel.f90
  model/des/cfslide.f90
  model/des/comp_mean_fields.f90
  model/des/des_bc_mod.f90
  model/des/des_init_namelist.f90
  model/des/des_time_march.f90
  model/des/desmpi_mod.f90
  model/des/discretelement_mod.f90
  model/des/drag_gp_des.f90
  model/des/drag_gs_des.f90
  model/des/gas_drag.f90
  model/des/layout_mi_dem.f90
  model/des/make_arrays_des.f90
  model/des/randomno_mod.f90
  model/des/read_particle_input.f90
  model/des/read_res0_des.f90
  model/des/read_res1_des_mod.f90
  model/des/set_geometry_des.f90
  model/des/set_phase_index.f90
  model/des/stl_functions_des_mod.f90
  model/des/vtp_mod.f90
  model/des/write_des_data.f90
  model/des/write_res0_des.f90
  model/des/write_res1_des_mod.f90
  model/discretization_mod.f90
  model/display_resid.f90
  model/dmp_modules/compar_mod.f90
  model/dmp_modules/gridmap_mod.f90
  model/drag_mod.f90
  model/eos_mod.f90
  model/error_manager_mod.f90
  model/exit.f90
  model/fld_constants_mod.f90
  model/flow_to_vel.f90
  model/functions_mod.f90
  model/geometry_mod.f90
  model/get_bc_area.f90
  model/get_data.f90
  model/get_ps.f90
  model/get_tunit.f90
  model/goal_seek_mflux.f90
  model/ic_mod.f90
  model/init_namelist.f90
  model/io/funits_mod.f90
  model/io/in_binary_512_mod.f90
  model/io/in_binary_512i_mod.f90
  model/io/open_file.f90
  model/io/open_files.f90
  model/io/out_array.f90
  model/io/out_array_k.f90
  model/io/out_array_kc.f90
  model/io/out_bin_512.f90
  model/io/out_bin_512i.f90
  model/io/out_bin_512r.f90
  model/io/out_bin_r.f90
  model/io/output_manager.f90
  model/io/read_res1.f90
  model/io/write_ab_m.f90
  model/io/write_ab_m_var.f90
  model/io/write_error.f90
  model/io/write_header.f90
  model/io/write_out0.f90
  model/io/write_out1.f90
  model/io/write_out3.f90
  model/io/write_res0.f90
  model/io/write_res1.f90
  model/io/write_table.f90
  model/location_check.f90
  model/machine_mod.f90
  model/matrix_mod.f90
  model/mfix.f90
  model/mod_bc.f90
  model/output_mod.f90
  model/param1_mod.f90
  model/param_mod.f90
  model/parse_line.f90
  model/parse_resid_string.f90
  model/physical_prop.f90
  model/ps_mod.f90
  model/read_namelist.f90
  model/remove_comment.f90
  model/residual_mod.f90
  model/run_mod.f90
  model/scales_mod.f90
  model/set_bc0.f90
  model/set_bc1.f90
  model/set_bc_flow.f90
  model/set_constprop.f90
  model/set_domain.f90
  model/set_flags.f90
  model/set_fluidbed_p.f90
  model/set_ic.f90
  model/set_icbc_flags.f90
  model/set_max2.f90
  model/set_outflow.f90
  model/set_param.f90
  model/set_ps.f90
  model/set_ro_g.f90
  model/set_wall_bc.f90
  model/solve_pp_g.f90
  model/solve_vel_star.f90
  model/solvers/BLAS.f90
  model/solvers/DGTSV.f90
  model/solvers/dgtsl.f90
  model/solvers/leq_bicgs.f90
  model/solvers/leqsol_mod.f90
  model/solvers/solve_lin_eq.f90
  model/solvers/xerbla.f90
  model/source_pp_g.f90
  model/source_u_g.f90
  model/source_v_g.f90
  model/source_w_g.f90
  model/tau_u_g.f90
  model/tau_v_g.f90
  model/tau_w_g.f90
  model/time_cpu_mod.f90
  model/toleranc_mod.f90
  model/ur_facs_mod.f90
  model/usr/usr0.f90
  model/usr/usr0_des.f90
  model/usr/usr1.f90
  model/usr/usr1_des.f90
  model/usr/usr2.f90
  model/usr/usr2_des.f90
  model/usr/usr3.f90
  model/usr/usr3_des.f90
  model/usr/usr_drag.f90
  model/usr/usr_init_namelist.f90
  model/usr/usr_mod.f90
  model/usr/usr_physical_prop.f90
  model/usr/usr_write_out0.f90
  model/usr/usr_write_out1.f90
  model/usr/write_usr0.f90
  model/usr/write_usr1.f90
  model/utilities_mod.f90
  model/vavg_g.f90
  model/xsi_mod.f90
  model/zero_norm_vel.f90
  model/boxlib_to_mfix.F90)
set(F77_sources  )
set(FPP_includes )
set(FPP_sources  )

include(PreprocessBoxLibFortran)
preprocess_boxlib_fortran(FPP_out ${FPP_sources})

FILE(GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
set(local_includes ${CXX_includes} ${FPP_includes})
set(local_sources  ${CXX_sources} ${F90_sources} ${F77_sources} ${FPP_out} )

add_executable(mfix model/main.cpp ${USR_OVERRIDES})
add_library(model STATIC ${local_sources} ${local_includes})

target_link_libraries(mfix model ${CCSE_LIBRARIES} ${MPI_LIBRARIES})


enable_testing()

set(TESTS
  tests/FLD01
  tests/FLD02
  tests/DEM01
  tests/DEM02
  tests/DEM03
  tests/DEM04
  tests/DEM05
  tests/DEM06
  #tests/DEM07
  #benchmarks/DEM01
  #benchmarks/DEM02
  #benchmarks/DEM03
  #benchmarks/DEM04
  #benchmarks/FLD01
  )

set (BUILD_TEST build_test)
add_custom_target(${BUILD_TEST})

foreach( rel_testdir ${TESTS})

  STRING(REPLACE "/" "_" TEST_NAME ${rel_testdir} )

  set(testdir ${CMAKE_CURRENT_SOURCE_DIR}/${rel_testdir})

  FILE(GLOB USR_OVERRIDES ${testdir}/*.f90)

  set(TEST_EXE mfix_${TEST_NAME})

  add_executable(${TEST_EXE} EXCLUDE_FROM_ALL model/main.cpp ${USR_OVERRIDES})

  target_include_directories(${TEST_EXE} PUBLIC ${testdir} ${CMAKE_Fortran_MODULE_DIRECTORY})
  set_target_properties(${TEST_EXE} PROPERTIES Fortran_MODULE_DIRECTORY ${testdir})

  add_test(NAME ${TEST_NAME}
    COMMAND ${CMAKE_COMMAND}
            -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
            -DCMAKE_COMMAND=${CMAKE_COMMAND}
            -DTEST_EXE=${TEST_EXE}
            -DTEST_NAME=${TEST_NAME}
            -Dtestdir=${testdir}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtests.cmake)


  add_dependencies ( ${BUILD_TEST} ${TEST_EXE})

  set_tests_properties ( ${TEST_NAME}
    PROPERTIES ENVIRONMENT "MFIX=${CMAKE_BINARY_DIR}/${TEST_EXE}"
    )

  target_link_libraries(${TEST_EXE} ${CCSE_LIBRARIES} ${MPI_LIBRARIES} model)

endforeach( rel_testdir ${TESTS})
