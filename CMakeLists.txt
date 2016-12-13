cmake_minimum_required(VERSION 2.8)
project("MFIX" Fortran)

include_directories(model/include)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)

option(MPI "MPI support" OFF)

if ( SMP )
  find_package ( OpenMP REQUIRED)
    message("openmp found")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ( SMP )

if ( MPI )
  find_package ( MPI REQUIRED)
  if ( MPI_Fortran_COMPILER )
    include_directories( ${MPI_INCLUDE_PATH} )
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    add_definitions(-DMPI)
  else ( MPI_Fortran_COMPILER )
    message (FATAL_ERROR " MPI support requested but MPI compiler not found ")
  endif ( MPI_Fortran_COMPILER )
endif ( MPI )

if (CMAKE_Fortran_COMPILER MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS "-fimplicit-none -fconvert=big-endian ${CMAKE_Fortran_FLAGS}")
  set (CMAKE_Fortran_FLAGS "-ffree-form ${CMAKE_Fortran_FLAGS}")
elseif (CMAKE_Fortran_COMPILER MATCHES "ifort.*")
  # ifort
  set (CMAKE_Fortran_FLAGS "-implicitnone -convert big_endian ${CMAKE_Fortran_FLAGS}")
  set (CMAKE_Fortran_FLAGS "-free ${CMAKE_Fortran_FLAGS}")
endif (CMAKE_Fortran_COMPILER MATCHES "gfortran.*")

set(SOURCES
  model/accum_resid.f90
  model/adjust_a_g.f90
  model/adjust_dt_mod.f90
  model/adjust_leq.f90
  model/allocate_arrays.f90
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
  model/des/calc_epg_des.f90
  model/des/calc_force_dem.f90
  model/des/calc_grad_des.f90
  model/des/calc_pg_grad.f90
  model/des/cfassign.f90
  model/des/cffctowall.f90
  model/des/cfnewvalues.f90
  model/des/cfrelvel.f90
  model/des/cfslide.f90
  model/des/check_cell_movement.f90
  model/des/comp_mean_fields.f90
  model/des/comp_mean_fields0.f90
  model/des/comp_mean_fields1.f90
  model/des/des_allocate_mod.f90
  model/des/des_bc_mod.f90
  model/des/des_init_namelist.f90
  model/des/des_time_march.f90
  model/des/desgrid_mod.f90
  model/des/desmpi_mod.f90
  model/des/discretelement_mod.f90
  model/des/drag_gp_des.f90
  model/des/drag_gs_des0.f90
  model/des/drag_gs_des1.f90
  model/des/gas_drag.f90
  model/des/generate_particles_mod.f90
  model/des/interpolation_mod.f90
  model/des/layout_mi_dem.f90
  model/des/make_arrays_des.f90
  model/des/mass_inflow_dem.f90
  model/des/mass_outflow_dem.f90
  model/des/mpi_comm_des_mod.f90
  model/des/mpi_funs_des_mod.f90
  model/des/mpi_init_des_mod.f90
  model/des/mpi_node_des_mod.f90
  model/des/mpi_pack_des_mod.f90
  model/des/mpi_unpack_des_mod.f90
  model/des/neighbour.f90
  model/des/particle_filter_mod.f90
  model/des/particles_in_cell.f90
  model/des/randomno_mod.f90
  model/des/read_particle_input.f90
  model/des/read_res0_des.f90
  model/des/read_res1_des_mod.f90
  model/des/set_bc_dem.f90
  model/des/set_bc_dem_mi.f90
  model/des/set_bc_dem_mo.f90
  model/des/set_filter_des.f90
  model/des/set_geometry_des.f90
  model/des/set_phase_index.f90
  model/des/stl_dbg_des_mod.f90
  model/des/stl_functions_des_mod.f90
  model/des/stl_mod.f90
  model/des/stl_preproc_des_mod.f90
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
  model/iterate.f90
  model/location_check.f90
  model/machine_mod.f90
  model/matrix_mod.f90
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
  model/set_flags.f90
  model/set_fluidbed_p.f90
  model/set_geometry.f90
  model/set_geometry1.f90
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
  model/time_march.f90
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
  model/zero_norm_vel.f90)

FILE(GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
add_executable(mfix model/mfix.f90 ${USR_OVERRIDES})
add_library(model STATIC ${SOURCES})

set_source_files_properties(model/mfix.f90 ${SOURCES} ${USR_OVERRIDES} PROPERTIES COMPILE_FLAGS -cpp)

set_target_properties(mfix PROPERTIES Fortran_FORMAT "FREE")
set_target_properties(model PROPERTIES Fortran_FORMAT "FREE")

target_link_libraries(mfix model)

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
  tests/DEM07
  benchmarks/DEM01
  benchmarks/DEM02
  benchmarks/DEM03
  benchmarks/FLD01
  # box3d
  # cyclone
  # fluidbed_DES
  )

foreach( testdir ${TESTS})

  STRING(REPLACE "/" "_" testname ${testdir} )

  FILE(GLOB USR_OVERRIDES ${testdir}/*.f90)
  add_executable(mfix_${testname} EXCLUDE_FROM_ALL model/mfix.f90 ${USR_OVERRIDES})

  set_target_properties(mfix_${testname} PROPERTIES Fortran_FORMAT "FREE")
  set_source_files_properties(model/mfix.f90 ${USR_OVERRIDES} PROPERTIES COMPILE_FLAGS -cpp)

  add_test(build_${testname} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target mfix_${testname})

  add_test(NAME ${testname}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${testdir}/runtests.sh ${CMAKE_CURRENT_BINARY_DIR}/mfix_${testname}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${testdir})

  set_tests_properties ( ${testname}
    PROPERTIES DEPENDS build_${testname}
    )

  set_tests_properties ( ${testname}
    PROPERTIES ENVIRONMENT "MFIX=${CMAKE_BINARY_DIR}/mfix_${testname}"
    )

  target_link_libraries(mfix_${testname} model)
endforeach( testdir ${TESTS})
