cmake_minimum_required(VERSION 2.8)
project("MFIX" Fortran)

include_directories(model/include)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)

option(MPI "MPI support" OFF)

if ( MPI )
  find_package ( MPI )
  if ( MPI_Fortran_COMPILER )
    include_directories( ${MPI_INCLUDE_PATH} )
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  else ( MPI_Fortran_COMPILER )
    message (FATAL_ERROR " MPI support requested but MPI compiler not found ")
  endif ( MPI_Fortran_COMPILER )
endif ( MPI )

set(SOURCES
  model/accum_resid.f
  model/adjust_a_g.f
  model/adjust_dt_mod.f
  model/adjust_leq.f
  model/allocate_arrays.f
  model/bc_mod.f
  model/bodyforce_mod.f
  model/boundfunijk_mod.f
  model/calc_cell.f
  model/calc_coeff.f
  model/calc_d_e.f
  model/calc_d_n.f
  model/calc_d_t.f
  model/calc_mflux.f
  model/calc_mu_g.f
  model/calc_outflow.f
  model/calc_resid.f
  model/calc_trd_g.f
  model/cartesian_grid/CG_set_bc0.f
  model/cartesian_grid/CG_set_outflow.f
  model/cartesian_grid/CG_source_u_g.f
  model/cartesian_grid/CG_source_v_g.f
  model/cartesian_grid/CG_source_w_g.f
  model/cartesian_grid/allocate_cut_cell_arrays.f
  model/cartesian_grid/allocate_dummy_cut_cell_arrays.f
  model/cartesian_grid/calc_vort_out.f
  model/cartesian_grid/cartesian_grid_init_namelist.f
  model/cartesian_grid/check_data_cartesian.f
  model/cartesian_grid/cut_cell_preprocessing.f
  model/cartesian_grid/cutcell_mod.f
  model/cartesian_grid/dashboard_mod.f
  model/cartesian_grid/deallocate_cut_cell_arrays.f
  model/cartesian_grid/define_quadrics.f
  model/cartesian_grid/dmp_cartesian.f
  model/cartesian_grid/eval_usr_fct.f
  model/cartesian_grid/get_alpha.f
  model/cartesian_grid/get_connectivity.f
  model/cartesian_grid/get_cut_cell_flags.f
  model/cartesian_grid/get_cut_cell_volume_area.f
  model/cartesian_grid/get_delh.f
  model/cartesian_grid/get_master.f
  model/cartesian_grid/get_poly_data.f
  model/cartesian_grid/get_stl_data.f
  model/cartesian_grid/polygon_mod.f
  model/cartesian_grid/progress_bar_mod.f
  model/cartesian_grid/quadric_mod.f
  model/cartesian_grid/set_Odxyz.f
  model/cartesian_grid/stl_mod.f
  model/cartesian_grid/update_dashboard.f
  model/cartesian_grid/vtk_mod.f
  model/cartesian_grid/vtk_out.f
  model/cartesian_grid/write_progress_bar.f
  model/check_ab_m.f
  model/check_bqend.f
  model/check_convergence.f
  model/check_data/check_axis.f
  model/check_data/check_bc_dem.f
  model/check_data/check_bc_geometry.f
  model/check_data/check_bc_inflow.f
  model/check_data/check_bc_outflow.f
  model/check_data/check_bc_walls.f
  model/check_data/check_boundary_conditions.f
  model/check_data/check_dmp_prereqs.f
  model/check_data/check_gas_phase.f
  model/check_data/check_geometry.f
  model/check_data/check_geometry_prereqs.f
  model/check_data/check_ic_common_discrete.f
  model/check_data/check_ic_dem.f
  model/check_data/check_initial_conditions.f
  model/check_data/check_numerics.f
  model/check_data/check_output_control.f
  model/check_data/check_point_sources.f
  model/check_data/check_run_control.f
  model/check_data/check_solids_common_all.f
  model/check_data/check_solids_common_discrete.f
  model/check_data/check_solids_dem.f
  model/check_data/check_solids_model_prereqs.f
  model/check_data/check_solids_phases.f
  model/check_data_20.f
  model/check_data_30.f
  model/check_plane.f
  model/constant_mod.f
  model/conv_dif_u_g.f
  model/conv_dif_v_g.f
  model/conv_dif_w_g.f
  model/conv_pp_g.f
  model/conv_rop.f
  model/corner.f
  model/correct_0.f
  model/dbg_mod.f
  model/deprecated.f
  model/des/calc_collision_wall_mod.f
  model/des/calc_drag_des.f
  model/des/calc_epg_des.f
  model/des/calc_force_dem.f
  model/des/calc_grad_des.f
  model/des/calc_interp_weights.f
  model/des/calc_pg_grad.f
  model/des/cfassign.f
  model/des/cffctowall.f
  model/des/cfnewvalues.f
  model/des/cfrelvel.f
  model/des/cfslide.f
  model/des/cfupdateold.f
  model/des/check_cell_movement.f
  model/des/comp_mean_fields.f
  model/des/comp_mean_fields0.f
  model/des/comp_mean_fields1.f
  model/des/des_allocate_mod.f
  model/des/des_bc_mod.f
  model/des/des_granular_temperature.f
  model/des/des_init_arrays.f
  model/des/des_init_namelist.f
  model/des/des_time_march.f
  model/des/desgrid_mod.f
  model/des/desmpi_mod.f
  model/des/desmpi_wrapper_mod.f
  model/des/discretelement_mod.f
  model/des/drag_gp_des.f
  model/des/drag_gs_des0.f
  model/des/drag_gs_des1.f
  model/des/gas_drag.f
  model/des/generate_particles_mod.f
  model/des/interpolation_mod.f
  model/des/layout_mi_dem.f
  model/des/make_arrays_des.f
  model/des/mass_inflow_dem.f
  model/des/mass_outflow_dem.f
  model/des/mpi_comm_des_mod.f
  model/des/mpi_funs_des_mod.f
  model/des/mpi_init_des_mod.f
  model/des/mpi_node_des_mod.f
  model/des/mpi_pack_des_mod.f
  model/des/mpi_unpack_des_mod.f
  model/des/neighbour.f
  model/des/nsquare.f
  model/des/particle_filter_mod.f
  model/des/particles_in_cell.f
  model/des/randomno_mod.f
  model/des/read_particle_input.f
  model/des/read_res0_des.f
  model/des/read_res1_des_mod.f
  model/des/sendrecvnode_mod.f
  model/des/set_bc_dem.f
  model/des/set_bc_dem_mi.f
  model/des/set_bc_dem_mo.f
  model/des/set_filter_des.f
  model/des/set_geometry_des.f
  model/des/set_phase_index.f
  model/des/solid_drag.f
  model/des/stl_dbg_des_mod.f
  model/des/stl_functions_des_mod.f
  model/des/stl_preproc_des_mod.f
  model/des/vtp_mod.f
  model/des/write_des_data.f
  model/des/write_res0_des.f
  model/des/write_res1_des_mod.f
  model/discretization_mod.f
  model/display_resid.f
  model/dmp_modules/compar_mod.f
  model/dmp_modules/dbg_util_mod.f
  model/dmp_modules/debug_mod.f
  model/dmp_modules/gridmap_mod.f
  model/dmp_modules/mpi_mod.f
  model/dmp_modules/mpi_utility_mod.f
  model/dmp_modules/parallel_mpi_mod.f
  model/dmp_modules/sendrecv3_mod.f
  model/dmp_modules/sendrecv_mod.f
  model/drag_mod.f
  model/eos_mod.f
  model/error_manager_mod.f
  model/exit.f
  model/fldvar_mod.f
  model/flow_to_vel.f
  model/fun_avg_mod.f
  model/functions_mod.f
  model/geometry_mod.f
  model/get_bc_area.f
  model/get_data.f
  model/get_ps.f
  model/ic_mod.f
  model/indices_mod.f
  model/init_fvars.f
  model/init_namelist.f
  model/io/funits_mod.f
  model/io/in_binary_512_mod.f
  model/io/in_binary_512i_mod.f
  model/io/open_file.f
  model/io/open_files.f
  model/io/out_array.f
  model/io/out_array_c.f
  model/io/out_array_k.f
  model/io/out_array_kc.f
  model/io/out_bin_512.f
  model/io/out_bin_512i.f
  model/io/out_bin_512r.f
  model/io/out_bin_r.f
  model/io/output_manager.f
  model/io/read_res1.f
  model/io/write_ab_m.f
  model/io/write_ab_m_var.f
  model/io/write_error.f
  model/io/write_header.f
  model/io/write_out0.f
  model/io/write_out1.f
  model/io/write_out3.f
  model/io/write_res0.f
  model/io/write_res1.f
  model/io/write_spx0.f
  model/io/write_spx1.f
  model/io/write_table.f
  model/iterate.f
  model/location_check.f
  model/machine_mod.f
  model/matrix_mod.f
  model/mflux_mod.f
  model/mod_bc_i.f
  model/mod_bc_j.f
  model/mod_bc_k.f
  model/output_mod.f
  model/param1_mod.f
  model/param_mod.f
  model/parse_line.f
  model/parse_resid_string.f
  model/physical_prop.f
  model/physprop_mod.f
  model/ps_mod.f
  model/read_namelist.f
  model/remove_comment.f
  model/reset_new.f
  model/residual_mod.f
  model/run_mod.f
  model/scales_mod.f
  model/set_bc0.f
  model/set_bc1.f
  model/set_bc_flow.f
  model/set_constants.f
  model/set_constprop.f
  model/set_flags.f
  model/set_fluidbed_p.f
  model/set_geometry.f
  model/set_geometry1.f
  model/set_ic.f
  model/set_icbc_flags.f
  model/set_increments.f
  model/set_index1.f
  model/set_index1a.f
  model/set_max2.f
  model/set_outflow.f
  model/set_param.f
  model/set_ps.f
  model/set_ro_g.f
  model/set_wall_bc.f
  model/shift_dxyz.f
  model/solve_pp_g.f
  model/solve_vel_star.f
  model/solvers/BLAS.f
  model/solvers/DGTSV.f
  model/solvers/dgtsl.f
  model/solvers/leq_bicgs.f
  model/solvers/leq_bicgst.f
  model/solvers/leq_cg.f
  model/solvers/leq_gmres.f
  model/solvers/leq_sor.f
  model/solvers/leqsol_mod.f
  model/solvers/solve_lin_eq.f
  model/solvers/test_lin_eq.f
  model/solvers/xerbla.f
  model/source_pp_g.f
  model/source_u_g.f
  model/source_v_g.f
  model/source_w_g.f
  model/tau_u_g.f
  model/tau_v_g.f
  model/tau_w_g.f
  model/time_cpu_mod.f
  model/time_march.f
  model/tmp_array1_mod.f
  model/tmp_array_mod.f
  model/toleranc_mod.f
  model/under_relax.f
  model/update_old.f
  model/ur_facs_mod.f
  model/usr/usr0.f
  model/usr/usr0_des.f
  model/usr/usr1.f
  model/usr/usr1_des.f
  model/usr/usr2.f
  model/usr/usr2_des.f
  model/usr/usr3.f
  model/usr/usr3_des.f
  model/usr/usr_drag.f
  model/usr/usr_init_namelist.f
  model/usr/usr_mod.f
  model/usr/usr_physical_prop.f
  model/usr/usr_write_out0.f
  model/usr/usr_write_out1.f
  model/usr/write_usr0.f
  model/usr/write_usr1.f
  model/utilities_mod.f
  model/vavg_u_g.f
  model/vavg_v_g.f
  model/vavg_w_g.f
  model/vf_gs.f
  model/xsi_array_mod.f
  model/xsi_mod.f
  model/zero_norm_vel.f)

FILE(GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f)
add_executable(mfix model/mfix.f ${USR_OVERRIDES})
add_library(model STATIC ${SOURCES})

set_source_files_properties(${SOURCES} ${USR_OVERRIDES} PROPERTIES COMPILE_FLAGS -cpp)

set_target_properties(mfix PROPERTIES Fortran_FORMAT "FREE")
set_target_properties(model PROPERTIES Fortran_FORMAT "FREE")

target_link_libraries(mfix model)

enable_testing()

set(TESTS
  DEM01
  DEM02
  DEM03
  DEM04
  DEM05
  DEM06
  FLD01
  FLD02
  # FLD03
  # box3d
  # cyclone
  # fluidbed_DES
  )

foreach( test ${TESTS})
  FILE(GLOB USR_OVERRIDES tests/${test}/*.f)
  add_executable(mfix_${test} EXCLUDE_FROM_ALL model/mfix.f ${USR_OVERRIDES})

  set_target_properties(mfix_${test} PROPERTIES Fortran_FORMAT "FREE")
  set_source_files_properties(${USR_OVERRIDES} PROPERTIES COMPILE_FLAGS -cpp)

  add_test(buildtest_${test} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target mfix_${test})

  add_test(NAME test_${test}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/${test}/runtests.sh ${CMAKE_CURRENT_BINARY_DIR}/mfix_${test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/${test})

  set_tests_properties ( test_${test} PROPERTIES DEPENDS buildtest_${test})

  target_link_libraries(mfix_${test} model)
endforeach( test ${TESTS})
