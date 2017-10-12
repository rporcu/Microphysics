set(CXX_sources
   main.cpp
   mfix_evolve.cpp
   mfix_init.cpp
   mfix_level.cpp
   mfix_linear_solve.cpp
   mfix_regrid.cpp
   solve_bicgstab.cpp
   MFIXParticleContainer.cpp
   io/mfixIO.cpp)

set(CXX_includes
   mfix_level.H
   MFIXParticleContainer.H
   include/mfix_F.H)

set(F90_sources
   adjust_a_g.f90
   adjust_dt_mod.f90
   bc_mod.f90
   calc_cell.f90
   calc_coeff.f90
   calc_d.f90
   calc_mflux.f90
   calc_mu_g.f90
   calc_ro_g.f90
   calc_trd_g.f90
   check_convergence.f90
   check_data/check_bc_flow.f90
   check_data/check_bc_geometry.f90
   check_data/check_bc_inflow.f90
   check_data/check_bc_outflow.f90
   check_data/check_bc_walls.f90
   check_data/check_boundary_conditions.f90
   check_data/check_collision_model.f90
   check_data/check_gas_properties.f90
   check_data/check_initial_conditions.f90
   check_data/check_inputs.f90
   check_data/check_particle_properties.f90
   check_data/check_point_sources.f90
   check_data/check_run_control.f90
   constant_mod.f90
   conv_dif_u_g.f90
   conv_dif_v_g.f90
   conv_dif_w_g.f90
   conv_pp_g.f90
   conv_rop.f90
   correct_0.f90
   deprecated.f90
   des/calc_collision_wall_mod.f90
   des/calc_drag.f90
   des/calc_force_dem.f90
   des/cffctowall.f90
   des/cfrelvel.f90
   des/cfslide.f90
   des/des_bc_mod.f90
   des/des_init_namelist.f90
   des/des_time_march.f90
   des/discretelement_mod.f90
   des/drag_gp_des.f90
   des/init_collision.f90
   des/set_particle_properties.f90
   des/stl_functions_des_mod.f90
   display_resid.f90
   drag_mod.f90
   eos_mod.f90
   error_manager_mod.f90
   fill_bc0.f90
   fld_constants_mod.f90
   flip_particle_vol.f90
   functions_mod.f90
   get_data.f90
   get_tunit.f90
   ic_mod.f90
   init_fluid.f90
   init_namelist.f90
   io/output_manager.f90
   io/write_error.f90
   io/write_out0.f90
   location_check.f90
   matrix_mod.f90
   output_mod.f90
   particle_mod.f90
   particle_pic.f90
   param_mod.f90
   parse_line.f90
   physical_prop.f90
   ps_mod.f90
   read_namelist.f90
   remove_comment.f90
   residual_mod.f90
   run_mod.f90
   scales_mod.f90
   set_bc0.f90
   set_bc1.f90
   set_bc_flow.f90
   set_bc_type.f90
   set_drag_bcs.f90
   set_ps.f90
   solve_pp_g.f90
   solve_vel_star.f90
   solver_params.f90
   solvers/leqsol_mod.f90
   solvers/matvec.f90
   source_pp_g.f90
   source_u_g.f90
   source_v_g.f90
   source_w_g.f90
   tau_u_g.f90
   tau_v_g.f90
   tau_w_g.f90
   ur_facs_mod.f90
   usr/usr0.f90
   usr/usr0_des.f90
   usr/usr1.f90
   usr/usr1_des.f90
   usr/usr2.f90
   usr/usr2_des.f90
   usr/usr3.f90
   usr/usr3_des.f90
   usr/usr_drag.f90
   usr/usr_mod.f90
   usr/write_usr0.f90
   usr/write_usr1.f90
   utilities_mod.f90
   zero_wall_norm_vel.f90
   amrex_to_mfix.F90)

#
# Add sources for projection method
# 
list ( APPEND F90_sources convection_mod.f90 diffusion_mod.f90 )


# Generate AMReX_buildInfo.cpp
add_custom_command(
   COMMAND  ${AMREX_INSTALL_PATH}/Tools/C_scripts/makebuildinfo_C.py
   --amrex_home "${AMREX_INSTALL_PATH}"
   --COMP ${CMAKE_C_COMPILER_ID} --COMP_VERSION ${CMAKE_C_COMPILER_VERSION}
   --FCOMP ${CMAKE_Fortran_COMPILER_ID} --FCOMP_VERSION ${CMAKE_C_COMPILER_VERSION}
   --GIT ". ${AMREX_GIT}"
   OUTPUT AMReX_buildInfo.cpp
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   COMMENT "Generating AMReX_buildInfo.cpp" )

add_library(mfixcore STATIC ${CXX_sources} ${F90_sources} ${CXX_includes} AMReX_buildInfo.cpp )

# Remove AMReX_buildInfo.cpp
add_custom_command(TARGET mfixcore POST_BUILD COMMAND rm AMReX_buildInfo.cpp
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   COMMENT "Removing AMReX_buildInfo.cpp" )

# Define sources for type checking
set (F90SRC ${F90_sources})
set (CXXINCLUDES ${CXX_includes})
include (MFIX_Typecheck)
