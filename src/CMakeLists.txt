set(CXX_sources
   main.cpp
   mfix_evolve.cpp
   mfix_init.cpp
   mfix_level.cpp
   mfix_regrid.cpp )

set(CXX_includes
   mfix_level.H
   include/mfix_F.H)

set(F90_sources
   bc_mod.f90
   calc_cell.f90
   calc_mu_g.f90
   calc_tau.f90
   compute_divu.f90
   compute_vort.f90
   constant_mod.f90
   eos_mod.f90
   fill_bc0.f90
   fld_constants_mod.f90
   get_data.f90
   ic_mod.f90
   init_fluid.f90
   init_namelist.f90
   output_mod.f90
   param_mod.f90
   run_mod.f90
   scales_mod.f90
   set_bc0.f90
   set_bc_flow.f90
   set_bc_type.f90
   set_p0.f90
   set_ps.f90
   set_velocity_bcs.f90
   ps_mod.f90
   utilities_mod.f90
   zero_wall_norm_vel.f90
   amrex_to_mfix.F90)

#
# Initialization routines
#
#list ( APPEND F90_sources init/

#
# EB stuff
#
list ( APPEND CXX_sources  eb/mfix_eb.cpp )
list ( APPEND F90_sources  eb/get_eb_walls.f90 eb/compute_normals.f90 )
# general EB geometry tools (projections onto edges, etc...)
list ( APPEND F90_sources  eb/eb_geometry.f90)
# level-set sources
list ( APPEND CXX_sources  eb/eb_levelset.cpp )
list ( APPEND CXX_includes eb/eb_levelset.H )
list ( APPEND CXX_includes eb/mfix_eb_F.H )
list ( APPEND F90_sources  eb/eb_levelset.f90)

#
# Particles
#
list ( APPEND CXX_sources des/MFIXParticleContainer.cpp )
list ( APPEND CXX_includes des/MFIXParticleContainer.H )
list ( APPEND F90_sources des/calc_drag.f90 des/calc_particle_collisions.f90)
list ( APPEND F90_sources des/calc_wall_collisions.f90 des/cffctowall.f90 )
list ( APPEND F90_sources des/cfrelvel.f90 des/cfslide.f90 )
list ( APPEND F90_sources des/des_bc_mod.f90 des/des_init_namelist.f90 )
list ( APPEND F90_sources des/des_time_march.f90 des/discretelement_mod.f90 )
list ( APPEND F90_sources des/drag_gp_des.f90 des/init_collision.f90 )
list ( APPEND F90_sources des/particle_generator.f90 des/set_particle_properties.f90 )
list ( APPEND F90_sources des/stl_functions_des_mod.f90 )
list ( APPEND F90_sources des/particle_pic.f90 des/particle_mod.f90 des/set_drag_bcs.f90 )
list ( APPEND F90_sources des/drag_mod.f90 des/flip_particle_vol.f90 )
list ( APPEND F90_sources des/construct_gradp.f90)

#
# SIMPLE routines
#
list ( APPEND CXX_sources simple/mfix_simple.cpp )
list ( APPEND CXX_sources simple/mfix_linear_solve.cpp )
list ( APPEND CXX_sources simple/solve_bicgstab.cpp )
list ( APPEND F90_sources simple/leqsol_mod.f90  simple/matvec.f90 simple/solver_params.f90 )
list ( APPEND F90_sources simple/source_u_g.f90 simple/source_v_g.f90 )
list ( APPEND F90_sources simple/source_w_g.f90 simple/source_pp_g.f90 )
list ( APPEND F90_sources simple/tau_u_g.f90 simple/tau_v_g.f90 simple/tau_w_g.f90 )
list ( APPEND F90_sources simple/conv_dif_u_g.f90 simple/conv_dif_v_g.f90 )
list ( APPEND F90_sources simple/conv_dif_w_g.f90 simple/conv_pp_g.f90 simple/conv_rop.f90 )
list ( APPEND F90_sources simple/solve_vel_star.f90 simple/solve_pp_g.f90 simple/adjust_a_g.f90 )
list ( APPEND F90_sources simple/matrix_mod.f90 simple/adjust_dt_mod.f90 simple/calc_trd_g.f90 )
list ( APPEND F90_sources simple/calc_mflux.f90 simple/calc_d.f90 simple/check_convergence.f90 )
list ( APPEND F90_sources simple/residual_mod.f90 simple/correct_0.f90 simple/ur_facs_mod.f90 )
list ( APPEND F90_sources simple/functions_mod.f90 )
list ( APPEND F90_sources simple/set_bc1.f90 )

#
# Data checks
#
list ( APPEND F90_sources check_data/check_bc_flow.f90  check_data/check_bc_geometry.f90 )
list ( APPEND F90_sources check_data/check_bc_inflow.f90 check_data/check_bc_outflow.f90 )
list ( APPEND F90_sources check_data/check_bc_walls.f90 check_data/check_boundary_conditions.f90 )
list ( APPEND F90_sources check_data/check_collision_model.f90 check_data/check_gas_properties.f90 )
list ( APPEND F90_sources check_data/check_initial_conditions.f90 check_data/check_inputs.f90 )
list ( APPEND F90_sources check_data/check_particle_properties.f90 check_data/check_point_sources.f90 )
list ( APPEND F90_sources check_data/check_run_control.f90 )

#
# I/O routines
#
list ( APPEND CXX_sources io/mfixIO.cpp )
list ( APPEND F90_sources io/eb_to_pvd.f90 io/output_manager.f90 )
list ( APPEND F90_sources io/write_error.f90 io/write_out0.f90 )
list ( APPEND F90_sources io/read_namelist.f90 io/display_resid.f90 io/remove_comment.f90 )
list ( APPEND F90_sources io/parse_line.f90 io/deprecated.f90 io/error_manager_mod.f90 )
list ( APPEND F90_sources io/write_user.f90 )

#
# Users hooks
#
list ( APPEND F90_sources usr/usr0.f90 usr/usr0_des.f90 usr/usr1.f90 )
list ( APPEND F90_sources usr/usr1_des.f90 usr/usr2.f90 usr/usr2_des.f90 )
list ( APPEND F90_sources usr/usr3.f90 usr/usr3_des.f90 usr/usr_drag.f90 )
list ( APPEND F90_sources usr/usr_mod.f90 usr/write_usr0.f90 usr/write_usr1.f90 )

#
# Projection
#
list ( APPEND F90_sources projection/convection_mod.f90 projection/diffusion_mod.f90 )
list ( APPEND F90_sources projection/projection_mod.f90 projection/projection_tests_mod.f90 )
list ( APPEND F90_sources projection/slopes_mod.f90 )
list ( APPEND F90_sources projection/set_pressure_bcs.f90 projection/set_projection_bcs.f90 )
list ( APPEND CXX_sources projection/mfix_projection.cpp )
list ( APPEND CXX_sources projection/compute_dt.cpp )



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
