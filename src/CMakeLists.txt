set ( CXXSRC )
set ( F90SRC )
set ( ALLHEADERS )


#
# Uncategorized
# 
list ( APPEND CXXSRC  mfix_evolve.cpp mfix_init.cpp )
list ( APPEND CXXSRC  mfix_level.cpp mfix_regrid.cpp )

list ( APPEND ALLHEADERS mfix_level.H mfix_F.H )

list ( APPEND F90SRC calc_cell.f90 calc_mu_g.f90 )
list ( APPEND F90SRC calc_tau.f90 compute_divu.f90 compute_vort.f90 )
list ( APPEND F90SRC fill_bc0.f90 )
list ( APPEND F90SRC get_data.f90 init_fluid.f90 init_namelist.f90 )
list ( APPEND F90SRC set_bc0.f90 set_bc_flow.f90 set_bc_type.f90 set_p0.f90 set_delp_dir.f90 )
list ( APPEND F90SRC set_ps.f90 set_velocity_bcs.f90 )
list ( APPEND F90SRC zero_wall_norm_vel.f90 amrex_to_mfix.F90 )

#
# Fortran Modules
# 
list ( APPEND F90SRC mods/bc_mod.f90 mods/constant_mod.f90 mods/eos_mod.f90 )
list ( APPEND F90SRC mods/fld_constants_mod.f90 mods/ps_mod.f90 mods/ic_mod.f90 )
list ( APPEND F90SRC mods/utilities_mod.f90 mods/scales_mod.f90 )
list ( APPEND F90SRC mods/output_mod.f90 mods/param_mod.f90 mods/run_mod.f90 )

#
# EB stuff
#
list ( APPEND CXXSRC  eb/mfix_eb.cpp )
list ( APPEND F90SRC  eb/get_eb_walls.f90 eb/compute_normals.f90 )
# general EB geometry tools (projections onto edges, etc...)
list ( APPEND F90SRC  eb/eb_geometry.f90)
# level-set sources
list ( APPEND CXXSRC  eb/eb_levelset.cpp )
list ( APPEND CXX_includes eb/eb_levelset.H )
list ( APPEND CXX_includes eb/mfix_eb_F.H )
list ( APPEND F90SRC  eb/eb_levelset_F.f90)

#
# Particles
#
list ( APPEND CXXSRC des/MFIXParticleContainer.cpp )
list ( APPEND CXX_includes des/MFIXParticleContainer.H )
list ( APPEND CXX_includes des/mfix_des_F.H )
list ( APPEND F90SRC des/calc_drag.f90 des/calc_particle_collisions.f90)
list ( APPEND F90SRC des/calc_wall_collisions.f90 des/cffctowall.f90 )
list ( APPEND F90SRC des/cfrelvel.f90 des/cfslide.f90 )
list ( APPEND F90SRC des/des_bc_mod.f90 des/des_init_namelist.f90 )
list ( APPEND F90SRC des/des_time_march.f90 des/discretelement_mod.f90 )
list ( APPEND F90SRC des/drag_gp_des.f90 des/init_collision.f90 )
list ( APPEND F90SRC des/particle_generator.f90 des/set_particle_properties.f90 )
list ( APPEND F90SRC des/stl_functions_des_mod.f90 )
list ( APPEND F90SRC des/particle_pic.f90 des/particle_mod.f90 des/set_drag_bcs.f90 )
list ( APPEND F90SRC des/drag_mod.f90 des/flip_particle_vol.f90 )
list ( APPEND F90SRC des/construct_gradp.f90)

#
# SIMPLE routines
#
list ( APPEND CXXSRC simple/mfix_simple.cpp )
list ( APPEND CXXSRC simple/mfix_linear_solve.cpp )
list ( APPEND CXXSRC simple/solve_bicgstab.cpp )
list ( APPEND F90SRC simple/leqsol_mod.f90  simple/matvec.f90 simple/solver_params.f90 )
list ( APPEND F90SRC simple/source_u_g.f90 simple/source_v_g.f90 )
list ( APPEND F90SRC simple/source_w_g.f90 simple/source_pp_g.f90 )
list ( APPEND F90SRC simple/tau_u_g.f90 simple/tau_v_g.f90 simple/tau_w_g.f90 )
list ( APPEND F90SRC simple/conv_dif_u_g.f90 simple/conv_dif_v_g.f90 )
list ( APPEND F90SRC simple/conv_dif_w_g.f90 simple/conv_pp_g.f90 simple/conv_rop.f90 )
list ( APPEND F90SRC simple/solve_vel_star.f90 simple/solve_pp_g.f90 simple/adjust_a_g.f90 )
list ( APPEND F90SRC simple/matrix_mod.f90 simple/adjust_dt_mod.f90 simple/calc_trd_g.f90 )
list ( APPEND F90SRC simple/calc_mflux.f90 simple/calc_d.f90 simple/check_convergence.f90 )
list ( APPEND F90SRC simple/residual_mod.f90 simple/correct_0.f90 simple/ur_facs_mod.f90 )
list ( APPEND F90SRC simple/functions_mod.f90 )
list ( APPEND F90SRC simple/set_bc1.f90 simple/get_solver_params.f90 )

#
# Data checks
#
list ( APPEND F90SRC check_data/check_bc_flow.f90  check_data/check_bc_geometry.f90 )
list ( APPEND F90SRC check_data/check_bc_inflow.f90 check_data/check_bc_outflow.f90 )
list ( APPEND F90SRC check_data/check_bc_walls.f90 check_data/check_boundary_conditions.f90 )
list ( APPEND F90SRC check_data/check_collision_model.f90 check_data/check_gas_properties.f90 )
list ( APPEND F90SRC check_data/check_initial_conditions.f90 check_data/check_inputs.f90 )
list ( APPEND F90SRC check_data/check_particle_properties.f90 check_data/check_point_sources.f90 )
list ( APPEND F90SRC check_data/check_run_control.f90 )

#
# I/O routines
#
list ( APPEND CXXSRC io/mfixIO.cpp )
list ( APPEND F90SRC io/eb_to_pvd.f90 io/output_manager.f90 )
list ( APPEND F90SRC io/write_error.f90 io/write_out0.f90 )
list ( APPEND F90SRC io/read_namelist.f90 io/display_resid.f90 io/remove_comment.f90 )
list ( APPEND F90SRC io/parse_line.f90 io/deprecated.f90 io/error_manager_mod.f90 )
list ( APPEND F90SRC io/write_user.f90 )

#
# Users hooks
#
list ( APPEND F90SRC usr/usr0.f90 usr/usr0_des.f90 usr/usr1.f90 )
list ( APPEND F90SRC usr/usr1_des.f90 usr/usr2.f90 usr/usr2_des.f90 )
list ( APPEND F90SRC usr/usr3.f90 usr/usr3_des.f90 usr/usr_drag.f90 )
list ( APPEND F90SRC usr/usr_mod.f90 usr/write_usr0.f90 usr/write_usr1.f90 )

#
# Projection
#
list ( APPEND F90SRC projection/convection_mod.f90 projection/diffusion_mod.f90 )
list ( APPEND F90SRC projection/projection_mod.f90 projection/projection_tests_mod.f90 )
list ( APPEND F90SRC projection/slopes_mod.f90 )
list ( APPEND F90SRC projection/set_pressure_bcs.f90 projection/set_projection_bcs.f90 )
list ( APPEND CXXSRC projection/mfix_projection.cpp )
list ( APPEND CXXSRC projection/compute_dt.cpp )


# # Accumulate sources
# set ( ALLSRC ${CXXSRC} ${F90SRC} )

# #
# # Add absolute path (needed by target_sources)
# # 
# prepend ( CXXSRC     ${CMAKE_CURRENT_LIST_DIR} )
# prepend ( F90SRC     ${CMAKE_CURRENT_LIST_DIR} )
# prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR} )


# Add MFIX core library
add_library ( mfixcore STATIC ${CXXSRC} ${F90SRC} ${ALLHEADERS} ) 


# Define sources for type checking
set ( CXXINCLUDES ${ALLHEADERS} )
include ( MFIX_Typecheck )
