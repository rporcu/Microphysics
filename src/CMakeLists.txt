#
# Define a macro to add sources.
# This will also accumulate include directories if
# a header file is spotted
# WARNING: this automagically prepends the correct
# path to the source file.
#
set ( MFIX_INCLUDE_PATHS  ${CMAKE_CURRENT_LIST_DIR}/include )

macro ( add_sources )
  foreach ( item IN ITEMS ${ARGN} )
    set ( source_name  ${CMAKE_CURRENT_LIST_DIR}/${item} )
    string ( REPLACE "//" "/" source_name ${source_name})
    target_sources ( ${MFIX_LIBNAME} PRIVATE ${source_name} )

    # Add folder to list of includes if necessary
    get_filename_component ( filetype ${source_name} EXT )
    if ( (${filetype} MATCHES ".h") OR (${filetype} MATCHES ".H" ) )
      list (APPEND MFIX_INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR})
    endif()
  endforeach()
endmacro()


#
# Fist, define the library we want to add
#
add_library ( ${MFIX_LIBNAME} "" )

#
# Now, one by one, we add all the sources from all the subdirs
# to the library defined above
# 

#
# Uncategorized
# 
add_sources ( mfix_evolve.cpp mfix_init.cpp )
add_sources ( mfix_level.cpp mfix_regrid.cpp )
add_sources ( mfix_level.H mfix_F.H )
add_sources ( calc_cell.f90 calc_mu_g.f90 )
add_sources ( calc_tau.f90 compute_divu.f90 compute_vort.f90 )
add_sources ( fill_bc0.f90 )
add_sources ( get_data.f90 init_fluid.f90 init_namelist.f90 )
add_sources ( set_bc0.f90 set_bc_flow.f90 set_bc_type.f90 set_p0.f90 set_delp_dir.f90 )
add_sources ( set_ps.f90 set_velocity_bcs.f90 )
add_sources ( zero_wall_norm_vel.f90 amrex_to_mfix.F90 )

#
# Fortran Modules
#
include ( mods/CMakeLists.txt )

#
# EB stuff
#
include ( eb/CMakeLists.txt )

#
# Particles
#
include ( des/CMakeLists.txt )
include ( des_fluid/CMakeLists.txt )

#
# SIMPLE routines
#
add_sources ( simple/mfix_simple.cpp )
add_sources ( simple/mfix_linear_solve.cpp )
add_sources ( simple/solve_bicgstab.cpp )
add_sources ( simple/leqsol_mod.f90  simple/matvec.f90 simple/solver_params.f90 )
add_sources ( simple/source_u_g.f90 simple/source_v_g.f90 )
add_sources ( simple/source_w_g.f90 simple/source_pp_g.f90 )
add_sources ( simple/tau_u_g.f90 simple/tau_v_g.f90 simple/tau_w_g.f90 )
add_sources ( simple/conv_dif_u_g.f90 simple/conv_dif_v_g.f90 )
add_sources ( simple/conv_dif_w_g.f90 simple/conv_pp_g.f90 simple/conv_rop.f90 )
add_sources ( simple/solve_vel_star.f90 simple/solve_pp_g.f90 simple/adjust_a_g.f90 )
add_sources ( simple/matrix_mod.f90 simple/adjust_dt_mod.f90 simple/calc_trd_g.f90 )
add_sources ( simple/calc_mflux.f90 simple/calc_d.f90 simple/check_convergence.f90 )
add_sources ( simple/residual_mod.f90 simple/correct_0.f90 simple/ur_facs_mod.f90 )
add_sources ( simple/functions_mod.f90 )
add_sources ( simple/set_bc1.f90 simple/get_solver_params.f90 )

#
# Data checks
#
add_sources ( check_data/check_bc_flow.f90  check_data/check_bc_geometry.f90 )
add_sources ( check_data/check_bc_inflow.f90 check_data/check_bc_outflow.f90 )
add_sources ( check_data/check_bc_walls.f90 check_data/check_boundary_conditions.f90 )
add_sources ( check_data/check_collision_model.f90 check_data/check_gas_properties.f90 )
add_sources ( check_data/check_initial_conditions.f90 check_data/check_inputs.f90 )
add_sources ( check_data/check_particle_properties.f90 check_data/check_point_sources.f90 )
add_sources ( check_data/check_run_control.f90 )

#
# I/O routines
#
add_sources ( io/mfixIO.cpp )
add_sources ( io/eb_to_pvd.f90 io/output_manager.f90 )
add_sources ( io/write_error.f90 io/write_out0.f90 )
add_sources ( io/read_namelist.f90 io/display_resid.f90 io/remove_comment.f90 )
add_sources ( io/parse_line.f90 io/deprecated.f90 io/error_manager_mod.f90 )
add_sources ( io/write_user.f90 )

#
# Users hooks
#
add_sources ( usr/usr0.f90 usr/usr0_des.f90 usr/usr1.f90 )
add_sources ( usr/usr1_des.f90 usr/usr2.f90 usr/usr2_des.f90 )
add_sources ( usr/usr3.f90 usr/usr3_des.f90 usr/usr_drag.f90 )
add_sources ( usr/usr_mod.f90 usr/write_usr0.f90 usr/write_usr1.f90 )

#
# Projection
#
add_sources ( projection/convection_mod.f90 projection/diffusion_mod.f90 )
add_sources ( projection/projection_mod.f90 projection/projection_tests_mod.f90 )
add_sources ( projection/slopes_mod.f90 )
add_sources ( projection/set_pressure_bcs.f90 projection/set_projection_bcs.f90 )
add_sources ( projection/mfix_projection.cpp )
add_sources ( projection/compute_dt.cpp )

# Remove duplicates from include lists
list( REMOVE_DUPLICATES MFIX_INCLUDE_PATHS )

# 
# Add the include paths
# By adding the include paths as a property of mfixcore, we make sure
# that the these same paths will be used anytime another target will be
# list mficore as a dependency
# 
target_include_directories ( ${MFIX_LIBNAME} PUBLIC
  ${MFIX_INCLUDE_PATHS} ${MFIX_EXTRA_INCLUDE_PATH} )
