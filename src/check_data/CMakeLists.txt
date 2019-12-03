target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   check_bc_flow.f90
   check_bc_inflow.f90
   check_bc_outflow.f90
   check_boundary_conditions.f90
   check_collision_model.f90
   check_gas_properties.f90
   check_initial_conditions.f90
   check_inputs.f90
   check_particle_properties.f90
   )
