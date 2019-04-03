target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   gas_phase.inc
   initial_conditions.inc
   physical_params.inc
   run_control.inc
   species_indices.inc
   usrnlst.inc
   desnamelist.inc
   geometry.inc
   numerical_params.inc
   point_sources.inc
   species.inc
   usr_hooks.inc
   )
