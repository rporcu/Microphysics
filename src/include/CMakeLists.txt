target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources( mfixcore
   PRIVATE
   desnamelist.inc
   gas_phase.inc
   geometry.inc
   initial_conditions.inc
   species.inc
   species_indices.inc
   usr_hooks.inc
   )
