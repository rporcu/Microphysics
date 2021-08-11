target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_interp_K.H
   mfix_eb_interp_K.H
   mfix_eb_interp_shepard_K.H
   mfix_eb_newton_solver_K.H
   mfix_eb_node_connectivity_K.H
   )

