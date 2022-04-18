target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources( mfixcore
   PRIVATE 
   mfix_calc_cell.H
   mfix_calc_cell.cpp
   mfix_init_fluid.H
   mfix_init_fluid.cpp
   mfix_init_redistribution.cpp
   mfix_set_ls_inflow.cpp
   mfix_set_p0.cpp
   mfix_set_gp0.cpp
   mfix_init.cpp
   mfix_arrays.cpp
   mfix_setup_solvers.cpp
   mfix_set_bc_type.cpp
   )
