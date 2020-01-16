target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources( mfixcore
   PRIVATE 
   amrex_to_mfix.F90
   calc_cell.f90
   calc_cell_F.H
   mfix_calc_mu_g.H
   mfix_calc_mu_g.cpp
   mfix_init_fluid.H
   mfix_init_fluid.cpp
   mfix_set_bc0.cpp
   mfix_set_ls_inflow.cpp
   mfix_set_p0.cpp
   mfix_set_gp0.cpp
   mfix_init.cpp
   mfix_arrays.cpp
   mfix_setup_solvers.cpp
   mfix_set_bc_type.cpp
   )
