target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources( mfixcore
   PRIVATE 
   amrex_to_mfix.F90
   calc_cell.f90 calc_mu_g.f90
   get_data.f90
   init_fluid.f90 init_namelist.f90
   set_bc0.f90
   set_ls_inflow.f90
   set_bc_type.f90
   set_delp_dir.f90 set_p0.f90 set_ps.f90
   set_ppe_bcs.f90
   mfix_init.cpp
   mfix_arrays.cpp
   mfix_setup_nodal_solver.cpp
   mfix_proj_F.H
   )
