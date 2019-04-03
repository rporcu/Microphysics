target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   bc_fill_nd.F90
   divop_mod.f90
   compute_vort.f90
   set_velocity_bcs.f90
   set_scalar_bcs.f90
   mfix_util_F.H
   mfix_util_f.f90
   mfix_util.cpp
   mfix_fillpatch.cpp
   reconstruct_velocity.f90
   interpolation_m.f90
   )
