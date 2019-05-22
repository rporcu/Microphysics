target_include_directories(mfixcore
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   bc_mod.f90
   bc_mod_F.H
   constant_mod.f90
   eos_mod.f90
   fld_constants_mod.f90
   ps_mod.f90
   ic_mod.f90
   utilities_mod.f90
   scales_mod.f90
   output_mod.f90
   param_mod.f90
   run_mod.f90
   )
