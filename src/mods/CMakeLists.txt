target_include_directories(mfixcore
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   bc_mod.f90
   bc_mod_F.H
   constant_mod.f90
   eos_mod.f90
   eos_mod.H
   eos_mod.cpp
   fld_constants_mod.f90
   fld_constants_mod_F.H
   ic_mod.f90
   ic_mod_F.H
   utilities_mod.f90
   output_mod.f90
   param_mod.f90
   param_mod_F.H
   param_mod_F.cpp
   run_mod.f90
   )
