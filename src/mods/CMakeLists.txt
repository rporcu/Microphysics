target_include_directories(mfixcore
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   bc_mod.f90
   bc_mod_F.H
   eos_mod.H
   eos_mod.cpp
   utilities_mod.f90
   param_mod.f90
   param_mod_F.H
   param_mod_F.cpp
   run_mod.f90
   MFIX_FLUID_Parms.cpp
   MFIX_FLUID_Parms.H
   MFIX_DEM_Parms.cpp
   MFIX_DEM_Parms.H
   MFIX_IC_Parms.cpp
   MFIX_IC_Parms.H
   MFIX_BC_Parms.cpp
   MFIX_BC_Parms.H
   MFIX_REGIONS_Parms.cpp
   MFIX_REGIONS_Parms.H
   )
