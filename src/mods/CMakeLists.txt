target_include_directories(mfixcore
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   eos_mod.H
   eos_mod.cpp
   MFIX_BC_List.H
   MFIX_EB_Parms.H
   MFIX_EB_Parms.cpp
   MFIX_FLUID_Parms.cpp
   MFIX_FLUID_Parms.H
   MFIX_SPECIES_Parms.cpp
   MFIX_SPECIES_Parms.H
   MFIX_DEM_Parms.cpp
   MFIX_DEM_Parms.H
   MFIX_PIC_Parms.cpp
   MFIX_PIC_Parms.H
   MFIX_IC_Parms.cpp
   MFIX_IC_Parms.H
   MFIX_BC_Parms.cpp
   MFIX_BC_Parms.H
   MFIX_REGIONS_Parms.cpp
   MFIX_REGIONS_Parms.H
   )
