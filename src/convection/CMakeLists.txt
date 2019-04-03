target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   convection_mod.f90 
   mfix_mac_F.H
   slopes_mod.f90
   set_mac_velocity_bcs.f90
   MacProjection.H
   MacProjection.cpp
   mfix_convection.cpp
   mfix_slopes.cpp
   convection_eb_mod.f90
   )
