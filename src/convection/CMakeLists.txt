target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_mac_F.H
   set_mac_velocity_bcs.f90
   mfix_compute_ugradu.cpp
   mfix_compute_MAC_velocities.cpp
   mfix_MAC_projection.cpp
   mfix_slopes.cpp
   )


