target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_set_mac_velocity_bcs.cpp
   mfix_compute_fluxes.cpp
   mfix_compute_convective_term.cpp
   mfix_correct_small_cells.cpp
   mfix_slopes.cpp
   mfix_compute_MAC_projected_velocities.cpp
   )
