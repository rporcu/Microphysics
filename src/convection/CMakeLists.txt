target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_divop_conv.cpp
   mfix_divop_conv.hpp
   mfix_set_mac_velocity_bcs.cpp
   mfix_compute_ugradu.cpp
   mfix_predict_vels_on_faces.cpp
   mfix_MAC_projection.cpp
   mfix_slopes.cpp
   )


