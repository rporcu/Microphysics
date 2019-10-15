target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_compute_divtau.cpp
   mfix_diffuse_velocity.cpp
   mfix_diff_F.H
   set_diff_bcs.f90
   )
