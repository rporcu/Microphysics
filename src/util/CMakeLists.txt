target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   bc_fill_nd.F90
   mfix_util_F.H
   mfix_util.cpp
   mfix_fillpatch.cpp
   interpolation_m.f90
   mfix_set_velocity_bcs.hpp
   mfix_set_velocity_bcs.cpp
   mfix_set_scalar_bcs.hpp
   mfix_set_scalar_bcs.cpp
   )
