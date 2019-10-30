target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_diffuse_eps.cpp
   mfix_diffuse_drag.cpp
   mfix_deposition_bcs.cpp
   MFIX_PC_Trilinear.cpp
   )
