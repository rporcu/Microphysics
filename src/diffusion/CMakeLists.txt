target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   DiffusionOp.cpp
   DiffusionOp.H
   mfix_diffuse_drag.cpp
   mfix_diffuse_velocity.cpp
   mfix_diffuse_scalar.cpp
   mfix_diffuse_volfrac.cpp
   )
