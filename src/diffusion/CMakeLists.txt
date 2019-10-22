target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   DiffusionOp.cpp
   DiffusionOp.H
   diffusion_F.H
   set_diff_bcs.f90
   )
