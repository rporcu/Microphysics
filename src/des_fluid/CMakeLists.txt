target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_set_gradp_bcs.cpp
   get_collision_coefficients.f90 
   mfix_des_F.H
   mfix_des_K.H
   mfix_drag_K.H
   MFIXParticleContainer.cpp
   MFIXParticleContainer.H
   MFIX_DEM_Parms.cpp
   MFIX_DEM_Parms.H
   MFIX_PC_Init.cpp
   MFIX_PC_RemoveOutOfRange.cpp
   MFIX_PC_Replicate.cpp
   )

