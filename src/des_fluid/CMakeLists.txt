target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   calc_particle_beta.f90
   flip_particle_vol.f90
   set_gradp_bcs.f90 
   get_collision_coefficients.f90 
   mfix_des_F.H
   mfix_des_K.H
   MFIXParticleContainer.cpp
   MFIXParticleContainer.H
   MFIX_DEM_Parms.cpp
   MFIX_DEM_Parms.H
   MFIX_PC_RemoveOutOfRange.cpp
   )
