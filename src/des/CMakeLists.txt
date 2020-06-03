target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   particle_generator.H
   particle_generator.cpp
   des_drag_K.H
   mfix_set_gradp_bcs.cpp
   mfix_des_K.H
   MFIXParticleContainer.cpp
   MFIXParticleContainer.H
   MFIX_PC_Init.cpp
   MFIX_PC_InitCollisionParams.cpp
   MFIX_PC_RemoveOutOfRange.cpp
   MFIX_PC_Replicate.cpp
   mfix_calc_transfer_coeffs.cpp
   mfix_calc_volume_fraction.cpp
   mfix_calc_txfr.cpp
   )

