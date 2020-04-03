target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   calc_particle_collisions.f90
   calc_wall_collisions.f90
   cfrelvel.f90
   des_time_march.f90
   discretelement_mod.f90
   init_collision.f90
   particle_generator.H
   particle_generator.cpp
   particle_mod.f90
   des_drag_K.H
   mfix_set_gradp_bcs.cpp
   get_collision_coefficients.f90 
   mfix_des_F.H
   mfix_des_K.H
   mfix_drag_K.H
   mfix_eb_interp_K.H
   MFIXParticleContainer.cpp
   MFIXParticleContainer.H
   MFIX_PC_Init.cpp
   MFIX_PC_RemoveOutOfRange.cpp
   MFIX_PC_Replicate.cpp
   mfix_calc_particle_beta.cpp
   mfix_calc_volume_fraction.cpp
   mfix_drag.cpp
   )

