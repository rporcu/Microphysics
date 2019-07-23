target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE 
   calc_particle_collisions.f90
   calc_wall_collisions.f90
   cfrelvel.f90
   des_init_namelist.f90
   des_time_march.f90
   discretelement_mod.f90
   discretelement_mod_F.H
   init_collision.f90
   particles_generator.hpp
   particles_generator.cpp
   particle_mod.f90
   des_drag_K.H
   )
