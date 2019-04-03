target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)


target_sources( mfixcore
   PRIVATE 
   diffusion_mod.f90
   diffusion_eb_mod.f90
   eb_wallflux_mod.f90
   mfix_compute_divtau.cpp
   mfix_diffuse_velocity.cpp
   mfix_diff_F.H
   )
