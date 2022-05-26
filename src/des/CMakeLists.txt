target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_particle_generator.H
   mfix_particle_generator.cpp
   mfix_des_drag_K.H
   mfix_des_rrates_K.H
   mfix_des_parts_gen_K.H
   mfix_set_gradp_bcs.cpp
   mfix_des_K.H
   mfix_pc.cpp
   mfix_pc.H
   mfix_pc_init.cpp
   mfix_pc_init_collision_params.cpp
   mfix_pc_remove_out_of_range.cpp
   mfix_pc_replicate.cpp
   mfix_pc_inflow.cpp
   mfix_calc_transfer_coeffs.cpp
   mfix_calc_volume_fraction.cpp
   mfix_calc_txfr.cpp
   )

