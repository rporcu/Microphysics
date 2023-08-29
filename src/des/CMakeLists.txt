target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/dem>)
target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/deposition>)
target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/interpolation>)
target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/pic>)

target_sources( mfixcore
   PRIVATE
   mfix_particle_generator.H
   mfix_particle_generator.cpp
   mfix_des_drag_K.H
   mfix_des_parts_gen_K.H
   mfix_des_K.H
   mfix_des_conv_coeff_K.H
   mfix_pc_updates_K.H
   mfix_pc.cpp
   mfix_pc.H
   mfix_pc_init.cpp
   mfix_pc_remove_out_of_range.cpp
   mfix_pc_replicate.cpp
   mfix_pc_inflow.cpp
   mfix_calc_transfer_coeffs.cpp
   mfix_calc_volume_fraction.cpp
   mfix_calc_txfr.cpp
   mfix_pc_mass_balance.cpp
   mfix_pc_impose_mean.cpp
   # dem ----------------------------------------------------------------------
   dem/mfix_des_coll_K.H
   dem/mfix_des_cond_K.H
   dem/mfix_evolve_particles.cpp
   dem/mfix_pc_init_collision_params.cpp
   dem/mfix_pc_interactions_K.H
   # deposition ---------------------------------------------------------------
   deposition/mfix_deposition_K.H
   deposition/mfix_deposition_bcs.cpp
   deposition/mfix_deposition_op.cpp
   deposition/mfix_deposition_op.H
   deposition/mfix_redistribute_deposition.cpp
   # interpolation ------------------------------------------------------------
   interpolation/mfix_interp_K.H
   interpolation/mfix_eb_interp_K.H
   interpolation/mfix_eb_interp_shepard_K.H
   interpolation/mfix_eb_newton_solver_K.H
   interpolation/mfix_eb_node_connectivity_K.H
   # pic ----------------------------------------------------------------------
   pic/mfix_evolve_parcels.cpp
   pic/mfix_calc_avg_solids_prop.cpp
   pic/mfix_calc_solids_stress.cpp
   pic/mfix_pic_iteration.cpp
   pic/mfix_pc_advance_pic.cpp
   pic/mfix_pc_deposit_pic.cpp
   pic/mfix_pic_K.H
   )
