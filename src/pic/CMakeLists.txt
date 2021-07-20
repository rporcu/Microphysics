target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_evolve_parcels.cpp
   mfix_calc_avg_solids_prop.cpp
   mfix_calc_solids_stress.cpp
   mfix_predict_volume_fraction.cpp
   mfix_pc_advance_pic.cpp
   mfix_pc_deposit_pic.cpp
   mfix_pc_impose_walls.cpp
   mfix_pic_K.H
   )

