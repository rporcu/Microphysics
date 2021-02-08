target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   Godunov.H
   mfix_godunov_predict.cpp
   mfix_godunov_advection.cpp
   mfix_godunov_plm.cpp
   mfix_godunov_ppm.cpp
   mfix_godunov_plm.H
   mfix_godunov_ppm.H
   mfix_godunov_corner_couple.H
   mfix_godunov_trans_bc.H
   )
