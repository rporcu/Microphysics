target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   EBGodunov.H
   mfix_ebgodunov_advection.cpp
   mfix_ebgodunov_corner_couple.H
   mfix_ebgodunov_plm.cpp
   mfix_ebgodunov_plm_fpu.cpp
   mfix_ebgodunov_predict.cpp
   mfix_ebgodunov_transverse_K.H
   mfix_ebgodunov_make_trans_vels.cpp
   )
