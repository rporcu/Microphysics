target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   redistribution.H
   mfix_redistribute.cpp
   mfix_flux_redistribute.cpp
   mfix_merge_redistribute.cpp
   mfix_state_redistribute.cpp
   )
