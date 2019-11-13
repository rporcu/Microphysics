target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_deposition_K.H
   mfix_deposition_bcs.cpp
   mfix_diffuse_array.cpp
   mfix_diffuse_scalar.cpp
   MFIX_PC_Deposition.cpp
   mfix_redistribute_deposition.cpp
   )
