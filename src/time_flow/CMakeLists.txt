target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE 
   mfix_advance.cpp
   mfix_compute_dt.cpp
   mfix_evolve.cpp
   )
