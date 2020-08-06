target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE 
   mfix_apply_predictor.cpp
   mfix_apply_corrector.cpp
   mfix_advance.cpp
   mfix_apply_predictor.cpp
   mfix_apply_corrector.cpp
   mfix_compute_dt.cpp
   mfix_fluid_equations_rhs.cpp
   mfix_evolve.cpp
   mfix_evolve_fluid.cpp
   )
