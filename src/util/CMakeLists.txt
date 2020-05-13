target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   MFIX_FilCC.H
   MFIX_FilCC.cpp
   mfix_util.cpp
   mfix_fillpatch.cpp
   mfix_set_velocity_bcs.cpp
   mfix_set_epg_bcs.cpp
   mfix_set_enthalpy_bcs.cpp
   mfix_set_temperature_bcs.cpp
   mfix_set_density_bcs.cpp
   mfix_set_scalar_bcs.cpp
   mfix_set_specific_heat_bcs.cpp
   mfix_set_tracer_bcs.cpp
   mfix_set_viscosity_bcs.cpp
   MFIX_MFHelpers.H
   MFIX_MFHelpers.cpp
   MFIX_gradient.H
   MFIX_gradient.cpp
   )
