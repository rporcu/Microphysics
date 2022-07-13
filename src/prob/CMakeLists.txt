target_include_directories(mfixcore
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_eb.H
   mfix_fluid.cpp
   mfix_fluid.H
   mfix_solids.cpp
   mfix_solids.H
   mfix_species.cpp
   mfix_species.H
   mfix_reactions.cpp
   mfix_reactions.H
   mfix_dem.cpp
   mfix_dem.H
   mfix_pic.cpp
   mfix_pic.H
   mfix_ic.cpp
   mfix_ic.H
   mfix_bc.cpp
   mfix_bc.H
   mfix_regions.cpp
   mfix_regions.H
   )
