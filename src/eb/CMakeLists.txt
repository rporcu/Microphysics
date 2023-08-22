target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_eb.cpp
   mfix_eb_box.cpp
   mfix_eb_cylinder.cpp
   mfix_eb_hopper.cpp
   mfix_eb_regular.cpp
   mfix_eb_generic.cpp
   mfix_eb_if.H
   mfix_get_walls.cpp
   mfix_write_eb_surface.cpp
   )

if (MFIX_CSG)
   target_sources(mfixcore PRIVATE mfix_eb_csg.cpp)
   if (MFIX_TESTS AND MFIX_CSG_TESTS)
      add_subdirectory(tests)
   endif()
   target_compile_definitions(mfixcore PUBLIC CSG_EB)
endif()
