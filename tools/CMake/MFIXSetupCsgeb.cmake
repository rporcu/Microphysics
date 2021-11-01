if(NOT MFIX_CSG)
  return()
endif()

if (NOT EXISTS  "${PROJECT_SOURCE_DIR}/.git")
   message(FATAL_ERROR
      "${PROJECT_SOURCE_DIR} is not a Git repo: missing .git")
endif ()

set(CSGEB_SRC_DIR "${PROJECT_SOURCE_DIR}/subprojects/csg-eb"
   CACHE INTERNAL "Path to CSG-EB source (submodule)")

message(STATUS "Initializing git submodule for csg-eb")

find_package(Git REQUIRED)

execute_process(
   COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive subprojects/csg-eb
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
   RESULT_VARIABLE GIT_SUBMOD_RESULT
   )

if ( NOT GIT_SUBMOD_RESULT EQUAL "0")
   message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
endif()

unset(GIT_SUBMOD_RESULT)

# Add subdirectory to the build
add_subdirectory(${CSGEB_SRC_DIR})
