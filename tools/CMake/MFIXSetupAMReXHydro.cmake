#
# ~~~~~~~~~~~~~~~~~~~~~~~ SUPERBUILD MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
if (NOT EXISTS  "${PROJECT_SOURCE_DIR}/.git")
   message(FATAL_ERROR
      "${PROJECT_SOURCE_DIR} is not a Git repo: missing .git")
endif ()

set(AMREX_HYDRO_SRC_DIR "${PROJECT_SOURCE_DIR}/subprojects/AMReX-Hydro"
   CACHE INTERNAL "Path to AMReX-Hydro source (submodule)")

if (NOT EXISTS "${AMREX_HYDRO_SRC_DIR}/.git")
   message(STATUS "Initializing git submodule for AMReX-Hydro")

   find_package(Git REQUIRED)

   execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive subprojects/AMReX-Hydro
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT
      )

   if ( NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
   endif()

   unset(GIT_SUBMOD_RESULT)

endif ()

# Options
set(HYDRO_SPACEDIM              3)
set(HYDRO_MPI                   ${MFIX_MPI})
set(HYDRO_OMP                   ${MFIX_OMP})
set(HYDRO_GPU_BACKEND           ${MFIX_GPU_BACKEND} CACHE STRING "" FORCE)
set(HYDRO_EB                    ON)

# Add subdirectory to the build
add_subdirectory(${AMREX_HYDRO_SRC_DIR})
