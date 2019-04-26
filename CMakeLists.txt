cmake_minimum_required(VERSION 3.14)

project(MFIX-Exa
   # VERSION <version>
   DESCRIPTION "A multiphase modeling tool for the exascale"
   HOMEPAGE_URL "https://amrex-codes.github.io/MFIX-Exa/"
   LANGUAGES  C CXX Fortran
   )


# Set the search path for cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/tools/CMake )

# Load utilities module
include( MFIX_Utils )

find_package(AMReX 19.03.1.243 CONFIG QUIET)

if (AMReX_FOUND)

   message(STATUS "AMReX found: configuration file located at ${AMReX_DIR}")
   
   # Load required modules
   include( MFIX_Options )

   echo_amrex_config_options ()

   if ( NOT ( AMREX_ENABLE_PARTICLES AND AMREX_ENABLE_AMRDATA
            AND AMREX_ENABLE_EB) )
      message ( FATAL_ERROR "AMReX installation missing required components!")
   endif ()

   if (ENABLE_OMP AND (NOT AMREX_ENABLE_OMP) )
      message(FATAL_ERROR "Option ENABLE_OMP requires AMReX to be built with OpenMP support")
   endif ()

   if (ENABLE_MPI AND (NOT AMREX_ENABLE_MPI) )
      message(FATAL_ERROR "Option ENABLE_MPI requires AMReX to be built with MPI support")
   endif ()

   if (ENABLE_HYPRE AND (NOT AMREX_ENABLE_HYPRE) )
      message(FATAL_ERROR "Option ENABLE_HYPRE requires AMReX to be built with HYPRE support")
   endif ()

   if (ENABLE_CUDA AND (NOT AMREX_ENABLE_CUDA) )
      message(FATAL_ERROR "Option ENABLE_CUDA requires AMReX to be built with CUDA support")
   endif ()

   if ( NOT ( "${CMAKE_BUILD_TYPE}" STREQUAL "${AMREX_BUILD_TYPE}" ) )
      message (WARNING "MFIX build type (${CMAKE_BUILD_TYPE}) type does not match AMReX build type (${AMREX_BUILD_TYPE})")
   endif ()

   if (ENABLE_CUDA)
      setup_cuda()
   endif ()
   
else ()

   message(STATUS "No AMReX installation found: cloning AMReX repo")
   
   set(AMREX_GIT_REPO "https://github.com/AMReX-Codes/amrex.git" )
   set(AMREX_GIT_COMMIT_MASTER   4eb4e7a25050ca83f02e551fcb9b8a591834395 )
   set(AMREX_GIT_COMMIT_DEVELOP 6452b0c189128831deb0886765de8405b874e194 )

   # Amrex as sub-project
   include(FetchContent)

   set(ENABLE_TUTORIALS OFF)
   set(ENABLE_PARTICLES ON)
   set(ENABLE_EB ON)
   set(ENABLE_3D_NODAL_MLMG ON)
   set(ENABLE_AMRDATA ON)
   
   FetchContent_Declare(
      amrex
      GIT_REPOSITORY         ${AMREX_GIT_REPO}
      GIT_TAG                mr/cmake     
      )

   FetchContent_MakeAvailable(amrex)
   FetchContent_GetProperties(amrex SOURCE_DIR amrex_src_dir)
   get_directory_property(amrex_targets DIRECTORY ${amrex_src_dir}/Src BUILDSYSTEM_TARGETS)

   foreach (_target IN LISTS amrex_targets )
      add_library(AMReX::${_target} ALIAS ${_target})
   endforeach ()
   

   list(APPEND CMAKE_MODULE_PATH ${amrex_src_dir}/Tools/CMake)

   if (ENABLE_CUDA)
      include(AMReX_SetupCUDA)
   endif ()
  
endif ()


get_git_info ( )
check_build_tree_path ()
set_default_config_flags ()

# Set some cache variable for the whole project
set ( MFIX_LIBNAME "mfixcore" CACHE INTERNAL "Name of MFIX core library")
set ( MFIX_EXENAME "mfix"     CACHE INTERNAL "Name of MFIX executable")



set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp )


#
# Core MFIX library
#
add_subdirectory(src)

# Print out config summary
print_mfix_configuration_summary (${MFIX_LIBNAME})

# We should check what this line does
file(GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)


add_executable (${MFIX_EXENAME} ${MAIN} ${USR_OVERRIDES} )
target_link_libraries ( ${MFIX_EXENAME} ${MFIX_LIBNAME} )

if ( ENABLE_CUDA)
   set_source_files_properties(${MAIN} PROPERTIES LANGUAGE CUDA )
   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++14 -Xcompiler=-std=c++14")
   set_target_properties( ${MFIX_EXENAME}
      PROPERTIES
      CUDA_SEPARABLE_COMPILATION ON  # This add -dc flag
      )
endif ()

# Tools
set (MFIX_TOOLS_EXE ${PROJECT_BINARY_DIR}/tools/)
add_subdirectory (tools)

# Tests
enable_testing ()
add_subdirectory ( tests )

# Tags
include (MFIX_Tags)
add_tags_targets ()

# Add a target to generate API documentation with Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND)
   set(DOC_DIR ${CMAKE_SOURCE_DIR}/doc)
   configure_file(${DOC_DIR}/Doxyfile.in ${DOC_DIR}/Doxyfile)
   add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOC_DIR}/Doxyfile
      WORKING_DIRECTORY ${DOC_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
endif (DOXYGEN_FOUND)

# Copy compile_commands.json from the build directory to the project root directory
if ( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
endif ()
