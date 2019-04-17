cmake_minimum_required (VERSION 3.14)

#
# Set the search path for cmake modules
#
set ( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/tools/CMake )

#
# Here we set the only option that decide if superbuild
# is enabled or not
#
set ( AMREX_INSTALL_DIR "" CACHE PATH
   "Path to installation directory (leave empty for superbuild)")

set ( SUPERBUILD  ON )

if ( AMREX_INSTALL_DIR )
   set ( SUPERBUILD OFF )
endif ()

#
# Now decide if invoking superbuild. If so, this listfile will be
# re-read by the MFIX_Superbuild.cmake
# Also add targets to clean sub-projects
#
if ( SUPERBUILD )
  if ( SUBMOD )
    include( MFIX_Superbuild2 )
  else()
    include( MFIX_Superbuild )
 endif ()

   add_custom_target(clean_mfix
      COMMAND make clean
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/mfix
      COMMENT "Cleaning mfix"
      )

   add_custom_target(clean_amrex
      COMMAND cd builddir && make clean
      COMMAND rm -rf installdir/*
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/amrex
      COMMENT "Cleaning amrex"
      )

   add_custom_target(realclean)
   add_dependencies(realclean clean_mfix clean_amrex)
   return ()
endif ()

#
# This is where the MFIX-standalone config starts
#
project (MFIX-Exa)

# Define the languages used by the project
enable_language (C)
enable_language (CXX)
enable_language (Fortran)

# Load required modules
include ( MFIX_Options )
include ( MFIX_Utils )

get_git_info ( )
check_build_tree_path ()
set_default_config_flags ()

# Set some cache variable for the whole project
set ( MFIX_LIBNAME "mfixcore" CACHE INTERNAL "Name of MFIX core library")
set ( MFIX_EXENAME "mfix"     CACHE INTERNAL "Name of MFIX executable")


# Import AMReX. Do it here so it is visibile to all targets
find_package(AMReX 19.03.1.348 REQUIRED HINTS ${AMREX_INSTALL_DIR})
echo_amrex_config_options ()
list(APPEND CMAKE_MODULE_PATH ${AMREX_MODULES_PATH})

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

set ( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp )

if (ENABLE_CUDA)
   setup_cuda()
endif ()

#
# Core MFIX library
#
add_subdirectory (src)

# Print out config summary
print_mfix_configuration_summary (${MFIX_LIBNAME})

# We should check what this line does
file (GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)


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

# Configure type-checking target
include ( MFIX_Typecheck )
add_typecheck_target()

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
