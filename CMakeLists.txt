cmake_minimum_required (VERSION 3.5)

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
#
if ( SUPERBUILD )
   include ( MFIX_Superbuild )
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
find_package (AMReX REQUIRED HINTS ${AMREX_INSTALL_DIR})
echo_amrex_config_options ()

if ( NOT ( AMREX_ENABLE_PARTICLES AND AMREX_ENABLE_AMRDATA
	 AND AMREX_ENABLE_EB) )
   message ( FATAL_ERROR "AMReX installation missing required components!")
endif ()

if ( NOT ( "${MFIX_BUILD_TYPE}" STREQUAL "${AMREX_BUILD_TYPE}" ) )
   message (WARNING "MFIX build type (${MFIX_BUILD_TYPE}) type does not match AMReX build type (${AMREX_BUILD_TYPE})")
endif ()

#
# This is the main file
#
if ( ENABLE_PROJCC )
  set ( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src_cc/main.cpp )
else ()
  set ( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp )
endif ()

#
# Core MFIX library
# For the time being, we keep to versions: one with sources in src/
# ( SIMPLE + Exact Projection ), one with sources in both src_cc/ and src/
# ( Approximate Projection only ).
#
if (ENABLE_PROJCC)
  message ( WARNING "Building MFIX in experimental APPROXIMATE PROJECTION mode (ENABLE_PROJCC=on)" )
  add_subdirectory (src_cc)
  set (ENABLE_PTESTS ON)
  set (ENABLE_STESTS OFF)
else ()
  add_subdirectory (src)
endif ()

# Print out config summary
print_mfix_configuration_summary (${MFIX_LIBNAME})

# We should check what this line does
file (GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
add_executable (${MFIX_EXENAME} ${MAIN} ${USR_OVERRIDES} )
target_link_libraries ( ${MFIX_EXENAME} ${MFIX_LIBNAME} )

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
