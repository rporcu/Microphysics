cmake_minimum_required(VERSION 3.14)

project(MFIX-Exa
   # VERSION <version>
   DESCRIPTION  "A multiphase modeling tool for the exascale"
   HOMEPAGE_URL "https://amrex-codes.github.io/MFIX-Exa/"
   LANGUAGES    C CXX Fortran
   )

# Set the search path for cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/tools/CMake )

# Load utilities and MFIX options module
include( MFIX_Utils )

# Creates `compile_commands.json` in the build build directory
#  `compile_commands.json` contains compiler flags used by plugins like YouCompleteMe
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# Set default build type to Release
if ( NOT CMAKE_BUILD_TYPE )
   message(STATUS "Setting build type to Release as none was specified.")
   set( CMAKE_BUILD_TYPE Release )
else ()
   message(STATUS "Build type set by user to '${CMAKE_BUILD_TYPE}'.")
endif()

# Setup AMReX
include(MFIXSetupAMReX)

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
file (GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90 ${CMAKE_CURRENT_BINARY_DIR}/*.cpp)


add_executable (${MFIX_EXENAME} ${MAIN} ${USR_OVERRIDES} )
target_link_libraries ( ${MFIX_EXENAME} ${MFIX_LIBNAME} )

if ( ENABLE_CUDA)
   set_source_files_properties(${MAIN} PROPERTIES LANGUAGE CUDA )
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
