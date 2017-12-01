cmake_minimum_required (VERSION 3.5)

# Creates `compile_commands.json` in the build build directory
#  `compile_commands.json` contains compiler flags used by plugins like YouCompleteMe
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )

IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()


#
# Set the search path for cmake modules
#
set ( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Tools/ )

#
# Here we set the only options that decide if superbuild
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

#
# Define the languages used by the project
#
enable_language (C)
enable_language (CXX)
enable_language (Fortran)

#
# Require C++11 standard
#
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON) 
set (CMAKE_CXX_EXTENSIONS OFF)

#
# Always verbose output during make
# 
# set ( CMAKE_VERBOSE_MAKEFILE on )


#
# Include module with custom macros and functions
#
include ( MFIX_Utils )

get_git_info ( var1 var2 )

#
# Stop if building in source
# 
check_build_tree_path ()

#
# Initialize all the global variables
#
include ( MFIX_CMakeVariables )

# 
# Load MFIX options
# 
include ( MFIX_Options )

#
# Config MFIX build
# 
include ( MFIX_Config )


# This is the main file
set ( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp )

include_directories (${MFIX_INCLUDE_PATH} ${MFIX_EXTRA_INCLUDE_PATH})
link_directories (${MFIX_EXTRA_LIBRARIES_PATH})

# Main
add_subdirectory (src)

file (GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
add_executable (mfix ${MAIN} ${USR_OVERRIDES} )
target_link_libraries ( mfix mfixcore ${MFIX_EXTRA_LINK_LINE} ) 

# Add rule to generate TAGS
add_custom_target ( tags
   COMMAND ctags-exuberant -R    --fortran-kinds=+i . 
   COMMAND ctags-exuberant -R -e --fortran-kinds=+i . 
   COMMENT "Generating tags files"
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )


# Tools
set (MFIX_TOOLS_EXE ${PROJECT_BINARY_DIR}/Tools/F90) 
add_subdirectory (Tools/F90)

# Tests
enable_testing ()
add_subdirectory ( tests )
#add_subdirectory( benchmarks )

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
