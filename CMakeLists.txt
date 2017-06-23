cmake_minimum_required (VERSION 3.6)

project (MFIX)

# Define the languages used by the project
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
set ( CMAKE_VERBOSE_MAKEFILE on )

#
# Set the search path for cmake modules
#
set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Tools/)

#
# Initialize all the global variables
#
include ( MFIX_CMakeVariables )

#
# Include modules with custom macros and functions
#
include ( MFIX_Utils )


# Additional build options
include ( MFIX_Options )
include ( MFIX_Config )


include(MFIXTools)
include(AMReXConfig)

check_build_tree_path()
set_defaults()
#set_compiler_flags()

# Config Amrex
if (ENABLE_SUPERBUILD)
   build_amrex()
else (ENABLE_SUPERBUILD)
   find_amrex()
endif (ENABLE_SUPERBUILD)

# Set build flags ( provides MFIX_INCLUDES,
# MFIX_EXT_LINK_STRING  and MFIX_EXT_LINK_DIRS )
set_build_flags()

# Print summary of config settings
summary()

project(MFIXEXA)

# This is the main file
set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp )

include_directories(${MFIX_INCLUDES})
link_directories(${MFIX_EXT_LINK_DIRS})

# Main
add_subdirectory( src )

file(GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
add_executable(mfix ${MAIN} ${USR_OVERRIDES})
target_link_libraries( mfix mfixcore ${MFIX_EXT_LINK_STRING} )

# Tools
set( MFIX_TOOLS_EXE ${PROJECT_BINARY_DIR}/Tools/F90) 
add_subdirectory( Tools/F90 )

# Tests
enable_testing()
add_subdirectory( tests )
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
