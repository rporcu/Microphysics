cmake_minimum_required (VERSION 3.5)

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
set ( CMAKE_VERBOSE_MAKEFILE on )

#
# Set the search path for cmake modules
#
set ( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Tools/ )


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
