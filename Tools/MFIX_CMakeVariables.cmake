###############################################

# Here we define all the global variable that #
# will be used in the CMake files

###############################################

#
# Check if project() has been called; abort if not
# 
if ( NOT PROJECT_NAME )
   message ( FATAL_ERROR "MFIX_CMakeVariables.cmake cannot be included\
before calling project()" )
endif ()

# 
# Check if this file has been loaded already
#
if (DEFINED __MFIX_CMAKEVARIABLES__)
   return ()
endif ()

set (__MFIX_CMAKEVARIABLES__ "")

# Set paths for build system
set ( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)

# The type of build ( will need to be uppercase )
set ( MFIX_BUILD_TYPE )

# Variable to hold amrex install path
set ( AMREX_INSTALL_PATH )

# Flags to accumulate preprocessor directives
set ( MFIX_DEFINES ) 

# Compiler flags
set ( MFIX_Fortran_FLAGS )
set ( MFIX_C_FLAGS )
set ( MFIX_CXX_FLAGS )

#
# Compile- and link-time variables 
#
set (MFIX_EXTRA_C_INCLUDE_PATH)
set (MFIX_EXTRA_CXX_INCLUDE_PATH)
set (MFIX_EXTRA_Fortran_INCLUDE_PATH)
set (MFIX_EXTRA_C_FLAGS)
set (MFIX_EXTRA_CXX_FLAGS)
set (MFIX_EXTRA_Fortran_FLAGS)
set (MFIX_EXTRA_C_LIBRARIES)
set (MFIX_EXTRA_CXX_LIBRARIES)
set (MFIX_EXTRA_Fortran_LIBRARIES)
set (MFIX_EXTRA_C_LIBS_DIR)
set (MFIX_EXTRA_CXX_LIBS_DIR)
set (MFIX_EXTRA_C_LINK_FLAGS)
set (MFIX_EXTRA_CXX_LINK_FLAGS)
set (MFIX_EXTRA_Fortran_LINK_FLAGS)
set (MFIX_EXTRA_C_LINK_LINE)
set (MFIX_EXTRA_CXX_LINK_LINE)
set (MFIX_EXTRA_Fortran_LINK_LINE)
set (MFIX_EXTRA_LINK_LINE)
set (MFIX_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/src
   ${PROJECT_SOURCE_DIR}/src/include)
set (MFIX_EXTRA_INCLUDE_PATH)
set (MFIX_EXTRA_LIBRARIES_PATH)

# MFIX git variables
set (MFIX_GIT_COMMIT)

