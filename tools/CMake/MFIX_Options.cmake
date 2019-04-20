#
# Here we define the default config options
# that can be overwritten by the user
# This file provides the following variables
#
# DEBUG
# MFIX_FFLAGS_OVERRIDES
# MFIX_CXXFLAGS_OVERRIDES
# ENABLE_FPE
#

if(DEFINED __MFIX_OPTIONS__)
  return()
endif()

# Define the following variable
# so that other included file can check if this file has been
# run already
set(__MFIX_OPTIONS__ "")

# Creates `compile_commands.json` in the build build directory
#  `compile_commands.json` contains compiler flags used by plugins like YouCompleteMe
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

#
# Populate the cache and check the value of the user-definable options
#
option( DEBUG "Build in debug mode" OFF )

if( DEBUG )
  set( CMAKE_BUILD_TYPE "Debug" )
else()
  set( CMAKE_BUILD_TYPE "Release" )
endif()

option( ENABLE_OMP "Enable OpenMP" NO )

option( ENABLE_MPI "Enable MPI" YES )

option( ENABLE_HYPRE "Enable HYPRE" NO )

option( ENABLE_CUDA "Enable CUDA" NO )

option( ENABLE_FPE "Enable Floating Point Exceptions checks" NO )

option( SUBMOD "Use submodules with Superbuild" NO )
