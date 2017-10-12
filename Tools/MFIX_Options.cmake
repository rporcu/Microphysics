#
# Here we define the default config options   
# that can be overwritten by the user         
# This file provides the followign variables  
#
# DEBUG
# MFIX_BUILD_TYPE
# MFIX_FFLAGS_OVERRIDES
# MFIX_CXXFLAGS_OVERRIDES
# ENABLE_FPE
#

if (DEFINED __MFIX_OPTIONS__)
   return ()
endif ()

# Define the following variable
# so that other included file can check if this file has been
# run already
set (__MFIX_OPTIONS__ "")


#
# Populate the cache and check the value of the user-definable options 
#
option ( DEBUG "Build in debug mode" OFF )

if ( DEBUG )
   set ( CMAKE_BUILD_TYPE "Debug" )
else ()
   set ( CMAKE_BUILD_TYPE "Release" )
endif ()

set ( MFIX_BUILD_TYPE ${CMAKE_BUILD_TYPE} )
string ( TOUPPER ${MFIX_BUILD_TYPE} MFIX_BUILD_TYPE) 


set (MFIX_FFLAGS_OVERRIDES "" CACHE STRING "User-defined Fortran compiler flags" )

set (MFIX_CXXFLAGS_OVERRIDES "" CACHE STRING "User-defined C++ compiler flags" )

option ( ENABLE_FPE "Enable Floating Point Exceptions checks" OFF )






