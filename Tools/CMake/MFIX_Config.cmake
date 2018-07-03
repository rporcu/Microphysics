###############################################

# Here we configure the build                 #

###############################################
if (DEFINED __MFIX_CONFIG__)
   return ()
endif ()
set ( __MFIX_CONFIG__ "" )

#
# Check if MFIX_Options.cmake and MFIX_CMakeVariables.cmake
# have been already processed
#
if (NOT (DEFINED __MFIX_OPTIONS__) )
   message ( FATAL_ERROR "MFIX_Options.cmake must be \
included before MFIX_Config.cmake" )
endif ()

# if (NOT (DEFINED __MFIX_CMAKEVARIABLES__) )
#    message ( FATAL_ERROR "MFIX_CMakeVariables.cmake must be \
# included before MFIX_Config.cmake" )
# endif ()

set (AMREX_INSTALL_PATH ${AMREX_INSTALL_DIR} )



# Set amrex typecheker path
set ( AMREX_TYPECHECKER  "${AMREX_INSTALL_PATH}/Tools/typechecker/typechecker.py" )

# 
# Set CMake compiler flags
# CMAKE_<LANG>_FLAGS_<BUILD_TYPE> is not used: this is
# NOT the standard CMake way of doing things.
#
set ( CMAKE_Fortran_FLAGS_${MFIX_BUILD_TYPE} "" )
set ( CMAKE_Fortran_FLAGS "${MFIX_Fortran_FLAGS}" ) 
set ( CMAKE_CXX_FLAGS_${MFIX_BUILD_TYPE} "" )
set ( CMAKE_CXX_FLAGS  "${MFIX_CXX_FLAGS}" )


#
# Config summary
#
message( STATUS "MFIX configuration summary: ")
message( STATUS "   Build type            = ${CMAKE_BUILD_TYPE}")
message( STATUS "   Preprocessor flags    = ${MFIX_DEFINES}")
message( STATUS "   C++ compiler          = ${CMAKE_CXX_COMPILER}")
message( STATUS "   Fortran compiler      = ${CMAKE_Fortran_COMPILER}")
message( STATUS "   C++ flags             = ${MFIX_CXX_FLAGS}")
message( STATUS "   Fortran flags         = ${MFIX_Fortran_FLAGS}")
message( STATUS "   MFIX extra link line  = ${MFIX_EXTRA_LINK_LINE}")
message( STATUS "   MFIX extra includes   = ${MFIX_EXTRA_INCLUDE_PATH}")
