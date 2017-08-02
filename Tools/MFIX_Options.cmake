###############################################

# Here we define the default config options   #
# that can be overwritten by the user         #

###############################################

if (DEFINED __MFIX_OPTIONS__)
   return ()
endif ()
# Define the following variable
# so that other included file can check if this file has been
# run already
set (__MFIX_OPTIONS__ "")


#
# Check weather the MFIX_CMakeVariables.cmake
# has been loaded; abort if not
#
if (NOT (DEFINED __MFIX_CMAKEVARIABLES__) )
   message ( FATAL_ERROR "MFIX_Options.cmake must be included \
after including MFIX_CMakeVariables.cmake" )
endif ()


#
# Define a macro to check the value of the inputs integer options 
# 
macro (check_option_value NAME VALUE RANGE_STARTS RANGE_ENDS )
   
   if ( ( ${VALUE} GREATER ${RANGE_ENDS} ) OR ( ${VALUE} LESS ${RANGE_STARTS} ) )
      message ( FATAL_ERROR "Variable ${NAME} has value ${VALUE}. \
Allowed range is [${RANGE_STARTS}:${RANGE_ENDS}]." )
   endif ()

   message ( STATUS "   ${NAME} = ${VALUE} (INT: ${RANGE_STARTS},${RANGE_ENDS})" )
endmacro ()


#
# Populate the cache and check the value of the user-definable options 
#
message (STATUS "MFIX configuration options: ")

if ( NOT CMAKE_BUILD_TYPE )
   # Default to debug if no other build type specified
   set ( CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE )
endif ()

set ( MFIX_BUILD_TYPE ${CMAKE_BUILD_TYPE} )

# Need to be uppercase so it can be used to refer to CMAKE variable names
string ( TOUPPER ${MFIX_BUILD_TYPE} MFIX_BUILD_TYPE) 

message ( STATUS "   CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE} (STRING:\
 Debug|Release)" )

set (AMREX_INSTALL_DIR "" CACHE PATH "Path to installation directory (leave empty for superbuild)")


set (MFIX_FFLAGS_OVERRIDES "" CACHE STRING "User-defined Fortran compiler flags" )
message (STATUS "   MFIX_FFLAGS_OVERRIDES = ${MFIX_FFLAGS_OVERRIDES}" )
set (MFIX_CXXFLAGS_OVERRIDES "" CACHE STRING "User-defined C++ compiler flags" )
message (STATUS "   MFIX_CXXFLAGS_OVERRIDES = ${MFIX_CXXFLAGS_OVERRIDES}" )

set (ENABLE_FPE 0 CACHE INT "Enable Floating Point Exceptions checks")
check_option_value ( "ENABLE_FPE" ${ENABLE_FPE} 0 1 )


#
# AMReX options are needed only if AMReX is built as part of MFIX (superbuild)
# 
if (NOT AMREX_INSTALL_DIR)  

   message (STATUS "Configuring AMReX with the following options: ")
   
   set (AMREX_ENABLE_PIC 0 CACHE INT
      "Compile with position-independent code enabled")
   check_option_value ( "ENABLE_PIC" ${AMREX_ENABLE_PIC} 0 1 )

   set (AMREX_ENABLE_MPI 1 CACHE INT "Enable MPI in AMReX build")
   check_option_value ( "ENABLE_MPI" ${AMREX_ENABLE_MPI} 0 1 )

   set (AMREX_ENABLE_OMP 0 CACHE INT "Enable OpenMP in AMReX build ")
   check_option_value ( "ENABLE_OMP" ${AMREX_ENABLE_OMP} 0 1 )

   set (AMREX_ENABLE_DP 1 CACHE INT "Enable double precision in AMReX build")
   check_option_value ( "ENABLE_DP" ${AMREX_ENABLE_DP} 0 1 )

   set (AMREX_ENABLE_DP_PARTICLES 1 CACHE INT "Enable double-precision for particles data\
in AMReX build") 
   check_option_value ( "ENABLE_DP_PARTICLES" ${AMREX_ENABLE_DP_PARTICLES} 0 1 )

   set (AMREX_ENABLE_PROFILING 0 CACHE INT "Include profiling information in AMReX build")
   check_option_value ( "ENABLE_PROFILING" ${AMREX_ENABLE_PROFILING} 0 1 )

   set (AMREX_ENABLE_TINY_PROFILING 0 CACHE INT "Include 'tiny'-profiling information in AMReX build")
   check_option_value ( "ENABLE_TINY_PROFILING" ${AMREX_ENABLE_TINY_PROFILING} 0 1 )

   set (AMREX_ENABLE_BACKTRACE 1 CACHE INT "Include backtrace information in AMReX build")
   check_option_value ( "ENABLE_BACKTRACE" ${AMREX_ENABLE_BACKTRACE} 0 1 )

   set (AMREX_ENABLE_ASSERTIONS 0 CACHE INT "Include assertions in AMReX build")
   check_option_value ( "ENABLE_ASSERTIONS" ${AMREX_ENABLE_ASSERTIONS} 0 1 )
   
   set (AMREX_GIT_COMMIT "" CACHE STRING "AMReX git commit to use in superbuild")
   
endif ()




