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

if (NOT (DEFINED __MFIX_CMAKEVARIABLES__) )
   message ( FATAL_ERROR "MFIX_CMakeVariables.cmake must be \
included before MFIX_Config.cmake" )
endif ()


# 
#  Setup amrex and its dependencies 
# 
check_path ( ${AMREX_INSTALL_DIR} FATAL_ERROR )
check_path ( ${AMREX_INSTALL_DIR}/cmake FATAL_ERROR )

set (AMREX_INSTALL_PATH ${AMREX_INSTALL_DIR} )


find_package (AMReX CONFIG REQUIRED HINTS ${AMREX_INSTALL_DIR}/cmake )   

# Check and print (if not superbuild) amrex options
if (NOT AMREX_ENABLE_PARTICLES)
   message ( FATAL_ERROR "AMReX must be configured with -DENABLE_PARTICLES=1" )
endif ()

if (NOT AMREX_ENABLE_EB)
   message ( FATAL_ERROR "AMReX must be configured with -DENABLE_EB=1" )
endif ()

# Check and print (if not superbuild) amrex options
if (NOT AMREX_ENABLE_AMRDATA)
   message ( FATAL_ERROR "AMReX must be configured with -DENABLE_AMRDATA=1" )
endif ()

# Check and print (if not superbuild) amrex options
if ( NOT (${AMREX_DIM} EQUAL 3) )
   message ( FATAL_ERROR "AMReX must be configured with -DDIM=3" )
endif ()


#
# Echo amrex config options ( this function is provided by the config file )
#
echo_amrex_config_options ()

#
# Load AMReX compile flags. Set MFIX flags to AMReX flags only if
# user did not set CMAKE_<LANG>_FLAGS 
#
if ( CMAKE_Fortran_FLAGS )
   set ( MFIX_Fortran_FLAGS
      "${CMAKE_Fortran_FLAGS} ${AMREX_EXTRA_Fortran_FLAGS} ${AMREX_Fortran_DEFINES}" )
else ()
   set ( MFIX_Fortran_FLAGS "${AMREX_Fortran_FLAGS}" )
endif ()

if ( CMAKE_CXX_FLAGS )
   set ( MFIX_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AMREX_EXTRA_CXX_FLAGS}" )
else ()
   set ( MFIX_CXX_FLAGS "${AMREX_CXX_FLAGS}" )
endif ()


# Need to add both amrex include directory + amrex external includes
list (APPEND MFIX_EXTRA_Fortran_INCLUDE_PATH  "${AMREX_EXTRA_Fortran_INCLUDE_PATH}")
list (APPEND MFIX_EXTRA_Fortran_INCLUDE_PATH  "${AMREX_INCLUDE_DIR}")
list (APPEND MFIX_EXTRA_C_INCLUDE_PATH "${AMREX_EXTRA_C_INCLUDE_PATH}")
list (APPEND MFIX_EXTRA_C_INCLUDE_PATH  "${AMREX_INCLUDE_DIR}")
list (APPEND MFIX_EXTRA_CXX_INCLUDE_PATH "${AMREX_EXTRA_CXX_INCLUDE_PATH}")
list (APPEND MFIX_EXTRA_CXX_INCLUDE_PATH  "${AMREX_INCLUDE_DIR}")

# Add amrex definitions (needed for amrex headers)
set ( MFIX_DEFINES  ${AMREX_DEFINES} )

# For external AMReX installs, check consistency of
# build type. If AMReX has been built in Release mode, but
# MFIX hasn't, remove NDEBUG defines from MFIX_DEFINES
if ( ( NOT AMREX_DEBUG ) AND DEBUG )
   list ( REMOVE_ITEM MFIX_DEFINES "-DNDEBUG")
endif ()
add_definitions ( ${MFIX_DEFINES} )

# Add amrex library path and extra link line
set (AMREX_LIB_FULLPATH ${AMREX_LIBRARY_DIR}/libamrex.a)
append_to_link_line ( AMREX_LIB_FULLPATH MFIX_EXTRA_LINK_LINE )
append_to_link_line ( AMREX_EXTRA_CXX_LINK_LINE MFIX_EXTRA_LINK_LINE )
append_to_link_line ( AMREX_EXTRA_Fortran_LINK_LINE MFIX_EXTRA_LINK_LINE )
list ( REMOVE_DUPLICATES MFIX_EXTRA_LINK_LINE )
list ( APPEND MFIX_EXTRA_LIBRARIES_PATH "${AMREX_LIBRARY_DIR}")


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

# Set total extra includes path
set ( MFIX_EXTRA_INCLUDE_PATH  ${MFIX_EXTRA_Fortran_INCLUDE_PATH}
   ${MFIX_EXTRA_C_INCLUDE_PATH} ${MFIX_EXTRA_CXX_INCLUDE_PATH} )

# Remove duplicates from lists
if ( MFIX_EXTRA_INCLUDE_PATH )
   list ( REMOVE_DUPLICATES MFIX_EXTRA_INCLUDE_PATH )
endif ()

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


# 
# Here we check if we can enable testing (i.e. check if compiler supports
# OpenMP > = 4.0 which is used to built plt_compare_diff_grid when ENABLE_OMP=on )
#
set ( MFIX_ENABLE_CTEST 1 )


# Just checking GNU and Intel: I doubt that ctest will be run on Cray machines
# And who uses PGI anyway :-D
if ( ENABLE_OMP )
   if ( ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" ) AND
	( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5.0 ) )
     set ( MFIX_ENABLE_CTEST 0 )
  elseif ( ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" ) AND
	( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 16.0 ) )
     set ( MFIX_ENABLE_CTEST 0 )
  endif ()
endif ()

if (NOT MFIX_ENABLE_CTEST)
   message (WARNING "The compiler does not support OpenMP >= 4.0: disabling testing suite")
endif ()
