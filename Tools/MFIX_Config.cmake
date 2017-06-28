###############################################

# Here we configure the build                 #

###############################################

#
#  Check if AMReX_Options.cmake has been already processed
#
if ( NOT MFIX_OPTIONS_SET )
   message ( FATAL_ERROR "MFIX_Options.cmake must be\
included before MFIX_Config.cmake" )
endif ()

# ------------------------------------------------------------- #
#  Setup core compiler flags 
# ------------------------------------------------------------- #
if ( MFIX_FFLAGS_OVERRIDES )
   set ( MFIX_Fortran_FLAGS ${MFIX_FFLAGS_OVERRIDES} )
else ()
   append ( MFIX_${FC_ID}_FFLAGS_${MFIX_BUILD_TYPE}
      MFIX_Fortran_FLAGS )
endif ()

if ( MFIX_CXXLAGS_OVERRIDES )
   set ( MFIX_CXX_FLAGS ${MFIX_CXXFLAGS_OVERRIDES} )
else ()
   append ( MFIX_${CXX_ID}_CXXFLAGS_${MFIX_BUILD_TYPE}
      MFIX_CXX_FLAGS )
endif ()


# ------------------------------------------------------------- #
#  Setup amrex and its dependencies 
# ------------------------------------------------------------- #
if (AMREX_INSTALL_DIR) # No superbuild
   check_path ( ${AMREX_INSTALL_DIR} FATAL_ERROR )
   check_path ( ${AMREX_INSTALL_DIR}/cmake FATAL_ERROR )
   set (ENABLE_SUPERBUILD 0)
   set (AMREX_INSTALL_PATH ${AMREX_INSTALL_DIR} )
   find_package (AMReX CONFIG REQUIRED HINTS ${AMREX_INSTALL_PATH}/cmake )   
   # Check and print (if not superbuild) amrex options
   if (NOT ENABLE_PARTICLES)
      message ( FATAL_ERROR "AMReX must be configured with -DENABLE_PARTICLES=1" )
   endif ()

   #
   # Echo amrex config options
   #
   message (STATUS "AMReX configuration options: ")
   print_option (AMREX_BUILD_TYPE ${AMREX_BUILD_TYPE})
   print_option (FORTRAN_ENABLE_MPI ${FORTRAN_ENABLE_MPI})
   print_option (ENABLE_FBASELIB ${ENABLE_FBASELIB})
   print_option (ENABLE_PIC ${ENABLE_PIC})
   print_option (BL_SPACEDIM ${BL_SPACEDIM})
   print_option (ENABLE_MPI ${ENABLE_MPI})
   print_option (ENABLE_OMP ${ENABLE_OMP})
   print_option (ENABLE_DP ${ENABLE_DP})
   print_option (ENABLE_PARTICLES ${ENABLE_PARTICLES})
   print_option (ENABLE_DP_PARTICLES ${ENABLE_DP_PARTICLES})
   print_option (ENABLE_PROFILING ${ENABLE_PROFILING})
   print_option (ENABLE_TINY_PROFILING ${ENABLE_TINY_PROFILING})
   print_option (ENABLE_BACKTRACE ${ENABLE_BACKTRACE})
   print_option (ENABLE_MG_BOXLIB ${ENABLE_MG_BOXLIB})
else ()
   set (ENABLE_SUPERBUILD 1)
   include (MFIX_ConfigSuperbuild)
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
add_definitions ( ${MFIX_DEFINES} )

# Add amrex library path and extra link line
append_to_link_line ( AMREX_LIBRARIES MFIX_EXTRA_LINK_LINE )
append_to_link_line ( AMREX_EXTRA_CXX_LINK_LINE MFIX_EXTRA_LINK_LINE )
list (APPEND MFIX_EXTRA_LIBRARIES_PATH "${AMREX_LIBRARY_DIR}")

# Add amrex extra compiler flags
append ( AMREX_EXTRA_Fortran_FLAGS MFIX_EXTRA_Fortran_FLAGS )
append ( AMREX_EXTRA_C_FLAGS MFIX_EXTRA_C_FLAGS )
append ( AMREX_EXTRA_CXX_FLAGS MFIX_EXTRA_CXX_FLAGS )

# ------------------------------------------------------------- #
#  Finalize configuration 
# ------------------------------------------------------------- #
append ( MFIX_EXTRA_Fortran_FLAGS MFIX_Fortran_FLAGS )
append ( MFIX_EXTRA_CXX_FLAGS MFIX_CXX_FLAGS )

# Add required flags
append ( MFIX_${FC_ID}_FFLAGS_REQUIRED MFIX_Fortran_FLAGS )
append ( MFIX_${CXX_ID}_CXXFLAGS_REQUIRED MFIX_CXX_FLAGS )

# Set CMake compiler flags
set ( CMAKE_Fortran_FLAGS_${MFIX_BUILD_TYPE} "${MFIX_Fortran_FLAGS}" ) 
set ( CMAKE_CXX_FLAGS_${MFIX_BUILD_TYPE} "${MFIX_CXX_FLAGS}" )

# Set total extra  includes path
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
message( STATUS "   C++ flags             = ${CMAKE_CXX_FLAGS_${MFIX_BUILD_TYPE}}")
# message( STATUS "   Fortran flags         = ${CMAKE_Fortran_FLAGS_${MFIX_BUILD_TYPE}}")
# message( STATUS "   C++ include paths     = ${MFIX_EXTRA_CXX_INCLUDE_PATH}") 
# message( STATUS "   Fortran include paths = ${MFIX_EXTRA_Fortran_INCLUDE_PATH}")
# message( STATUS "   C++ external libs     = ${MFIX_EXTRA_CXX_LIBRARIES}") 
# message( STATUS "   Fortran external libs = ${MFIX_EXTRA_Fortran_LIBRARIES}")
# message( STATUS "   C++ link flags        = ${MFIX_EXTRA_CXX_LINK_FLAGS}") 
message( STATUS "   Fortran link flags    = ${MFIX_EXTRA_LINK_FLAGS}")
message( STATUS "   MFIX extra link line  = ${MFIX_EXTRA_LINK_LINE}")
message( STATUS "   MFIX extra includes   = ${MFIX_EXTRA_INCLUDE_PATH}")
