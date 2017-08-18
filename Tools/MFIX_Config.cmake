###############################################

# Here we configure the build                 #

###############################################
if (DEFINED __MFIX_CONFIG__)
   return ()
endif ()
set ( __MFIX_CONFIG__ "" )

#
#  Check if MFIX_Options.cmake has been already processed
#
if (NOT (DEFINED __MFIX_OPTIONS__) )
   message ( FATAL_ERROR "MFIX_Options.cmake must be \
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

if ( MFIX_CXXFLAGS_OVERRIDES )
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

   # Check and print (if not superbuild) amrex options
   if (NOT ENABLE_FBASELIB)
      message ( FATAL_ERROR "AMReX must be configured with -DENABLE_FBASELIB=1" )
   endif ()

  
   # Check and print (if not superbuild) amrex options
   if ( NOT (${AMREX_DIM} EQUAL 3) )
      message ( FATAL_ERROR "AMReX must be configured with -DDIM=3" )
   endif ()

   
   #
   # Echo amrex config options
   #
   echo_amrex_config_options ()
   # message (STATUS "AMReX configuration options: ")
   # print_option (AMREX_DEBUG ${AMREX_DEBUG})
   # print_option (AMRED_DIM ${AMREX_DIM})
   # print_option (ENABLE_PIC ${ENABLE_PIC})
   # print_option (ENABLE_MPI ${ENABLE_MPI})
   # print_option (ENABLE_OMP ${ENABLE_OMP})
   # print_option (ENABLE_DP ${ENABLE_DP})
   # print_option (ENABLE_FORTRAN_INTERFACES ${ENABLE_FORTRAN_INTERFACES})
   # print_option (ENABLE_LINEAR_SOLVERS ${ENABLE_LINEAR_SOLVERS})
   # print_option (ENABLE_FBASELIB ${ENABLE_FBASELIB})
   # print_option (ENABLE_AMRDATA ${ENABLE_AMRDATA})
   # print_option (ENABLE_PARTICLES ${ENABLE_PARTICLES})
   # print_option (ENABLE_DP_PARTICLES ${ENABLE_DP_PARTICLES})
   # print_option (ENABLE_PROFILING ${ENABLE_PROFILING})
   # print_option (ENABLE_TINY_PROFILING  ${ENABLE_TINY_PROFILING})
   # print_option (ENABLE_TRACE_PROFILING ${ENABLE_TRACE_PROFILING})
   # print_option (ENABLE_COMM_PROFILING  ${ENABLE_COMM_PROFILING})   
   # print_option (ENABLE_BACKTRACE ${ENABLE_BACKTRACE})

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
set (AMREX_LIB_FULLPATH ${AMREX_LIBRARY_DIR}/libamrex.a)
append_to_link_line ( AMREX_LIB_FULLPATH MFIX_EXTRA_LINK_LINE )
append_to_link_line ( AMREX_EXTRA_CXX_LINK_LINE MFIX_EXTRA_LINK_LINE )
list (APPEND MFIX_EXTRA_LIBRARIES_PATH "${AMREX_LIBRARY_DIR}")

# Add amrex extra compiler flags
append ( AMREX_EXTRA_Fortran_FLAGS MFIX_EXTRA_Fortran_FLAGS )
append ( AMREX_EXTRA_C_FLAGS MFIX_EXTRA_C_FLAGS )
append ( AMREX_EXTRA_CXX_FLAGS MFIX_EXTRA_CXX_FLAGS )

# Set amrex typecheker path
set ( AMREX_TYPECHECKER  "${AMREX_INSTALL_PATH}/Tools/typechecker/typechecker.py" )


# ------------------------------------------------------------- #
#  Finalize configuration 
# ------------------------------------------------------------- #
append ( MFIX_EXTRA_Fortran_FLAGS MFIX_Fortran_FLAGS )
append ( MFIX_EXTRA_CXX_FLAGS MFIX_CXX_FLAGS )

# Add required flags
append ( MFIX_${FC_ID}_FFLAGS_REQUIRED MFIX_Fortran_FLAGS )
append ( MFIX_${CXX_ID}_CXXFLAGS_REQUIRED MFIX_CXX_FLAGS )

# Add FPE flags if required 
if (ENABLE_FPE)
   append ( MFIX_${FC_ID}_FFLAGS_FPE MFIX_Fortran_FLAGS )
   append ( MFIX_${CXX_ID}_CXXFLAGS_FPE MFIX_CXX_FLAGS )
endif ()


# Set CMake compiler flags
set ( CMAKE_Fortran_FLAGS_${MFIX_BUILD_TYPE}
   "${MFIX_Fortran_FLAGS} ${AMREX_Fortran_DEFINITIONS}" ) 
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
message( STATUS "   Fortran flags         = ${CMAKE_Fortran_FLAGS_${MFIX_BUILD_TYPE}}")
# message( STATUS "   C++ include paths     = ${MFIX_EXTRA_CXX_INCLUDE_PATH}") 
# message( STATUS "   Fortran include paths = ${MFIX_EXTRA_Fortran_INCLUDE_PATH}")
# message( STATUS "   C++ external libs     = ${MFIX_EXTRA_CXX_LIBRARIES}") 
# message( STATUS "   Fortran external libs = ${MFIX_EXTRA_Fortran_LIBRARIES}")
# message( STATUS "   C++ link flags        = ${MFIX_EXTRA_CXX_LINK_FLAGS}") 
#message( STATUS "   Fortran link flags    = ${MFIX_EXTRA_LINK_FLAGS}")
message( STATUS "   MFIX extra link line  = ${MFIX_EXTRA_LINK_LINE}")
message( STATUS "   MFIX extra includes   = ${MFIX_EXTRA_INCLUDE_PATH}")
