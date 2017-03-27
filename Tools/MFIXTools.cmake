###################################################################
# Check if dir or file given by path exists and issue a warning or
#  error if not
##################################################################
function ( check_path  path  message_type )
   if ( EXISTS ${path} )
   else ()
      message(${message_type} ${path} " does not exist!")
   endif ( EXISTS ${path} )
endfunction ()

#################################
# Set defaults options
#################################
macro(set_defaults)

   # Defaults for build type
   if ( NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE "Debug"
	 CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
	 FORCE )
   endif ()
 
   # For compiler options, use lists instead of plain string so duplicates can be removed later on
   # GNU compiler flags 
   set (CXX_FLAGS_RELEASE_GNU "-std=c++11 -g -Wall -O2" CACHE STRING "" FORCE)  
   set (CXX_FLAGS_DEBUG_GNU   "-fmax-errors=2 -std=c++11 -g -Wall -Wextra -O0" CACHE STRING "" FORCE) 
   set (Fortran_FLAGS_DEBUG_GNU "-fmax-errors=2 -g -std=f2008 -fcheck=all -Wall -Wextra -Wimplicit-interface" 
      CACHE STRING "" FORCE)
   set (Fortran_FLAGS_RELEASE_GNU "-g -std=f2008 -fcheck=all -Wall -O2" CACHE STRING "" FORCE)

   # Intel compiler flags
   set (CXX_FLAGS_RELEASE_Intel "" CACHE STRING "" FORCE)  
   set (CXX_FLAGS_DEBUG_Intel   "" CACHE STRING "" FORCE) 
   set (Fortran_FLAGS_DEBUG_Intel "" 
      CACHE STRING "" FORCE)
   set (Fortran_FLAGS_RELEASE_Intel "" CACHE STRING "" FORCE)

   # User-defined flags
   set (CXXFLAGS "" CACHE STRING "User-defined C++ flags" )
   set (FFLAGS   "" CACHE STRING "User-defined Fortran flags" )

   # Set MPI flags (only required for C++)
   set( MPI_INCLUDES "")
   set( MPI_LFLAGS "")
   set( MPI_LIBS "" )
 
   # Set some defaults
   set(BL_SPACEDIM 3 CACHE INT "Dimension of AMReX build")
   set(ENABLE_MPI 0 CACHE INT "Enable build with MPI")
   set(ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
   set(BL_PRECISION "DOUBLE" CACHE INT "Precision of AMReX build")
   set(BL_USE_PARTICLES 1 CACHE INT "Include Particles classes in AMReX build")
   set(ENABLE_PROFILING 0 CACHE INT "Include profiling information in AMReX build")
   set(ENABLE_BACKTRACE 1 CACHE INT "Include backtrace information in AMReX build")
   set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)
   set(ENABLE_SUPERBUILD 1 CACHE INT  "Allow CMake to build and install AMReX")
      
endmacro(set_defaults)

###############################
# Set build flags
###############################
macro( set_build_flags )

   # Includes 
   set( MFIX_EXT_INCLUDE_PATH  ${AMREX_INCLUDES} ${MPI_CXX_INCLUDE_PATH} )
   set( MFIX_INCLUDES  ${MFIX_EXT_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/src
      ${CMAKE_CURRENT_SOURCE_DIR}/src/include )

   # Set MFIX external includes
   set( MFIX_EXT_INCLUDE_PATH  ${AMREX_INCLUDES} ${MPI_CXX_INCLUDE_PATH} )
   set( MFIX_EXT_LINK_STRING   "${AMREX_LIBRARIES} ${MPI_LFLAGS} ${MPI_LIBS}" )
   set( MFIX_EXT_LINK_DIRS     ${AMREX_LIB_DIR} )

   string( STRIP "${MFIX_EXT_LINK_STRING}" MFIX_EXT_LINK_STRING )

endmacro( set_build_flags )


#################################
# Set compiler flags
#################################
macro(set_compiler_flags)

   if (${CMAKE_Fortran_COMPILER_ID} STREQUAL ${CMAKE_CXX_COMPILER_ID} )
      set(COMPILER_ID   ${CMAKE_Fortran_COMPILER_ID})
   else ()
      message(FATAL_ERROR "C++ and Fortran compilers do not have matching IDs.")
   endif()

   # Convert plain string to list
   separate_arguments(CXXFLAGS)
   separate_arguments(FFLAGS)

   # Merge defaults + user-defined flags
   string( TOUPPER ${CMAKE_BUILD_TYPE}  BUILD_TYPE )
   set( CMAKE_CXX_FLAGS_${BUILD_TYPE}      ${CXX_FLAGS_${BUILD_TYPE}_${COMPILER_ID}} ${CXXFLAGS} )
   set( CMAKE_Fortran_FLAGS_${BUILD_TYPE}  ${Fortran_FLAGS_${BUILD_TYPE}_${COMPILER_ID}} ${FFLAGS} )

   # Remove duplicates
   list( REMOVE_DUPLICATES CMAKE_CXX_FLAGS_${BUILD_TYPE} )

   # Convert list to plain string
   string (REPLACE ";" " "  CMAKE_CXX_FLAGS_${BUILD_TYPE} "${CMAKE_CXX_FLAGS_${BUILD_TYPE}}")
   string (REPLACE ";" " "  CMAKE_Fortran_FLAGS_${BUILD_TYPE} "${CMAKE_Fortran_FLAGS_${BUILD_TYPE}}")

endmacro(set_compiler_flags)


#####################################
# Print build config values
#####################################
macro(summary)
   message( STATUS "MFIX build settings summary: ")
   message( STATUS "   Build type        = ${CMAKE_BUILD_TYPE}")
   message( STATUS "   C++ compiler      = ${CMAKE_CXX_COMPILER}")
   message( STATUS "   Fortran compiler  = ${CMAKE_Fortran_COMPILER}")
   message( STATUS "   C++ flags         = ${CMAKE_CXX_FLAGS_${BUILD_TYPE}}")
   message( STATUS "   Fortran flags     = ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}}")
   message( STATUS "   Include paths     = ${MFIX_EXT_INCLUDE_PATH}") 
   message( STATUS "   Link paths        = ${MFIX_EXT_LINK_DIRS}")
   message( STATUS "   Link string       = ${MFIX_EXT_LINK_STRING}")
endmacro(summary)