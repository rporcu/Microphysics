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


#
# Decide whether or not to use PIC 
#
if (ENABLE_PIC)
   set (CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()


# No idea why we need this.
# I think it was required for Franklin build. -- lpritch
# if(PREFER_STATIC_LIBRARIES)
#   # Prefer static libraries, but don't require that everything must be static. 
#   # This appears to be necessary on Franklin at NERSC right now.  --RTM
#   set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
# endif(PREFER_STATIC_LIBRARIES)



# #
# # Detect Fortran name mangling scheme for C/Fortran interface 
# #
# include(FortranCInterface)
# include(${FortranCInterface_BINARY_DIR}/Output.cmake)

# if ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "" )  AND
#       ( FortranCInterface_GLOBAL_CASE STREQUAL "UPPER") )
#    message(STATUS "Fortran name mangling scheme to UPPERCASE \
# (upper case, no append underscore)")
#    add_define ( BL_FORT_USE_UPPERCASE MFIX_DEFINES )
# elseif ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "") AND
#       ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER") )
#    message(STATUS "Fortran name mangling scheme to LOWERCASE \
# (lower case, no append underscore)")
#    add_define ( BL_FORT_USE_LOWERCASE MFIX_DEFINES )
# elseif ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "_" ) AND
#       ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER") )
#    message(STATUS "Fortran name mangling scheme to UNDERSCORE \
# (lower case, append underscore)")
#    add_define ( BL_FORT_USE_UNDERSCORE MFIX_DEFINES )
# else ()
#    message(AUTHOR_WARNING "Fortran to C mangling not backward\
#  compatible with older style BoxLib code") 
# endif ()


# ------------------------------------------------------------- #
#    Set preprocessor flags 
# ------------------------------------------------------------- #

#
# This defines were always on in older version
# of AMReX/CMake. Need more details on these???
#
# add_define (BL_NOLINEVALUES MFIX_DEFINES)
# add_define (BL_PARALLEL_IO MFIX_DEFINES)
# add_define (BL_SPACEDIM=${BL_SPACEDIM} MFIX_DEFINES)
# add_define (BL_${CMAKE_SYSTEM_NAME} MFIX_DEFINES)


# if ( ENABLE_DP )
#    add_define (BL_USE_DOUBLE MFIX_DEFINES)
# else ()
#    add_define (BL_USE_FLOAT MFIX_DEFINES)
# endif ()

# add_define (USE_PARTICLES MFIX_DEFINES ENABLE_PARTICLES)
# add_define (BL_PROFILING MFIX_DEFINES ENABLE_PROFILING) 
# add_define (BL_COMM_PROFILING MFIX_DEFINES ENABLE_COMM_PROFILING)

# if ( ENABLE_PARTICLES AND ( NOT ENABLE_DP_PARTICLES ) ) 
#    add_define ( BL_SINGLE_PRECISION_PARTICLES MFIX_DEFINES )
# endif ()

# if (ENABLE_FORTRAN_MPI AND ENABLE_MPI)
#    add_define ( BL_USE_FORTRAN_MPI=1 MFIX_DEFINES )
# endif ()

# add_define (MG_USE_FBOXIB=1 MFIX_DEFINES ENABLE_MG_FBOXLIB)
# add_define (BL_USE_F_BASELIB=1  MFIX_DEFINES ENABLE_FBASELIB)
# add_define ( BL_USE_MPI MFIX_DEFINES ENABLE_MPI)
# add_define (BL_USE_OMP MFIX_DEFINES ENABLE_OMP)

# #
# # Add all preprocessor definitions to compile string
# # 
# add_definitions ( ${MFIX_DEFINES} )

# ------------------------------------------------------------- #
#    Setup third party packages 
# ------------------------------------------------------------- #


if ( NOT ENABLE_SUPERBUILD )
   find_amrex ()
   list (APPEND MFIX_EXTRA_Fortran_INCLUDE_PATH "${AMREX_INCLUDES}")
   list (APPEND MFIX_EXTRA_C_INCLUDE_PATH "${AMREX_INCLUDES}")
   list (APPEND MFIX_EXTRA_CXX_INCLUDE_PATH "${AMREX_INCLUDES}")   
   # message (  "CCSE_INCLUDE_DIR = ${AMREX_INCLUDES}" )
   # message (  "CCSE_LIBRARIES   = ${AMREX_LIBRARIES}" )
   # message (  "CCSE_LIBRARY_DIR = ${AMREX_LIB_DIR}"  )
endif()


if (ENABLE_MPI)
   find_package (MPI REQUIRED)
   # Includes
   list (APPEND MFIX_EXTRA_Fortran_INCLUDE_PATH "${MPI_Fortran_INCLUDE_PATH}")
   list (APPEND MFIX_EXTRA_C_INCLUDE_PATH "${MPI_C_INCLUDE_PATH}")
   list (APPEND MFIX_EXTRA_CXX_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}")
   # Compile flags
   append ( MPI_Fortran_COMPILE_FLAGS MFIX_EXTRA_Fortran_FLAGS ) 
   append ( MPI_C_COMPILE_FLAGS MFIX_EXTRA_C_FLAGS )
   append ( MPI_CXX_COMPILE_FLAGS MFIX_EXTRA_CXX_FLAGS )
   # Libraries
   list (APPEND MFIX_EXTRA_Fortran_LIBRARIES "${MPI_Fortran_LIBRARIES}")
   list (APPEND MFIX_EXTRA_C_LIBRARIES "${MPI_C_LIBRARIES}")
   list (APPEND MFIX_EXTRA_CXX_LIBRARIES "${MPI_CXX_LIBRARIES}")
   # Link flags
   list (APPEND MFIX_EXTRA_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
   list (APPEND MFIX_EXTRA_C_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
   list (APPEND MFIX_EXTRA_CXX_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif ()

if (ENABLE_OMP)
   find_package (OpenMP REQUIRED)
   # Compile flags
   append ( OpenMP_Fortran_FLAGS MFIX_EXTRA_Fortran_FLAGS ) 
   append ( OpenMP_C_FLAGS MFIX_EXTRA_C_FLAGS )
   append ( OpenMP_CXX_FLAGS MFIX_EXTRA_CXX_FLAGS )
endif()




# ------------------------------------------------------------- #
#    Setup compiler flags 
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

append ( MFIX_EXTRA_Fortran_FLAGS MFIX_Fortran_FLAGS )
append ( MFIX_EXTRA_CXX_FLAGS MFIX_CXX_FLAGS )

# Add required flags
append ( MFIX_${FC_ID}_FFLAGS_REQUIRED MFIX_Fortran_FLAGS )
append ( MFIX_${CXX_ID}_CXXLAGS_REQUIRED MFIX_CXX_FLAGS )

# Set CMake compiler flags
set ( CMAKE_Fortran_FLAGS_${MFIX_BUILD_TYPE} "${MFIX_Fortran_FLAGS}" ) 
set ( CMAKE_CXX_FLAGS_${MFIX_BUILD_TYPE} "${MFIX_CXX_FLAGS}" )


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
message( STATUS "   C++ include paths     = ${MFIX_EXTRA_CXX_INCLUDE_PATH}") 
message( STATUS "   Fortran include paths = ${MFIX_EXTRA_Fortran_INCLUDE_PATH}")
message( STATUS "   C++ external libs     = ${MFIX_EXTRA_CXX_LIBRARIES}") 
message( STATUS "   Fortran external libs = ${MFIX_EXTRA_Fortran_LIBRARIES}")
message( STATUS "   C++ link flags        = ${MFIX_EXTRA_CXX_LINK_FLAGS}") 
message( STATUS "   Fortran link flags    = ${MFIX_EXTRA_LINK_FLAGS}")
