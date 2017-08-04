#
# This file configures the superbuild 
# 
if (AMREX_GIT_COMMIT)
   set (AMREX_GIT_TAG ${AMREX_GIT_COMMIT})
else ()
   if (MFIX_GIT_COMMIT MATCHES "master")
      set (AMREX_GIT_TAG ${AMREX_GIT_COMMIT_MASTER})
   else ()
      set (AMREX_GIT_TAG ${AMREX_GIT_COMMIT_DEVELOP})
   endif()
endif ()

message (STATUS "SUPERBUILD mode: AMReX will be built as part of MFIX")
message (STATUS "AMReX commit: ${AMREX_GIT_TAG}")

# Include cmake config files to build external projects
include(ExternalProject)


set (AMREX_INSTALL_PATH ${CMAKE_BINARY_DIR}/ThirdParty/amrex/installdir)

# Check that MFIX_CXX_FLAGS and MFIX_Fortran_FLAGS
# has been defined already
if ( (NOT MFIX_CXX_FLAGS) OR (NOT MFIX_Fortran_FLAGS) )
message (FATAL_ERROR "MFIX_CXX_FLAGS and MFIX_Fortran_FLAGS must be defined \
before including MFIX_ConfigSuperbuild.cmake.")
endif ()


ExternalProject_Add ( amrex
   PREFIX          ${AMREX_SUPERBUILD_DIR}
   INSTALL_DIR     ${AMREX_INSTALL_PATH}
   GIT_REPOSITORY  ${AMREX_GIT_REPO}
   GIT_TAG         ${AMREX_GIT_TAG}
   CMAKE_ARGS
   -DENABLE_MPI=${ENABLE_MPI}
   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
   -DENABLE_PIC=${AMREX_ENABLE_PIC}
   -DENABLE_OMP=${AMREX_ENABLE_OMP}
   -DENABLE_MPI=${AMREX_ENABLE_MPI}
   -DENABLE_DP=${AMREX_ENABLE_DP}
   -DENABLE_PARTICLES=1
   -DENABLE_DP_PARTICLES=${AMREX_ENABLE_DP_PARTICLES}      
   -DBL_SPACEDIM=3
   -DENABLE_LINEAR_SOLVERS=0
   -DENABLE_FBASELIB=1 # Needed for test utilities
   -DENABLE_FORTRAN_INTERFACES=0
   -DENABLE_PROFILING=${AMREX_ENABLE_PROFILING}
   -DENABLE_TINY_PROFILING=${AMREX_ENABLE_TINY_PROFILING}
   -DENABLE_COMM_PROFILING=${AMREX_ENABLE_COMM_PROFILING}      
   -DENABLE_TRACE_PROFILING=${AMREX_ENABLE_TRACE_PROFILING}      
   -DENABLE_BACKTRACE=${AMREX_ENABLE_BACKTRACE}
   -DENABLE_FPE=${ENABLE_FPE}
   -DENABLE_ASSERTIONS=${AMREX_ENABLE_ASSERTIONS}
   -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
   -DAMREX_FFLAGS_OVERRIDES=${MFIX_Fortran_FLAGS}
   -DAMREX_CXXFLAGS_OVERRIDES=${MFIX_CXX_FLAGS}
   # LOG_CONFIGURE 1
   # LOG_BUILD 1
   # LOG_INSTALL 1
   USES_TERMINAL_DOWNLOAD 1
   USES_TERMINAL_CONFIGURE 1
   USES_TERMINAL_BUILD 1
   USES_TERMINAL_INSTALL 1
   )

# Since we cannot rely on a cmake config file to discover
# external libraries and defines (amrex has not been installed
# yet), we need to do this manually.

# Pass amrex flags to internal flags
set (ENABLE_PIC ${AMREX_ENABLE_PIC})
set (ENABLE_MPI ${AMREX_ENABLE_MPI})
set (ENABLE_OMP ${AMREX_ENABLE_OMP})
set (ENABLE_DP  ${AMREX_ENABLE_DP})
set (ENABLE_DP_PARTICLES  ${AMREX_ENABLE_DP_PARTICLES})
set (ENABLE_PROFILING  ${AMREX_ENABLE_PROFILING})
set (ENABLE_TINY_PROFILING  ${AMREX_ENABLE_TINY_PROFILING})
set (ENABLE_TRACE_PROFILING ${AMREX_ENABLE_TRACE_PROFILING})
set (ENABLE_COMM_PROFILING  ${AMREX_ENABLE_COMM_PROFILING})
set (ENABLE_BACKTRACE  ${AMREX_ENABLE_BACKTRACE})
set (ENABLE_ASSERTIONS  ${AMREX_ENABLE_ASSERTIONS})
set (ENABLE_FORTRAN_INTERFACES  0)
set (ENABLE_LINEAR_SOLVERS 0)
set (ENABLE_FBASELIB 1)
set (AMREX_BUILD_TYPE  ${MFIX_BUILD_TYPE})
set (BL_SPACEDIM 3)
set (ENABLE_FBASELIB  0 )
set (ENABLE_PARTICLES 1 )


# ------------------------------------------------------------- #
#    Set preprocessor flags 
# ------------------------------------------------------------- #
add_define (BL_NOLINEVALUES AMREX_DEFINES)
add_define (AMREX_NOLINEVALUES AMREX_DEFINES)
add_define (BL_PARALLEL_IO AMREX_DEFINES)
add_define (AMREX_PARALLEL_IO AMREX_DEFINES)
add_define (BL_SPACEDIM=${BL_SPACEDIM} AMREX_DEFINES)
add_define (AMREX_SPACEDIM=${BL_SPACEDIM} AMREX_DEFINES)
add_define (BL_${CMAKE_SYSTEM_NAME} AMREX_DEFINES)
add_define (AMREX_${CMAKE_SYSTEM_NAME} AMREX_DEFINES)

if ( ENABLE_DP )
   add_define (BL_USE_DOUBLE AMREX_DEFINES)
   add_define (AMREX_USE_DOUBLE AMREX_DEFINES)
else ()
   add_define (BL_USE_FLOAT AMREX_DEFINES)
   add_define (AMREX_USE_FLOAT AMREX_DEFINES)
endif ()

add_define (USE_PARTICLES AMREX_DEFINES ENABLE_PARTICLES)

add_define (BL_PROFILING AMREX_DEFINES ENABLE_PROFILING)
add_define (AMREX_PROFILING AMREX_DEFINES ENABLE_PROFILING) 

add_define (BL_TINY_PROFILING AMREX_DEFINES ENABLE_TINY_PROFILING)
add_define (AMREX_TINY_PROFILING AMREX_DEFINES ENABLE_TINY_PROFILING)

add_define (BL_COMM_PROFILING AMREX_DEFINES ENABLE_COMM_PROFILING)
add_define (AMREX_COMM_PROFILING AMREX_DEFINES ENABLE_COMM_PROFILING)

if ( ENABLE_PARTICLES AND ( NOT ENABLE_DP_PARTICLES ) ) 
   add_define ( BL_SINGLE_PRECISION_PARTICLES AMREX_DEFINES )
   add_define ( AMREX_SINGLE_PRECISION_PARTICLES AMREX_DEFINES )
endif ()

if (ENABLE_FBASELIB)
   add_define ( BL_USE_FORTRAN_MPI AMREX_DEFINES ENABLE_MPI)
   add_define ( AMREX_USE_FORTRAN_MPI AMREX_DEFINES ENABLE_MPI)
   add_define ( BL_USE_F_BASELIB  AMREX_DEFINES )
   add_define ( AMREX_USE_F_BASELIB  AMREX_DEFINES )
endif ()

add_define (BL_USE_MPI AMREX_DEFINES ENABLE_MPI)
add_define (AMREX_USE_MPI AMREX_DEFINES ENABLE_MPI)

add_define (BL_USE_OMP AMREX_DEFINES ENABLE_OMP)
add_define (AMREX_USE_OMP AMREX_DEFINES ENABLE_OMP)

add_define (BL_USE_F_INTERFACES AMREX_DEFINES ENABLE_FORTRAN_INTERFACES)

add_define (BL_USE_ASSERTION AMREX_DEFINES ENABLE_ASSERTIONS) 




#
# Add all preprocessor definitions to compile string
# 
add_definitions ( ${AMREX_DEFINES} )


#
# Detect Fortran name mangling scheme for C/Fortran interface 
#
include(FortranCInterface)
include(${FortranCInterface_BINARY_DIR}/Output.cmake)

if ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "" )  AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "UPPER") )
   message(STATUS "Fortran name mangling scheme to UPPERCASE \
(upper case, no append underscore)")
   add_define ( BL_FORT_USE_UPPERCASE AMREX_DEFINES )
   add_define ( AMREX_FORT_USE_UPPERCASE AMREX_DEFINES )
elseif ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "") AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER") )
   message(STATUS "Fortran name mangling scheme to LOWERCASE \
(lower case, no append underscore)")
   add_define ( BL_FORT_USE_LOWERCASE AMREX_DEFINES )
   add_define ( AMREX_FORT_USE_LOWERCASE AMREX_DEFINES )
elseif ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "_" ) AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER") )
   message(STATUS "Fortran name mangling scheme to UNDERSCORE \
(lower case, append underscore)")
   add_define ( BL_FORT_USE_UNDERSCORE AMREX_DEFINES )
   add_define ( AMREX_FORT_USE_UNDERSCORE AMREX_DEFINES )
else ()
   message(AUTHOR_WARNING "Fortran to C mangling not backward\
 compatible with older style BoxLib code") 
endif ()

#
# We set the AMReX preprocessor flags according to the options
# set by user and amrex defaults, i.e. mfix user can modify only
# a subset of all the amrex options so for all the others we use defaults
#  and here we set defines accordingly.
add_define (BL_NOLINEVALUES AMREX_DEFINES)
add_define (AMREX_NOLINEVALUES AMREX_DEFINES)
add_define (BL_PARALLEL_IO AMREX_DEFINES)
add_define (AMREX_PARALLEL_IO AMREX_DEFINES)
add_define (BL_SPACEDIM=${BL_SPACEDIM} AMREX_DEFINES)
add_define (AMREX_SPACEDIM=${BL_SPACEDIM} AMREX_DEFINES)
add_define (BL_${CMAKE_SYSTEM_NAME} AMREX_DEFINES)
add_define (AMREX_${CMAKE_SYSTEM_NAME} AMREX_DEFINES)

if ( ENABLE_DP )
   add_define (BL_USE_DOUBLE AMREX_DEFINES)
   add_define (AMREX_USE_DOUBLE AMREX_DEFINES)
else ()
   add_define (BL_USE_FLOAT AMREX_DEFINES)
   add_define (AMREX_USE_FLOAT AMREX_DEFINES)
endif ()

add_define (USE_PARTICLES AMREX_DEFINES ENABLE_PARTICLES)
add_define (BL_PROFILING AMREX_DEFINES ENABLE_PROFILING)
add_define (AMREX_PROFILING AMREX_DEFINES ENABLE_PROFILING) 
add_define (BL_TINY_PROFILING AMREX_DEFINES ENABLE_TINY_PROFILING)
add_define (AMREX_TINY_PROFILING AMREX_DEFINES ENABLE_TINY_PROFILING)
add_define (AMREX_ENABLE_ASSERTIONS AMREX_DEFINES ENABLE_)

if ( ENABLE_PARTICLES AND ( NOT ENABLE_DP_PARTICLES ) ) 
   add_define ( BL_SINGLE_PRECISION_PARTICLES AMREX_DEFINES )
   add_define ( AMREX_SINGLE_PRECISION_PARTICLES AMREX_DEFINES )
endif ()

add_define (BL_USE_FORTRAN_MPI=1 AMREX_DEFINES ENABLE_MPI)
add_define (AMREX_USE_FORTRAN_MPI=1 AMREX_DEFINES ENABLE_MPI)
add_define (BL_USE_MPI AMREX_DEFINES ENABLE_MPI)
add_define (BL_USE_OMP AMREX_DEFINES ENABLE_OMP)
add_define (AMREX_USE_MPI AMREX_DEFINES ENABLE_MPI)
add_define (AMREX_USE_OMP AMREX_DEFINES ENABLE_OMP)
add_define (BL_USE_ASSERTION AMREX_DEFINES ENABLE_ASSERTIONS)

add_definitions ( ${AMREX_DEFINES} )

# ------------------------------------------------------------- #
#    Setup third party packages 
# ------------------------------------------------------------- #
set (AMREX_INCLUDE_DIR ${AMREX_INSTALL_PATH}/include)


if (ENABLE_MPI)
   find_package (MPI REQUIRED)
   # Includes
   list (APPEND AMREX_EXTRA_Fortran_INCLUDE_PATH "${MPI_Fortran_INCLUDE_PATH}")
   list (APPEND AMREX_EXTRA_C_INCLUDE_PATH "${MPI_C_INCLUDE_PATH}")
   list (APPEND AMREX_EXTRA_CXX_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}")
   # Compile flags
   append ( MPI_Fortran_COMPILE_FLAGS AMREX_EXTRA_Fortran_FLAGS ) 
   append ( MPI_C_COMPILE_FLAGS AMREX_EXTRA_C_FLAGS )
   append ( MPI_CXX_COMPILE_FLAGS AMREX_EXTRA_CXX_FLAGS )
   # Libraries
   list (APPEND AMREX_EXTRA_Fortran_LIBRARIES "${MPI_Fortran_LIBRARIES}")
   list (APPEND AMREX_EXTRA_C_LIBRARIES "${MPI_C_LIBRARIES}")
   list (APPEND AMREX_EXTRA_CXX_LIBRARIES "${MPI_CXX_LIBRARIES}")
   # Link flags
   list (APPEND AMREX_EXTRA_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
   list (APPEND AMREX_EXTRA_C_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
   list (APPEND AMREX_EXTRA_CXX_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
   # Link Line
   append_to_link_line ( AMREX_EXTRA_Fortran_LIBRARIES
      AMREX_EXTRA_Fortran_LINK_LINE AMREX_EXTRA_Fortran_LINK_FLAGS )
   append_to_link_line ( AMREX_EXTRA_C_LIBRARIES
      AMREX_EXTRA_C_LINK_LINE AMREX_EXTRA_C_LINK_FLAGS )
   append_to_link_line ( AMREX_EXTRA_CXX_LIBRARIES
      AMREX_EXTRA_CXX_LINK_LINE AMREX_EXTRA_CXX_LINK_FLAGS )
endif ()

if (ENABLE_OMP)
   find_package (OpenMP REQUIRED)
   # Compile flags
   append ( OpenMP_Fortran_FLAGS AMREX_EXTRA_Fortran_FLAGS ) 
   append ( OpenMP_C_FLAGS AMREX_EXTRA_C_FLAGS )
   append ( OpenMP_CXX_FLAGS AMREX_EXTRA_CXX_FLAGS )
endif()

append ( AMREX_EXTRA_Fortran_FLAGS AMREX_Fortran_FLAGS )
append ( AMREX_EXTRA_CXX_FLAGS AMREX_CXX_FLAGS )

# Set amrex libs
set (AMREX_LIBRARIES     amrex)
set (AMREX_LIBRARY_DIR   ${AMREX_INSTALL_PATH}/lib)
