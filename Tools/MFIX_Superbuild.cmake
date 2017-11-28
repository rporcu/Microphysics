#
# This file perform the Superbuild
# Both amrex and mfix are treated as external project
# This allows to build amrex and only after to handle
# mfix ( and read the amrex config file )
# 
project ( MFIX-Exa_Superbuild )

#
# Load Utilities
# 
include ( MFIX_Utils )

#
# Define the languages used by the project
#
enable_language (C)
enable_language (CXX)
enable_language (Fortran)

#
# Require C++11 standard
#
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON) 
set (CMAKE_CXX_EXTENSIONS OFF)


#
# Amrex-related variables
# 

# AMReX Git variables
set (AMREX_GIT_REPO "https://github.com/AMReX-Codes/amrex.git" )
set (AMREX_GIT_COMMIT_MASTER  3506f5aea50d27237dda43df3ba4611fd4eda638 )
set (AMREX_GIT_COMMIT_DEVELOP 79d1e0c0257e85f89d9c55d5a109136a2d70ae08 )
set (AMREX_GIT_TAG)  # The commit id or branch to download 

# AMReX Superbuild variables
set (AMREX_SUPERBUILD_DIR   ${PROJECT_BINARY_DIR})

#
# MFIX-related options
#
include ( MFIX_CMakeVariables )
include ( MFIX_Options )

#
#  Setup core compiler flags 
#
include ( MFIX_Compilers )

if ( MFIX_FFLAGS_OVERRIDES )
   set ( MFIX_Fortran_FLAGS ${MFIX_FFLAGS_OVERRIDES} )
else ()
   append ( MFIX_FFLAGS_${MFIX_BUILD_TYPE}
      MFIX_Fortran_FLAGS )
endif ()

if ( MFIX_CXXFLAGS_OVERRIDES )
   set ( MFIX_CXX_FLAGS ${MFIX_CXXFLAGS_OVERRIDES} )
else ()
   append ( MFIX_CXXFLAGS_${MFIX_BUILD_TYPE}
      MFIX_CXX_FLAGS )
endif ()


#
# AMReX-related config options
# 
option ( AMREX_ENABLE_PIC "Build position-independent code" OFF)

option ( AMREX_ENABLE_MPI  "Enable MPI"  ON)

option ( AMREX_ENABLE_OMP  "Enable OpenMP" OFF)

option ( AMREX_ENABLE_DP "Enable double precision" ON)

option ( AMREX_ENABLE_ASSERTION "Enable assertions" OFF)

option ( AMREX_ENABLE_BASE_PROFILE "Enable basic profiling" OFF )

option ( AMREX_ENABLE_TINY_PROFILE "Enable tiny-profiling" OFF)

option ( AMREX_ENABLE_TRACE_PROFILE "Enable trace-profiling" OFF )

option ( AMREX_ENABLE_MEM_PROFILE   "Enable memory profiling" OFF )

option ( AMREX_ENABLE_COMM_PROFILE  "Enable communicator-profiling" OFF )

option ( AMREX_ENABLE_BACKTRACE "Enable backtracing" OFF)

option ( AMREX_ENABLE_PROFPARSER "Enable profile parser" OFF)

option ( AMREX_ENABLE_DP_PARTICLES "Enable double-precision particle data" ON )

set ( AMREX_GIT_COMMIT "" CACHE STRING "AMReX git commit to use in superbuild")

#
# Set the git commit to use for amrex
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

set (AMREX_INSTALL_PATH ${CMAKE_BINARY_DIR}/amrex/installdir)


ExternalProject_Add ( amrex
   PREFIX          ${AMREX_SUPERBUILD_DIR}
   INSTALL_DIR     ${AMREX_INSTALL_PATH}
   GIT_REPOSITORY  ${AMREX_GIT_REPO}
   GIT_TAG         ${AMREX_GIT_TAG}
   CMAKE_ARGS
   -DENABLE_PIC=${AMREX_ENABLE_PIC}
   -DENABLE_OMP=${AMREX_ENABLE_OMP}
   -DENABLE_MPI=${AMREX_ENABLE_MPI}
   -DENABLE_DP=${AMREX_ENABLE_DP}
   -DENABLE_PARTICLES=ON
   -DENABLE_DP_PARTICLES=${AMREX_ENABLE_DP_PARTICLES}      
   -DDIM=3
   -DDEBUG=${DEBUG}
   -DENABLE_LINEAR_SOLVERS=ON
   -DENABLE_FBASELIB=ON # Needed for test utilities
   -DENABLE_FORTRAN_INTERFACES=OFF
   -DENABLE_BASE_PROFILE=${AMREX_ENABLE_BASE_PROFILE}
   -DENABLE_TINY_PROFILE=${AMREX_ENABLE_TINY_PROFILE}
   -DENABLE_COMM_PROFILE=${AMREX_ENABLE_COMM_PROFILE}      
   -DENABLE_TRACE_PROFILE=${AMREX_ENABLE_TRACE_PROFILE}
   -DENABLE_MEM_PROFILE=${AMREX_ENABLE_MEM_PROFILE}      
   -DENABLE_BACKTRACE=${AMREX_ENABLE_BACKTRACE}
   -DENABLE_FPE=${ENABLE_FPE}
   -DENABLE_ASSERTIONS=${AMREX_ENABLE_ASSERTION}
   -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
   -DAMREX_FFLAGS_OVERRIDES=${MFIX_Fortran_FLAGS}
   -DAMREX_CXXFLAGS_OVERRIDES=${MFIX_CXX_FLAGS}
   UPDATE_COMMAND ""
   # LOG_CONFIGURE 1
   # LOG_BUILD 1
   # LOG_INSTALL 1
   USES_TERMINAL_DOWNLOAD 1
   USES_TERMINAL_CONFIGURE 1
   USES_TERMINAL_BUILD 1
   USES_TERMINAL_INSTALL 1
   )

#
# Now have Cmake call itself to set up mfix
# 
ExternalProject_Add ( mfix-superbuild
   DEPENDS amrex 
   CMAKE_ARGS
   -DDEBUG=${DEBUG}
   -DMFIX_FFLAGS_OVERRIDES=${MFIX_FFLAGS_OVERRIDES}
   -DMFIX_CXXFLAGS_OVERRIDES=${MFIX_CXXFLAGS_OVERRIDES}
   -DENABLE_FPE=${ENABLE_FPE}
   -DAMREX_INSTALL_DIR=${AMREX_INSTALL_PATH}
   SOURCE_DIR ${PROJECT_SOURCE_DIR}
   BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}
   USES_TERMINAL_CONFIGURE 1
   USES_TERMINAL_BUILD 1
   INSTALL_COMMAND ""
   )
