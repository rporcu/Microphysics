# # This file perform the Superbuild
# Both amrex and mfix are treated as external project
# This allows to build amrex and only after to handle
# mfix ( and read the amrex config file )
#
project ( MFIX-Exa_Superbuild )

message (STATUS "SUPERBUILD mode: AMReX will be built as part of MFIX")

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
# Amrex-related variables
#

# AMReX Git variables
set (AMREX_GIT_REPO "https://github.com/AMReX-Codes/amrex.git" )
set (AMREX_GIT_COMMIT_MASTER  777a7c4241cb82437e06aabd7dfe0a32d4 )
set (AMREX_GIT_COMMIT_DEVELOP 01a4ac4249302e6a8d19402b4b1c70ca0b19fc70 )

set (AMREX_GIT_TAG)  # The commit id or branch to download

#
# MFIX-related options
#
include ( MFIX_Options )

#
# AMReX-related config options
#
set ( AMREX_Fortran_FLAGS "" CACHE STRING "User-defined Fortran compiler flags for AMReX (Superbuild only)" )

set ( AMREX_CXX_FLAGS "" CACHE STRING "User-defined C++ compiler flags for AMReX (Superbuild only)" )

option ( AMREX_ENABLE_EB "Build EB code" ON)

option ( AMREX_ENABLE_PIC "Build position-independent code" OFF)

option ( AMREX_ENABLE_MPI  "Enable MPI"  ON)

option ( AMREX_ENABLE_OMP  "Enable OpenMP" OFF)

option ( AMREX_ENABLE_DP "Enable double precision" ON)

option ( AMREX_ENABLE_ASSERTION "Enable assertions" OFF)

option ( AMREX_ENABLE_BASE_PROFILE "Enable basic profiling" OFF )

option ( AMREX_ENABLE_TINY_PROFILE "Enable tiny-profiling" ON)

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
get_git_info ( MFIX_GIT_BRANCH MFIX_GIT_COMMIT )

if (AMREX_GIT_COMMIT)
   set (AMREX_GIT_TAG ${AMREX_GIT_COMMIT})
else ()
   if (MFIX_GIT_BRANCH MATCHES "master")
      set (AMREX_GIT_TAG ${AMREX_GIT_COMMIT_MASTER})
   else ()
      set (AMREX_GIT_TAG ${AMREX_GIT_COMMIT_DEVELOP})
   endif()
endif ()

message (STATUS "AMReX commit: ${AMREX_GIT_TAG}")

# Include cmake config files to build external projects
include(ExternalProject)

# This modify the dir structure in external project directory
set_directory_properties ( PROPERTIES EP_BASE ${CMAKE_BINARY_DIR}/amrex )

set (AMREX_SUPERBUILD_INSTALLDIR ${CMAKE_BINARY_DIR}/amrex/installdir)
set (AMREX_SUPERBUILD_BUILDDIR   ${CMAKE_BINARY_DIR}/amrex/builddir)
set (AMREX_SUPERBUILD_SOURCEDIR  ${CMAKE_BINARY_DIR}/amrex/sourcedir)

ExternalProject_Add ( amrex
   INSTALL_DIR     ${AMREX_SUPERBUILD_INSTALLDIR}
   BINARY_DIR      ${AMREX_SUPERBUILD_BUILDDIR}
   SOURCE_DIR      ${AMREX_SUPERBUILD_SOURCEDIR}
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
   -DENABLE_AMRDATA=ON # Needed to compile postprocessing tools
   -DENABLE_LINEAR_SOLVERS=ON
   -DENABLE_EB=${AMREX_ENABLE_EB}
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
   -DCMAKE_Fortran_FLAGS=${AMREX_Fortran_FLAGS}
   -DCMAKE_CXX_FLAGS=${AMREX_CXX_FLAGS}
   -DCMAKE_EXPORT_COMPILE_COMMANDS=${CMAKE_EXPORT_COMPILE_COMMANDS}
   UPDATE_COMMAND ""
   BUILD_ALWAYS 1
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
set (MFIX_SUPERBUILD_BUILDDIR   ${CMAKE_BINARY_DIR}/mfix)


ExternalProject_Add ( mfix
   PREFIX          ${MFIX_SUPERBUILD_BUILDDIR}
   DEPENDS amrex
   CMAKE_ARGS
   -DDEBUG=${DEBUG}
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
   -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
   -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
   -DENABLE_FPE=${ENABLE_FPE}
   -DENABLE_PROJCC=${ENABLE_PROJCC}
   -DENABLE_PTESTS=${ENABLE_PTESTS}
   -DENABLE_STESTS=${ENABLE_STESTS}
   -DAMREX_INSTALL_DIR=${AMREX_SUPERBUILD_INSTALLDIR}
   -DCMAKE_EXPORT_COMPILE_COMMANDS=${CMAKE_EXPORT_COMPILE_COMMANDS}
   SOURCE_DIR ${PROJECT_SOURCE_DIR}
   BINARY_DIR ${MFIX_SUPERBUILD_BUILDDIR}
   USES_TERMINAL_CONFIGURE 1
   USES_TERMINAL_BUILD 1
   INSTALL_COMMAND ""
   )

# When using superbuild, the compile commands databases do not exist before
# compile-time.  Hence create a new build target (compile_db) which collects
# both mfix's and amrex's compile_commands.json, concatenates them in the
# project's root directory
add_custom_target( compile_db
    # First take mfix's compile database (compile_commands.json) and:
    # 1.  $d     => deletes (d) last line ($) of a file (in this case the trailing ])
    # 2.a $      => jumps to end ($) of file, then apply:
    # 2.b s/$/,/ => string search-and-replace (s/.../.../) for the and of the line ($) there add a comma (,)
    COMMAND sed "$d" ${MFIX_SUPERBUILD_BUILDDIR}/compile_commands.json | sed "$ s/$/,/" > ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
    # Now take amrex's compile database (compile_commands.json) and:
    # 1. 1d => deletes the first line (in this case the leading [)
    # 2. >> => appends the the end of mfix's compile database in the root-dir
    COMMAND sed "1d" ${AMREX_SUPERBUILD_BUILDDIR}/compile_commands.json >> ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
    VERBATIM
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/
    )
