cmake_minimum_required(VERSION 3.14)

project(MFIX-Exa
   # VERSION <version>
   DESCRIPTION  "A multiphase modeling tool for the exascale"
   HOMEPAGE_URL "https://amrex-codes.github.io/MFIX-Exa/"
   LANGUAGES    C CXX
   )

# Set the search path for cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/tools/CMake )

#
# This sets CMAKE_CXX_FLAGS_<CONFIG> to a default value
# if the variable is empty
#
if( NOT CMAKE_CXX_FLAGS_DEBUG )
   set(CMAKE_CXX_FLAGS_DEBUG "-g")
endif()

if( NOT CMAKE_CXX_FLAGS_RELEASE )
   set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
endif()

# Load utilities and MFIX options module
include( MFIX_Utils )

# Creates `compile_commands.json` in the build build directory
#  `compile_commands.json` contains compiler flags used by plugins like YouCompleteMe
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# Set default build type to Release
if ( NOT CMAKE_BUILD_TYPE )
   message(STATUS "Setting build type to Release as none was specified.")
   set( CMAKE_BUILD_TYPE Release )
else ()
   message(STATUS "Build type set by user to '${CMAKE_BUILD_TYPE}'.")
endif()

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set( CMAKE_CXX_COMPILER_LAUNCHER ccache )
  set( CMAKE_CUDA_COMPILER_LAUNCHER ccache )
endif()

# We want to offer the user some core options, i.e. MPI, OpenMP, CUDA, HYPRE
# and floating-point exceptions, which require related AMReX components to
# be enabled in the installation.
# In STANDALONE mode, find_package will make sure these options are available
# IN SUPERBUILD mode, These options will supersede AMReX options with the same name
option( MFIX_OMP    "Enable OpenMP" NO )
option( MFIX_MPI    "Enable MPI"   YES )

#
# Temporary hack
# We need AMReX_FORTRAN=on (default, equivalent to BL_NO_FORT=FALSE)  for amrex.
# In STANDALONE mode, this implies that that the AMReX config file will try to
# search for the Fortran OpenMP library but this will fail unless the
# Fortran language is enabled on mfix' side.
#
if (MFIX_OMP)
   enable_language(Fortran)
endif ()

include(CMakeDependentOption)
cmake_dependent_option( MFIX_MPI_THREAD_MULTIPLE
   "whether to initialize MPI so that multiple threads can make MPI calls at the same time"  OFF
   "MFIX_MPI" OFF)

option( MFIX_HYPRE  "Enable HYPRE"  NO )
option( MFIX_FPE    "Enable Floating Point Exceptions checks" NO )

#
# GPU backends    =============================================================
#
set(MFIX_GPU_BACKEND_VALUES NONE SYCL CUDA HIP)
set(MFIX_GPU_BACKEND NONE CACHE STRING "On-node, accelerated GPU backend: <NONE,SYCL,CUDA,HIP>")
set_property(CACHE MFIX_GPU_BACKEND PROPERTY STRINGS ${MFIX_GPU_BACKEND_VALUES})
if (NOT MFIX_GPU_BACKEND IN_LIST MFIX_GPU_BACKEND_VALUES)
   message(FATAL_ERROR "MFIX_GPU_BACKEND=${MFIX_GPU_BACKEND} is not allowed."
      " Must be one of ${MFIX_GPU_BACKEND_VALUES}")
endif ()

if (NOT MFIX_GPU_BACKEND STREQUAL NONE)
   message( STATUS "   MFIX_GPU_BACKEND = ${MFIX_GPU_BACKEND}")
endif ()

if (MFIX_GPU_BACKEND STREQUAL SYCL)
   set(MFIX_DPCPP ON  )
   set(MFIX_CUDA  OFF )
   set(MFIX_HIP   OFF )
elseif (MFIX_GPU_BACKEND STREQUAL CUDA)
   set(MFIX_DPCPP OFF )
   set(MFIX_CUDA  ON  )
   set(MFIX_HIP   OFF )
elseif (MFIX_GPU_BACKEND STREQUAL HIP)
   set(MFIX_DPCPP OFF )
   set(MFIX_CUDA  OFF )
   set(MFIX_HIP   ON  )
elseif ()
   set(MFIX_DPCPP OFF )
   set(MFIX_CUDA  OFF )
   set(MFIX_HIP   OFF )
endif ()

# CSG support
option( MFIX_CSG "Build with CSG support" OFF)

if (MFIX_CUDA)
    enable_language(CUDA)
endif ()

# Setup AMReX
include(MFIXSetupAMReX)

get_git_info( )

# Set some cache variable for the whole project
set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp CACHE INTERNAL "Main source file")

if(MFIX_CSG)
  add_subdirectory(subprojects/csg-eb)
endif()

#
# Core MFIX library
#
add_subdirectory(src)

# Print out config summary
print_mfix_configuration_summary(mfixcore)

add_executable(mfix ${MAIN})
target_link_libraries( mfix mfixcore )

# We should check what this line does
file( GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.cpp)
add_udfs_to_target(mfix ${USR_OVERRIDES})

if (MFIX_CUDA)
   setup_target_for_cuda_compilation(mfix)
endif ()

# Tools
add_subdirectory(tools)

# Tests
enable_testing()
add_subdirectory( tests )

# Tags
include(MFIX_Tags)
add_tags_targets()

# Add a target to generate API documentation with Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND)
   set(DOC_DIR ${CMAKE_SOURCE_DIR}/doc)
   configure_file(${DOC_DIR}/Doxyfile.in ${DOC_DIR}/Doxyfile)
   add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOC_DIR}/Doxyfile
      WORKING_DIRECTORY ${DOC_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
endif (DOXYGEN_FOUND)

# Copy compile_commands.json from the build directory to the project root directory
if ( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
endif ()

if(MFIX_CSG)
  target_link_libraries(mfixcore PRIVATE csg)
endif()
