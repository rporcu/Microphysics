cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
enable_language(Fortran)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if (CMAKE_Fortran_COMPILER MATCHES ".*gfortran.*")
   set (CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -g -Wall -O2" CACHE STRING "" FORCE)
   set (CMAKE_CXX_FLAGS_DEBUG   "-fmax-errors=2 -std=c++11 -g -Wall -Wextra -O0" CACHE STRING "" FORCE)
   set (CMAKE_Fortran_FLAGS_RELEASE "-g -std=f2008 -fcheck=all -Wall -O2" CACHE STRING "" FORCE)
   set (CMAKE_Fortran_FLAGS_DEBUG   "-fmax-errors=2 -g -std=f2008 -fcheck=all -Wall -Wextra -Wimplicit-interface" 
      CACHE STRING "" FORCE)

endif (CMAKE_Fortran_COMPILER MATCHES ".*gfortran.*")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set some defaults
set(BL_SPACEDIM 3 CACHE INT "Dimension of AMReX build")
set(ENABLE_MPI 0 CACHE INT "Enable build with MPI")
set(ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
set(BL_PRECISION "DOUBLE" CACHE INT "Precision of AMReX build")
set(BL_USE_PARTICLES 1 CACHE INT "Include Particles classes in AMReX build")
set(ENABLE_PROFILING 0 CACHE INT "Include profiling information in AMReX build")
set(ENABLE_BACKTRACE 1 CACHE INT "Include backtrace information in AMReX build")
set(ENABLE_SUPERBUILD 1 CACHE INT  "Allow CMake to download, build and install AMReX")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)


# Load MFIX Cmake Tools
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Tools/)
include(MFIXTools)
include(AMReX)

project(MFIXEXA)

# This is the main file
set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp )


# Find CCSE check for mpi libraries if ENABLE_MPI is on.
# Thus all the variables defined by FindMPI are available
# Notice that only need c++-related variables are needed, since
# mfix itself does not call mpi function explicitly 
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${CCSE_INCLUDE_DIR})
include_directories(${AMREX_INCLUDES})
# include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(src)
include_directories(src/include)
link_directories(${AMREX_LIB_DIR})
# link_directories(${CCSE_LIBRARY_DIR})



# Main 
add_subdirectory( src )

FILE(GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
add_executable(mfix ${MAIN} ${USR_OVERRIDES})
target_link_libraries( mfix mfixcore ${AMREX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES} )

# Tests
enable_testing()
add_subdirectory( tests )

