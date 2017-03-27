cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_language(Fortran)

# Load MFIX Cmake Tools
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Tools/)
include(MFIXTools)
include(AMReXConfig)

set_defaults()
set_compiler_flags()

# Config Amrex
if (ENABLE_SUPERBUILD)
   build_amrex()
else (ENABLE_SUPERBUILD)
   find_amrex()
endif (ENABLE_SUPERBUILD)

# Set build flags ( provides MFIX_INCLUDES, 
# MFIX_EXT_LINK_STRING  and MFIX_EXT_LINK_DIRS )
set_build_flags()

# Print summary of config settings
summary()

project(MFIXEXA)

# This is the main file
set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp )

include_directories(${MFIX_INCLUDES})
link_directories(${MFIX_EXT_LINK_DIRS})

# Main 
add_subdirectory( src )

FILE(GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
add_executable(mfix ${MAIN} ${USR_OVERRIDES})
target_link_libraries( mfix mfixcore ${MFIX_EXT_LINK_STRING} )

# Tests
enable_testing()
add_subdirectory( tests )

