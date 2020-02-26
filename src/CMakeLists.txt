add_library(mfixcore)

#
# Path where to store module files
#
set( MFIX_Fortran_MOD_PATH ${CMAKE_CURRENT_BINARY_DIR}/mod_files
   CACHE INTERNAL "Fortran modules location")

set_target_properties( mfixcore
   PROPERTIES
   Fortran_MODULE_DIRECTORY
   ${MFIX_Fortran_MOD_PATH}
   )

target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${MFIX_Fortran_MOD_PATH}>)


#
# Add sources and header from current dir and all sub dirs
#

# Uncategorized
target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

# IMPORTANT: DO NOT ADD main.cpp HERE!
target_sources(mfixcore
   PRIVATE
   mfix_compute_diveu.cpp
   mfix_regrid.cpp
   mfix.cpp
   #  main.cpp <---- DO NOT ADD THIS HERE
   mfix_F.H
   mfix.H
   )

# Fortran include files
add_subdirectory(include)

# I/0
add_subdirectory(io)

# Convection
add_subdirectory(convection)

# DES
add_subdirectory(deposition)
add_subdirectory(des)

# Diffusion
add_subdirectory(diffusion)

# LevelData
add_subdirectory(level_data)

# Projection
add_subdirectory(projection)

# Setup
add_subdirectory(setup)

# Utility functions
add_subdirectory(util)

# Global module files
add_subdirectory(mods)

# User hooks (UDFs)
add_subdirectory(usr)

# Embedded boundary functions
add_subdirectory(eb)

# Time evolution functions
add_subdirectory(time_flow)

#
# Link to imported AMReX
#
target_link_libraries( mfixcore PUBLIC AMReX::amrex )


#
# If ENABLE_CUDA, set CUDA-specific properties and
# make C++ files be compiled as CUDA sources
#
if (ENABLE_CUDA)
   set_target_properties( mfixcore
      PROPERTIES
      CUDA_SEPARABLE_COMPILATION ON      # This adds -dc
      CUDA_RESOLVE_DEVICE_SYMBOLS OFF
      )

   get_target_property(MFIX_CUDA_SOURCES mfixcore SOURCES)
   list(FILTER MFIX_CUDA_SOURCES INCLUDE REGEX "\\.cpp")
   set_source_files_properties(${MFIX_CUDA_SOURCES} PROPERTIES LANGUAGE CUDA)
   unset(MFIX_CUDA_SOURCES)

   #
   # Add manually curand
   # CMake >= 3.17 provides the module FindCUDAToolkit to do this natively.
   #
   find_library(CURAND curand PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
   target_link_libraries(mfixcore PUBLIC ${CURAND})
endif ()

#
# Set the default flags
#
if ( NOT CMAKE_Fortran_FLAGS )
   target_link_libraries(mfixcore PUBLIC AMReX::Flags_Fortran)
endif ()

if ( NOT CMAKE_CXX_FLAGS )
   target_link_libraries(mfixcore PUBLIC AMReX::Flags_CXX)
endif ()

if (ENABLE_FPE)
   target_link_libraries(mfixcore PUBLIC AMReX::Flags_FPE)
endif ()


#
# Add typecheck target for mfixcore
#
include( AMReXTypecheck )
add_typecheck_target(mfixcore)


#
# Add AMReX "build info" capabilities
#
include(AMReXBuildInfo)
generate_buildinfo(mfixcore ${CMAKE_SOURCE_DIR} REQUIRED)
