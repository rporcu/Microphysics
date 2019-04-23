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
   compute_dt.cpp
   mfix_regrid.cpp
   mfix_advance.cpp
   compute_diveu.cpp
   mfix_drag.cpp
   mfix_evolve.cpp 
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
add_subdirectory(des_fluid)
add_subdirectory(src_des)

# Diffusion
add_subdirectory(diffusion)

# Projection
add_subdirectory(projection)

# Setup
add_subdirectory(setup)

# Utility functions
add_subdirectory(util)

# Setup check functions
add_subdirectory(check_data)

# Global module files
add_subdirectory(mods)

# User hooks (UDFs)
add_subdirectory(usr)

# Embedded boundary functions
add_subdirectory(eb)

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
include( MFIX_Typecheck )
add_typecheck_target(mfixcore)

