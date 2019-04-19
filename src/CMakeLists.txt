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

# Uncategorized IMPORTANT: DO NOT ADD main.cpp HERE!
target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

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
   
   get_target_properties(MFIX_CUDA_SOURCES mfixcore SOURCES)
   list(FILTER MFIX_CUDA_SOURCES INCLUDE REGEX "\\.cpp")
   set_source_files_properties(${MFIX_CUDA_SOURCES} PROPERTIES LANGUAGE CUDA)
   unset(MFIX_CUDA_SOURCES)  
endif ()

#
# Set the default flags 
# 
if ( NOT CMAKE_Fortran_FLAGS )
   target_compile_options ( ${MFIX_LIBNAME}
      PUBLIC
      # GNU Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
      -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -finit-real=snan -finit-integer=2147483647>>>>
      # GNU Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
      >>>>
      # Intel Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:Fortran>:
      -O0 -traceback -check bounds,uninit,pointers>>>>
      # Intel Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:Fortran>:
      -ip -qopt-report=5 -qopt-report-phase=vec>>>>
      # Cray Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
      -O0 -e i>>>>
      # Cray Release 
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
      >>>>
      # PGI Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:Fortran>:
      -O0 -Mbounds -Mchkptr>>>>
      # PGI Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:Fortran>:
      -gopt -fast>>>>
      )	  
endif ()

if ( NOT CMAKE_CXX_FLAGS )
   target_compile_options ( ${MFIX_LIBNAME}
      PUBLIC
      # GNU Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
      -O0 -fno-inline -ggdb -Wall -Wno-sign-compare>>>>
      # GNU Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
      >>>>
      # Intel Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
      -O0 -traceback -Wcheck>>>>
      # Intel Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
      -ip -qopt-report=5 -qopt-report-phase=vec>>>>
      # Cray Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
      -O0>>>>
      # Cray Release 
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
      >>>>
      # PGI Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
      -O0 -Mbounds>>>>
      # PGI Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
      -gopt -fast>>>>
      )	  
endif ()

if (DEFINED ENABLE_FPE)
   if (ENABLE_FPE)
      target_compile_options ( ${MFIX_LIBNAME}
   	 PUBLIC
	 # GNU
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -ffpe-trap=invalid,zero -ftrapv>>
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -ftrapv>>
	 # Intel
	 $<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -fpe3>>
	 $<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
	 -fpe3>>
	 # Cray
	 $<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -K trap=fp>>
	 $<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
	 -K trap=fp>>
	 #  PGI
	 $<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -Ktrap=divz,inv>>
	 $<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
	 >> )
   endif ()
else ()
   message (AUTHOR_WARNING "Variable ENABLE_FPE is not defined")
endif ()


#
# Add typecheck target for mfixcore
#
include( MFIX_Typecheck )
add_typecheck_target(mfixcore)

