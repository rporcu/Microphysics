#
# Define a macro to add sources.
# This will also accumulate include directories if
# a header file is spotted
# WARNING: this automagically prepends the correct
# path (relative to the file invoking the macro)
# to the source file.
#
set ( MFIX_Fortran_MOD_PATH ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod )
set ( MFIX_INCLUDE_PATHS  ${PROJECT_SOURCE_DIR}/src/include ${MFIX_Fortran_MOD_PATH})

macro ( add_sources )
  foreach ( item IN ITEMS ${ARGN} )
    set ( source_name  ${CMAKE_CURRENT_LIST_DIR}/${item} )
    string ( REPLACE "//" "/" source_name ${source_name})
    target_sources ( ${MFIX_LIBNAME} PRIVATE ${source_name} )

    # Add folder to list of includes if necessary
    get_filename_component ( filetype ${source_name} EXT )
    if ( (${filetype} MATCHES ".h") OR (${filetype} MATCHES ".H" ) )
      list (APPEND MFIX_INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR})
    endif()
  endforeach()
endmacro()

#
# Fist, define the library we want to add
#
add_library ( ${MFIX_LIBNAME} "" )

set_target_properties ( ${MFIX_LIBNAME}
   PROPERTIES
   Fortran_MODULE_DIRECTORY
   ${MFIX_Fortran_MOD_PATH}
   )

#
# Now, one by one, we add all the sources from all the subdirs
# to the library defined above
# 

#
# Uncategorized
# 
add_sources ( amrex_to_mfix.F90 fill_bc0.f90 init_namelist.f90 mfix_init.cpp  mfix_regrid.cpp mfix_arrays.cpp)
add_sources ( set_bc_flow.f90   set_p0.f90   zero_wall_norm_vel.f90 calc_cell.f90 )
add_sources ( compute_vort.f90  get_data.f90   mfix_evolve.cpp  mfix_level.cpp  set_bc_type.f90 )
add_sources ( set_ps.f90 calc_mu_g.f90  init_fluid.f90  main.cpp  mfix_F.H mfix_level.H )
add_sources ( set_epg.f90 set_bc0.f90   set_delp_dir.f90  set_velocity_bcs.f90 )

#
# I/0
# 
include ( io/CMakeLists.txt )

#
# DES
#
include ( des_fluid/CMakeLists.txt )

#
# Projection
#
include ( projection/CMakeLists.txt )


#
# Add folders from src ( this is the part of code base in commons between
# src and src_cc )
#
set ( SRCDIR ${PROJECT_SOURCE_DIR}/src )

include ( ${SRCDIR}/check_data/CMakeLists.txt )
include ( ${SRCDIR}/des/CMakeLists.txt )
include ( ${SRCDIR}/eb/CMakeLists.txt )
include ( ${SRCDIR}/mods/CMakeLists.txt )
include ( ${SRCDIR}/usr/CMakeLists.txt )


# Remove duplicates from include lists
list ( REMOVE_DUPLICATES MFIX_INCLUDE_PATHS )
target_include_directories ( ${MFIX_LIBNAME} PUBLIC ${MFIX_INCLUDE_PATHS} )

#
# Link to imported AMReX
#
target_link_libraries ( ${MFIX_LIBNAME} AMReX::amrex )


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
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -fpe3>>
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -fpe3>>
	 # Cray
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -K trap=fp>>
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -K trap=fp>>
	 #  PGI
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -Ktrap=divz,inv>>
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 >> )
   endif ()
else ()
   message (AUTHOR_WARNING "Variable ENABLE_FPE is not defined")
endif ()
