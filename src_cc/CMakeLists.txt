#
# Define a macro to add sources.
# This will also accumulate include directories if
# a header file is spotted
# WARNING: this automagically prepends the correct
# path (relative to the file invoking the macro)
# to the source file.
#
set ( MFIX_INCLUDE_PATHS  ${PROJECT_SOURCE_DIR}/src/include )

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
add_sources ( set_bc0.f90   set_delp_dir.f90  set_velocity_bcs.f90 )

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

# 
# Add the include paths
# By adding the include paths as a property of mfixcore, we make sure
# that these same paths will be used whenever another target target
# lists mficore as a dependency
# 
target_include_directories ( ${MFIX_LIBNAME} PUBLIC
  ${MFIX_INCLUDE_PATHS} ${MFIX_EXTRA_INCLUDE_PATH} )
