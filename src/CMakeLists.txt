#
# Define a macro to add sources.
# This will also accumulate include directories if
# a header file is spotted
# WARNING: this automagically prepends the correct
# path (relative to the file invoking the macro)
# to the source file.
#
set ( MFIX_INCLUDE_PATHS  ${CMAKE_CURRENT_LIST_DIR}/include )

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
add_sources ( mfix_evolve.cpp mfix_init.cpp mfix_arrays.cpp )
add_sources ( mfix_level.cpp mfix_regrid.cpp )
add_sources ( mfix_level.H mfix_F.H )
add_sources ( calc_cell.f90 calc_mu_g.f90 )
add_sources ( calc_tau.f90 compute_divu.f90 compute_vort.f90 )
add_sources ( fill_bc0.f90 )
add_sources ( get_data.f90 init_fluid.f90 init_namelist.f90 )
add_sources ( set_bc0.f90 set_bc_flow.f90 set_bc_type.f90 set_p0.f90 set_delp_dir.f90 )
add_sources ( set_ps.f90 set_velocity_bcs.f90 )
add_sources ( zero_wall_norm_vel.f90 amrex_to_mfix.F90 )

#
# Fortran Modules
#
include ( mods/CMakeLists.txt )

#
# EB stuff
#
include ( eb/CMakeLists.txt )

#
# Particles
#
include ( des/CMakeLists.txt )
include ( des_fluid/CMakeLists.txt )

#
# SIMPLE routines
#
include ( simple/CMakeLists.txt )

#
# Data checks
#
include ( check_data/CMakeLists.txt )

#
# I/O routines
#
include ( io/CMakeLists.txt )

#
# Users hooks
#
include ( usr/CMakeLists.txt )

#
# Projection
#
include ( projection/CMakeLists.txt )

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
