#
# Check if dir or file given by path exists and issue a warning or
# error if not
#
function ( check_path  path  message_type )
   if ( EXISTS ${path} )
   else ()
      message(${message_type} ${path} " does not exist!")
   endif ( EXISTS ${path} )
endfunction ()

#
# This function turns a list into a string
# 
function ( list_to_string list )
   string (REPLACE ";" " " tmp "${${list}}")
   set ( ${list} "${tmp}" PARENT_SCOPE)
endfunction ()


#
# Append new_var to all_var
# 
function ( append new_var all_var )
   if ( ${new_var} )
      set ( tmp  "${${all_var}} ${${new_var}}" )

      # Since I am OCD, remove the double spaces.
      string ( REPLACE "  " " " tmp ${tmp} )
      set ( ${all_var}  ${tmp} PARENT_SCOPE )
   endif ()
endfunction ()

#
# Function to append to link line
#
function ( append_to_link_line libs link_line )

   if ( ${ARGC} EQUAL 3 )  # Third one is optional flags
      set ( flags  ${ARGV2} )
   else ()
      set ( flags )
   endif ()
   
   set ( tmp "${${link_line}} ${${flags}} ${${libs}} " )
   string ( STRIP "${tmp}" tmp )
   set ( ${link_line} ${tmp} PARENT_SCOPE )
   
endfunction ()



#
# Function to accumulate preprocessor directives
#
function ( add_define new_define all_defines )
   
   set ( condition  1 )
   
   if ( ${ARGC} EQUAL 3 ) #
      set ( condition ${${ARGV2}} )
   elseif ( ${ARGC} GREATER 3 )
      message ( AUTHOR_WARNING "Function add_define accept AT MOST 3 args" )
   endif ()
   
   if ( ${condition} )
      set ( ${all_defines} "${${all_defines}} -D${new_define}" PARENT_SCOPE )
   endif ()
   
endfunction ()

#
# Stop if in-source build
#
macro ( check_build_tree_path )
   string ( FIND "${CMAKE_BINARY_DIR}" "${PROJECT_BINARY_DIR}/src" RESULT )
   if ( NOT "${RESULT}" STREQUAL "-1")
      message ( FATAL_ERROR  "ERROR: in-source builds are not allowed!")
   endif ()
endmacro( check_build_tree_path )





#
# To be removed as soon as find amrex is fixed
#
macro(set_defaults)

   # Set some defaults
   set(BL_SPACEDIM 3 CACHE INT "Dimension of AMReX build")
   set(ENABLE_MPI 0 CACHE INT "Enable build with MPI")
   set(ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
   set(BL_PRECISION "DOUBLE" CACHE INT "Precision of AMReX build")
   set(BL_USE_PARTICLES 1 CACHE INT "Include Particles classes in AMReX build")
   set(ENABLE_PROFILING 0 CACHE INT "Include profiling information in AMReX build")
   set(ENABLE_BACKTRACE 1 CACHE INT "Include backtrace information in AMReX build")
   set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)
   set(ENABLE_SUPERBUILD 1 CACHE INT  "Allow CMake to build and install AMReX")
      
endmacro(set_defaults)
