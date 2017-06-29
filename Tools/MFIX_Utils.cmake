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
# Print variable (useful for debug) 
#
function (print var)
   message (STATUS "   ${var} = ${${var}}")
endfunction ()


#
# Print option 
#
function (print_option name value)
   message (STATUS "   ${name} = ${value}")
endfunction ()

#
# Find git info
#
function ( get_git_info branch commit )

 
   # Find branch
   execute_process (
      COMMAND git branch
      COMMAND grep \*
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} 
      OUTPUT_VARIABLE out
      ERROR_VARIABLE  err
      )

   if (err)
      message (WARNING "Failing to retrieve MFIX Git branch")
   else ()
      string ( REPLACE "*" "" out ${out} )
      string ( STRIP ${out} out)
      message (STATUS "MFIX branch: ${out}" )
      set ( ${branch} ${out} PARENT_SCOPE)
   endif ()

   
   # Find commit
   execute_process (
      COMMAND git rev-parse HEAD
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} 
      OUTPUT_VARIABLE out
      ERROR_VARIABLE  err
      )

   if (err)
      message (WARNING "Failing to retrieve MFIX Git commit")
   else ()
      string (STRIP ${out} out)
      message (STATUS "MFIX commit: ${out}" )
      set ( ${commit} ${out} PARENT_SCOPE)
   endif ()
  
endfunction ()


