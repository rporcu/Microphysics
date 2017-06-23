###################################################################
# Check if dir or file given by path exists and issue a warning or
#  error if not
##################################################################
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

