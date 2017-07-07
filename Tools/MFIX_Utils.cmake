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
      set ( tmp  ${${all_var}} ${${new_var}} )

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
      #set ( ${all_defines} ${${all_defines}} -D${new_define} PARENT_SCOPE )
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


#
# Create target for typecheck
# ATTENTION: only works with GNU compiler
#
function ( create_typecheck_target f90src cxxheaders typecheck_exe )

   # If not GNU compiler, return
   if ( NOT (CMAKE_Fortran_COMPILER_ID MATCHES GNU) )
      return ()
   endif ()
   if ( NOT (CMAKE_CXX_COMPILER_ID MATCHES GNU) )
      return ()
   endif ()

   # Filter out %_f.H and %_F.H from list of includes
   # Note that this can be achieved natively in CMake >= 3.6
   # by using " list (FILTER ...) "
   foreach ( item ${${cxxheaders}})
      string ( REGEX MATCH "_f.H" tmp1 ${item})
      string ( REGEX MATCH "_F.H" tmp2 ${item})
      if ( (NOT tmp1) AND (NOT tmp2) )
	 list (REMOVE_ITEM ${cxxheaders} ${item})	
      endif ()
   endforeach ()

   # Lists of fortran origin files
   string (REPLACE ".f90" ".f90.orig" f90orig "${${f90src}}" )
   string (REPLACE ";" " " f90orig "${f90orig}")
   string (REPLACE ";" " " f90src "${${f90src}}")

   # Generate include line
   get_property (dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES )
   set (include_line ) #"" )
   foreach (dir ${dirs})
      set ( dir  -I${dir} )
      list (APPEND include_line ${dir} )
   endforeach ()

   # Generate list of defines
   string (STRIP ${MFIX_DEFINES} defines) 
   string (REPLACE " " ";" defines ${defines}) 


   print (defines)
   # Create typecheck dir
   set ( tempdir ${CMAKE_BINARY_DIR}/TypeCheckTemp )

   add_custom_command ( OUTPUT ${tempdir}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${tempdir}
      COMMENT "Making typecheck temporary directory" )

   # add_custom_command ( OUTPUT  ${tempdir}
   #    COMMAND mkdir ARGS ${tempdir}
   #    COMMENT "Making typecheck temporary directory" )
   
   print (MFIX_DEFINES)
   
   # Preprocess fortran files
   set ( cmd gfortran ${defines} ${include_line} -fsyntax-only -fdump-fortran-original )
   print (cmd)
   add_custom_command(
      OUTPUT adjust_a_g.f90.orig 
      COMMAND ${cmd} 
      ARGS adjust_a_g.f90 > ${tempdir}/adjust_a_g.f90.orig
      DEPENDS adjust_a_g.f90 ${tempdir} 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Preprocessing ${f}..."
      # VERBATIM
      # USES_TERMINAL
      )
   


   
  # add_custom_command ( OUTPUT  adjust_a_g.f90.orig  #${f90orig}
  #    COMMAND gfortran
  #    ARGS ${MFIX_DEFINES} ${include_line} -fsyntax-only -fdump-fortran-original $< > $@
  #    #[COMMAND command2 [ARGS] [args2...] ...]
  #    MAIN_DEPENDENCY   adjust_a_g.f90 #${f90src} 
  #    DEPENDS ${tempdir}
  #    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} #${CMAKE_BINARY_DIR}/TypeCheckTemp
  #    COMMENT "Generating origin files " )

   add_custom_target ( typecheck DEPENDS adjust_a_g.f90.orig ) #${f90orig} )

   
   # add_custom_target ( [ALL] [command1 [args1...]]
   #    [COMMAND command2 [args2...] ...]
   #    [DEPENDS depend depend depend ... ]
   #    [BYPRODUCTS [files...]]
   #    [WORKING_DIRECTORY dir]
   #    [COMMENT comment]
   #    [VERBATIM] [USES_TERMINAL]
   #                      [SOURCES src1 [src2...]])



   
endfunction  ()


