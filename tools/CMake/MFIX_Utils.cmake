#
# Stop if in-source build
#
macro( check_build_tree_path )
  string( FIND "${CMAKE_BINARY_DIR}" "${PROJECT_BINARY_DIR}/src" RESULT )
  if( NOT "${RESULT}" STREQUAL "-1" )
    message( FATAL_ERROR  "ERROR: in-source builds are not allowed!" )
  endif()
endmacro()


#
# Print variable (useful for debug)
#
function(print var)
  message(STATUS "   ${var} = ${${var}}")
endfunction()


#
# Print option
#
function(print_option name value)
  message(STATUS "   ${name} = ${value}")
endfunction()

#
# Find git info
#
macro( get_git_info ) # EXTRA ARGS: branch commit

  # Find branch
  execute_process(
    COMMAND git branch
    COMMAND grep \*
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE  err
    )

  if(err)
    message(WARNING "Failing to retrieve MFIX Git branch")
  else()
    string( REPLACE "*" "" out ${out} )
    string( STRIP ${out} out )
    message( STATUS "MFIX branch: ${out}" )
    if( ${ARGC} GREATER 0 ) # branch
      set( ${ARGV0} ${out} )
    endif()
  endif()

  # Find commit
  execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE  err
    )

  if(err)
    message(WARNING "Failing to retrieve MFIX Git commit")
  else()
    string(STRIP ${out} out)
    message( STATUS "MFIX commit: ${out}" )
    if( ${ARGC} EQUAL 2 ) # commit
      set( ${ARGV1} ${out} )
    endif()
  endif()

  unset(out)

endmacro()


#
# This sets CMake_<LANG>_FLAGS_<CONFIG> to default values
# if the variable is empty
#
macro( set_default_config_flags )

  if( NOT CMAKE_Fortran_FLAGS_DEBUG )
    set(CMAKE_Fortran_FLAGS_DEBUG "-g")
  endif()

  if( NOT CMAKE_Fortran_FLAGS_RELEASE )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  endif()

  if( NOT CMAKE_CXX_FLAGS_DEBUG )
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
  endif()

  if( NOT CMAKE_CXX_FLAGS_RELEASE )
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
  endif()

endmacro()


#
# Print Configuration Summary
#
function( print_mfix_configuration_summary mfix_libname ) 
   
  if(NOT TARGET ${mfix_libname})
    message(AUTHOR_WARNING "Target ${mfix_libname} is not defined.")
    return()
  endif()

  if(NOT TARGET AMReX::amrex)
    message(AUTHOR_WARNING "Target ${mfix_libname} is not defined.")
    return()
  endif()

  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPERCASE)


  # Get AMReX cmake functions 
  include(AMReXGenexHelpers)
  include(AMReXTargetHelpers)

  get_target_properties_flattened(${mfix_libname} _includes _defines _flags _link_line)

  set(_lang CXX Fortran) 
  set(_prop _includes _defines _flags _link_line )


  # Loop over all combinations of language and property and extract 
  # what you need 
  foreach( _p IN LISTS _prop )
     foreach( _l IN LISTS _lang )

        string(TOLOWER ${_l} _ll) # Lower case language name

        # _${_ll}${_p} is a variable named as _lang_property,
        # both lower case. 
        evaluate_genex(${_p} _${_ll}${_p}
           LANG ${_l}
           COMP ${CMAKE_${_l}_COMPILER_ID}
           CONFIG ${CMAKE_BUILD_TYPE}
           INTERFACE BUILD)

        if (_${_ll}${_p})

           list(REMOVE_DUPLICATES _${_ll}${_p})
           
           if ("${_p}" STREQUAL "_defines")
              string(REPLACE ";" " -D" _${_ll}${_p} "-D${_${_ll}${_p}}")
           elseif ("${_p}" STREQUAL "_includes")
              string(REPLACE ";" " -I" _${_ll}${_p} "-I${_${_ll}${_p}}")
           else()
              string(REPLACE ";" " " _${_ll}${_p} "${_${_ll}${_p}}")
           endif ()              

        endif ()
        
     endforeach()
  endforeach ()

  
  #
  # Config summary
  #
  message( STATUS "MFIX configuration summary: " )
  message( STATUS "   C++ defines           = ${_cxx_defines}" )
  message( STATUS "   Fortran defines       = ${_fortran_defines}" )
  message( STATUS "   C++ compiler          = ${CMAKE_CXX_COMPILER}" )
  message( STATUS "   Fortran compiler      = ${CMAKE_Fortran_COMPILER}" )
  message( STATUS "   C++ flags             = ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${CMAKE_CXX_FLAGS} ${_cxx_flags}" )
  message( STATUS "   Fortran flags         = ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${CMAKE_Fortran_FLAGS} ${_fortran_flags}" )
  if (ENABLE_CUDA)
     message( STATUS "   CUDA flags            = ${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${CMAKE_CUDA_FLAGS}" )
  endif ()
  message( STATUS "   MFIX includes         = ${_cxx_includes}" )
  message( STATUS "   MFIX link line        = ${_cxx_link_line}" )

endfunction()



