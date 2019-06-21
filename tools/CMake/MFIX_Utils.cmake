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

  #
  # Get preprocessor flags
  #
  get_target_property( MFIX_DEFINES AMReX::amrex INTERFACE_COMPILE_DEFINITIONS )
  replace_genex( MFIX_DEFINES MFIX_Fortran_DEFINES LANGUAGE Fortran )
  replace_genex( MFIX_DEFINES MFIX_CXX_DEFINES LANGUAGE CXX )
  string(REPLACE " " ";-D" MFIX_Fortran_DEFINES "-D${MFIX_Fortran_DEFINES}")
  string(REPLACE " " ";-D" MFIX_CXX_DEFINES "-D${MFIX_CXX_DEFINES}")

  #
  # Get compiler flags flags
  #
  get_target_property( MFIX_FLAGS ${mfix_libname} INTERFACE_COMPILE_OPTIONS )
  replace_genex( MFIX_FLAGS MFIX_Fortran_FLAGS LANGUAGE Fortran )
  replace_genex( MFIX_FLAGS MFIX_CXX_FLAGS LANGUAGE CXX )

  if(NOT MFIX_Fortran_FLAGS)
    set(MFIX_Fortran_FLAGS ${CMAKE_Fortran_FLAGS})
  endif()

  if(NOT MFIX_CXX_FLAGS)
    set(MFIX_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  endif()

  string( REPLACE ";" " " MFIX_CXX_FLAGS "${MFIX_CXX_FLAGS}" )
  string( REPLACE ";" " " MFIX_Fortran_FLAGS "${MFIX_Fortran_FLAGS}" )

  #
  # Get extra includes
  #
  get_target_property( MFIX_INCLUDES AMReX::amrex INTERFACE_INCLUDE_DIRECTORIES )
  list( REMOVE_DUPLICATES MFIX_INCLUDES )

  #
  # Get extra libraries
  #
  get_target_property( TMP AMReX::amrex INTERFACE_LINK_LIBRARIES )
  replace_genex( TMP MFIX_LINK_LINE )
  string(REPLACE ";" " " MFIX_LINK_LINE "${MFIX_LINK_LINE}")

  #
  # Config summary
  #
  message( STATUS "MFIX configuration summary: " )
  message( STATUS "   C++ defines           = ${MFIX_CXX_DEFINES}" )
  message( STATUS "   Fortran defines       = ${MFIX_Fortran_DEFINES}" )
  message( STATUS "   C++ compiler          = ${CMAKE_CXX_COMPILER}" )
  message( STATUS "   Fortran compiler      = ${CMAKE_Fortran_COMPILER}" )
  message( STATUS "   C++ flags             = ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${MFIX_CXX_FLAGS}" )
  message( STATUS "   Fortran flags         = ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${MFIX_Fortran_FLAGS}" )
  message( STATUS "   MFIX includes         = ${MFIX_INCLUDES}" )
  message( STATUS "   MFIX extra link line  = ${MFIX_LINK_LINE}" )

endfunction()


#
#  Helper macro to replace genex in list
#
macro( replace_genex input_list output_list )

  cmake_parse_arguments( ARG "" "LANGUAGE" "" ${ARGN} )

  set(tmp_list ${${input_list}})

  # Replace all ; with a place holder (*)
  string( REPLACE ";" "*" tmp_list "${tmp_list}" )

  # Add tmp_list delimiter only where it suits us
  string( REPLACE ">*" ">;" tmp_list "${tmp_list}" )
  string( REPLACE "*$" ";$" tmp_list "${tmp_list}" )
  string( REPLACE "*/" ";/" tmp_list "${tmp_list}" )
  string( REPLACE "*" " "   tmp_list "${tmp_list}" )

  #
  # First remove entries related to:
  # 1) a compiler other than the one currently in use
  # 2) a build type other than the current one
  #
  foreach( item IN ITEMS ${tmp_list} )
    string( REPLACE "$<" "" item ${item} )
    string( REPLACE ">" "" item ${item} )
    string( REPLACE ":" "" item ${item} )

    # Accept build interface generator expressions
    string(REPLACE "BUILD_INTERFACE" "" item ${item})

    # Skip genex for compilers other than the one in use
    string( FIND ${item} "C_COMPILER_ID" idx1 )
    if( ${idx1} GREATER -1 )
      string( FIND ${item} "C_COMPILER_ID${CMAKE_C_COMPILER_ID}" idx2 )
      if( ${idx2} GREATER -1 )
        string( REPLACE "C_COMPILER_ID${CMAKE_C_COMPILER_ID}" "" item ${item} )
      else()
        continue()
      endif()
    endif()

    string( FIND ${item} "STREQUAL" idx1 )
    if( ${idx1} GREATER -1 )
      string( FIND ${item} "\"${CMAKE_Fortran_COMPILER_ID}\",\"${CMAKE_Fortran_COMPILER_ID}\"" idx2 )
      if( ${idx2} GREATER -1 )
        string( REPLACE "STREQUAL\"${CMAKE_Fortran_COMPILER_ID}\",\"${CMAKE_Fortran_COMPILER_ID}\"" "" item ${item} )
      else()
        continue()
      endif()
    endif()

    string( FIND ${item} "CXX_COMPILER_ID" idx1 )
    if( ${idx1} GREATER -1 )
      string( FIND ${item} "CXX_COMPILER_ID${CMAKE_CXX_COMPILER_ID}" idx2 )
      if( ${idx2} GREATER -1 )
        string( REPLACE "CXX_COMPILER_ID${CMAKE_CXX_COMPILER_ID}" "" item ${item} )
      else()
        continue()
      endif()
    endif()

    string( FIND ${item} "CONFIG" idx3 )
    if( ${idx3} GREATER -1 )
      string(FIND ${item} "${CMAKE_BUILD_TYPE}" idx4)
      if( ${idx4} GREATER -1 )
        string( REPLACE "CONFIG${CMAKE_BUILD_TYPE}" "" item ${item} )
      else()
        continue()
      endif()
    endif()

    # Extract by Language part
    if( ARG_LANGUAGE )
      string( FIND ${item} "COMPILE_LANGUAGE" idx1 )
      if(${idx1} GREATER -1)
        if( ${ARG_LANGUAGE} STREQUAL Fortran )
          string( FIND ${item} "Fortran" idx2 )
          if( ${idx2} GREATER -1 )
            string( REPLACE "COMPILE_LANGUAGEFortran" "" item ${item} )
          else()
            continue()
          endif()
        elseif(${ARG_LANGUAGE} STREQUAL CXX)
          string( FIND ${item} "CXX" idx2 )
          if( ${idx2} GREATER -1 )
            string( REPLACE "COMPILE_LANGUAGECXX" "" item ${item} )
          else()
            continue()
          endif()
        endif()
      endif()
    endif()

    # Now item should be ok to be added to final list
    list( APPEND ${output_list} ${item} )

  endforeach()

  if(${output_list})
    list( REMOVE_DUPLICATES ${output_list} )
  endif()

endmacro()


#
# Strip string from trailing and leading whitespace
# after veryfing it is not empty
#
macro(strip var)
  if(${var})
    string( STRIP "${${var}}" ${var} )
  endif()
endmacro()
