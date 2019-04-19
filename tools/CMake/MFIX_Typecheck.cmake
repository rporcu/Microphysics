# 
# This module  provides function add_typecheck_target().
# 
# add_typecheck_target() adds a target to check C++ signatures of Fortran routines are
# consistet with the Fortran interface of the routine itself.
# Works only with GNU compiler so it returns if the compiler id is not GNU.
#
function( add_typecheck_target _target)

   # 
   # Check if we have all we need to define the typecheck target
   # 
   if ( NOT (CMAKE_Fortran_COMPILER_ID MATCHES GNU) OR NOT (CMAKE_C_COMPILER_ID MATCHES GNU) )
      message(WARNING "Typecheck disabled because compiler ID is not GNU")
      return ()
   endif ()

   if (NOT (TARGET AMReX::amrex) )
      message(FATAL_ERROR "Target AMReX::amrex does not exist")
   endif ()

   if ( NOT (TARGET ${_target}) )
      message(AUTHOR_WARNING "Skipping type-checking on target ${_target} bacause it does not exist")
      return ()
   endif ()

   find_package(Python3 COMPONENTS Interpreter Development QUIET)
   if (NOT Python3_FOUND)
      message(WARNING "Typecheck disabled because Python 3 was not found")
      return ()     
   endif ()

   # 
   # Set directory for typecheck 
   # 
   set( TYPECHECK_DIR  "${CMAKE_CURRENT_BINARY_DIR}/TypeCheckTemp/${_target}" )

   # 
   # Get the fortran sources and the fortran-interfaces headers from _target   
   # 
   get_target_property( _sources ${_target} SOURCES )

   set(_fsources ${_sources})
   list(FILTER _fsources INCLUDE REGEX "(\.f90|\.f|\.F90|\.F)$" )
   
   set(_fheaders ${_sources})
   list(FILTER _fheaders INCLUDE REGEX "(\_f\.H|\_F\.H)$")

   # 
   # Find includes and defines required for those sources
   # Must be done manually since we will use them in a custom command
   #
   set(_includes ${TYPECHECK_DIR})

   get_target_property( _target_includes ${_target}    INCLUDE_DIRECTORIES )
   get_target_property( _amrex_includes  AMReX::amrex  INTERFACE_INCLUDE_DIRECTORIES )

   foreach (_item IN LISTS _target_includes _amrex_includes)
      if (_item)
         list(APPEND _includes  ${_item})
      endif ()
   endforeach ()
   
   list( REMOVE_DUPLICATES _includes )

   set(_defines)

   get_target_property( _target_defines ${_target}    COMPILE_DEFINITIONS )
   get_target_property( _amrex_defines  AMReX::amrex  INTERFACE_COMPILE_DEFINITIONS )

   foreach (_item IN LISTS _target_defines _amrex_defines)
      if (_item)
         list(APPEND _defines  ${_item})
      endif ()
   endforeach ()
   
   list( REMOVE_DUPLICATES _defines )


   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 1: create fortran modules
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   #
   # Define the typecheckobjs library.
   # Basically, we only want to generate the module files
   # associated with fortran otherwise the other two steps
   # of the type check (see below) will fail because
   # it will not be able to include symbols
   # if modules are not there.
   #

   # ---------------------------------------------------------------------------------
   # NOTE:
   # ---------------------------------------------------------------------------------
   # We could skip generating the objs file and create only the modules
   # by setting the following compile option:
   # target_compile_options(typecheckobjs PRIVATE -fsyntax-only)
   # However, since this would only generate ".mod" files and no ".o" files, the target
   # would always be out of date and repeat itself.
   # A work around would be to create the module and the orig file at the same time.
   # this could be achieved with a add_custom_command which is aware of dependecies info.
   # To this end, we could use IMPLICIT_DEPENDS. Unfortunately this option is supported
   # for C and C++ only for now.

   set(_typecheckobjs  typecheckobjs_${_target})
   
   add_library( ${_typecheckobjs} OBJECT EXCLUDE_FROM_ALL )

   target_sources( ${_typecheckobjs} 
      PRIVATE
      ${_fsources}
      )

   set_target_properties( ${_typecheckobjs} 
      PROPERTIES
      Fortran_MODULE_DIRECTORY  ${TYPECHECK_DIR}
      )

   target_include_directories( ${_typecheckobjs} 
      PUBLIC
      ${_target_includes} ${TYPECHECK_DIR}
      )
   
   target_link_libraries( ${_typecheckobjs} 
      PRIVATE
      AMReX::amrex
      )

   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 2: create CPPD files from C++ headers
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   if (AMREX_ENABLE_DP)
      set (AMREX_REAL double)
   else ()
      set (AMREX_REAL  float)
   endif ()

   if (AMREX_ENABLE_DP_PARTICLES)
      set (AMREX_PARTICLE_REAL double)
   else ()
      set (AMREX_PARTICLE_REAL  float)
   endif ()

   #
   # Manually setup includes and defines
   # 
   if (_includes)
      string(REPLACE ";" ";-I" _includes "-I${_includes}")
   endif ()

   # get rid of genex in define list
   string( GENEX_STRIP "${_defines}" _defines )
   if (_defines)
      string(REPLACE ";" ";-D" _defines "-D${_defines}")
   endif ()
   
   set (_cppd)
   foreach ( _file IN LISTS _fheaders )
      get_filename_component( _fname    ${_file} NAME )     # This strips away the path
      get_filename_component( _fullname ${_file} ABSOLUTE ) # This add the absolute path to fname
      set( _cppd_file ${TYPECHECK_DIR}/${_fname}-cppd.h )
      add_custom_command(
         OUTPUT  ${_cppd_file}
         COMMAND ${CMAKE_C_COMPILER}
	 ARGS    ${_defines} ${_includes} -E -P -x c -std=c99 ${_fullname} > ${_cppd_file}
	 COMMAND sed
	 ARGS -i -e 's/amrex::Real/${AMREX_REAL}/g' ${_cppd_file}
	 COMMAND sed
	 ARGS -i -e 's/amrex_real/${AMREX_REAL}/g' ${_cppd_file}
	 COMMAND sed
	 ARGS -i -e 's/amrex_particle_real/${AMREX_PARTICLE_REAL}/g' ${_cppd_file}
	 COMMAND sed
	 ARGS -i -e '/typedef\\s*${AMREX_REAL}/d' ${_cppd_file}
         COMMAND sed
         ARGS -i -e 's/AMREX_GPU_DEVICE/ /g' ${_cppd_file}
         COMMAND sed
         ARGS -i -e 's/AMREX_GPU_HOST_DEVICE/ /g' ${_cppd_file}
	 COMMAND sed
	 ARGS -i -e 's/\\&/*/g' ${_cppd_file}
         DEPENDS ${_file} 
	 WORKING_DIRECTORY ${TYPECHECK_DIR}
	 COMMENT "Generating ${_fname}-cppd.h" )
      list(APPEND _cppd ${_cppd_file})
   endforeach ()

   
   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 3: generate orig files from fortran sources
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   set(_orig)
   foreach ( _file IN LISTS _fsources )
      get_filename_component( _fname    ${_file} NAME )     # This strips away the path
      get_filename_component( _fullname ${_file} ABSOLUTE ) # This add the absolute path to fname
      set( _orig_file ${TYPECHECK_DIR}/${_fname}.orig )     # Use fullpath otherwise it rebuilds everytime 
      add_custom_command(
         OUTPUT   ${_orig_file}
         COMMAND  ${CMAKE_Fortran_COMPILER}
         ARGS     ${_defines} ${_includes} -fsyntax-only -fdump-fortran-original ${_fullname} > ${_orig_file}
         DEPENDS  ${_file} ${_typecheckobjs}
         WORKING_DIRECTORY    ${TYPECHECK_DIR} 
         COMMENT  "Generating ${_fname}.orig" )
      list(APPEND _orig ${_orig_file}) 
   endforeach ()


   # 
   # Add typecheck target
   #
   set(_outfile  "${TYPECHECK_DIR}/${_target}_typecheck.ou" )
   
   add_custom_target( typecheck_${_target}
      COMMAND python3  ${AMREX_TYPECHECKER}
      --workdir ${TYPECHECK_DIR} --output ${_outfile}
      DEPENDS ${_cppd} ${_orig}
      WORKING_DIRECTORY ${TYPECHECK_DIR}
      )

   # Add rules to remove typecheck outfile when cleaning 
   set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${_outfile})
   
endfunction ()
