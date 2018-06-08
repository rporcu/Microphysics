#
# This file adds a target to typecheck C++ calls to Fortran routines.
# This works only with GNU compiler so it is not executed if the
# compiler id is not GNU.
# This file will check for Fortran/C++ headers in the dir of inclusion,
# and in all dirs below.
# 
if (NOT (CMAKE_Fortran_COMPILER_ID MATCHES GNU))
   return ()
endif ()

#
# Check that MFIX_Config.cmake has been loaded already
# 
if (NOT (DEFINED __MFIX_CONFIG__))
   message (FATAL "MFIX_Config.cmake must be included before \
MFIX_Typecheck.cmake")
endif ()


#
# Set the name of the directory from where the file is called
#
message ("Typecheck in ${CMAKE_CURRENT_SOURCE_DIR}")


#
# Directory for typecheck 
#
set ( TYPECHECK_DIR  ${CMAKE_BINARY_DIR}/TypeCheckTemp )

#
# Get the fortran sources and the fortran headers for MFIXCORE
#
get_target_property ( ALLSRC ${MFIXCORE} SOURCES )

set ( F90SRC )
set ( F90HEADERS )
set ( HEXT ".H" )
set ( FEXT ".f90;.F90;.f;.F")

foreach ( item ${ALLSRC} )
  
  # Get the file extension
  get_filename_component ( FILETYPE ${item} EXT )

  # Add to F90 sources
  if ( FILETYPE IN_LIST FEXT )
    list ( APPEND F90SRC ${item} )
  endif()

  # Add to F90 Headers
  if ( FILETYPE IN_LIST HEXT)
    string ( REGEX MATCH "_f.H" COND1 ${item})
    string ( REGEX MATCH "_F.H" COND2 ${item})

    if ( COND1 OR COND2 )
      print(item)
      list ( APPEND F90HEADERS ${item})	
    endif ()

  endif ()
  
endforeach ()

#
# Object library to generate before typecheck
# 
add_library ( typecheckobjs OBJECT EXCLUDE_FROM_ALL ${F90SRC} )
set_target_properties ( typecheckobjs
   PROPERTIES
   Fortran_MODULE_DIRECTORY ${TYPECHECK_DIR} ) 

#
# Find includes needed for typecheck
#
set (INCLUDES)
get_target_property ( TMP typecheckobjs INCLUDE_DIRECTORIES )
list (APPEND TMP ${TYPECHECK_DIR})

foreach (item ${TMP})
   list ( APPEND INCLUDES -I${item} )
endforeach ()

#
# Generate cppd files
#
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

set (CPPDHEADERS)
foreach ( file ${F90HEADERS} )
   get_filename_component ( fname ${file} NAME ) # This strips away the path
   set ( CPPD_FILE ${fname}-cppd.h )
   get_filename_component ( fullname ${file} ABSOLUTE ) # This add the absolute path to fname
   add_custom_command ( OUTPUT  ${CPPD_FILE} COMMAND ${CMAKE_C_COMPILER}
      ARGS ${MFIX_DEFINES} ${INCLUDES} -E -P -x c -std=c99 ${fullname} > ${CPPD_FILE}
      COMMAND sed
      ARGS -i -e 's/amrex::Real/${AMREX_REAL}/g' ${CPPD_FILE} 
      COMMAND sed
      ARGS -i -e 's/amrex_real/${AMREX_REAL}/g' ${CPPD_FILE} 
      COMMAND sed
      ARGS -i -e 's/amrex_particle_real/${AMREX_PARTICLE_REAL}/g' ${CPPD_FILE}
      COMMAND sed
      ARGS -i -e '/typedef\\s*${AMREX_REAL}/d' ${CPPD_FILE} 
      COMMAND sed
      ARGS -i -e 's/\\&/*/g' ${CPPD_FILE} 
      DEPENDS ${file} typecheckobjs # Leave dependency to typecheck so typecheckdir is created
      WORKING_DIRECTORY ${TYPECHECK_DIR}  
      COMMENT "Generating ${CPPD_FILE} " )
   list (APPEND CPPDHEADERS ${CPPD_FILE})
endforeach ()

#
# Generate Origin Files
#
set (F90ORIG)
foreach ( file ${F90SRC} )
   get_filename_component ( fname ${file} NAME ) # This strips away the path
   set ( ORIG_FILE ${fname}.orig )
   get_filename_component ( fullname ${file} ABSOLUTE ) # This add the absolute path to fname
   add_custom_command ( OUTPUT  ${ORIG_FILE} COMMAND ${CMAKE_Fortran_COMPILER}
      ARGS ${DEFINES} ${INCLUDES} -fsyntax-only -fdump-fortran-original ${fullname} > ${ORIG_FILE}
      DEPENDS ${file} typecheckobjs
      WORKING_DIRECTORY ${TYPECHECK_DIR} 
      COMMENT "Generating ${ORIG_FILE} " )
   list (APPEND F90ORIG ${ORIG_FILE}) 
endforeach ()

# 
# Add typecheck target
#
add_custom_target ( typecheck
   COMMAND python3  ${AMREX_TYPECHECKER}
   --workdir ${TYPECHECK_DIR} --output ${TYPECHECK_DIR}/amrex_typecheck.ou
   DEPENDS ${F90ORIG} ${CPPDHEADERS}
   WORKING_DIRECTORY ${TYPECHECK_DIR}
   COMMENT "Type-checking") 

#
# Undefine variables so this file may be included elsewhere
# ( this must be tested )
#
set (F90SRC "")

