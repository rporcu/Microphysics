#
# This file adds a target to typecheck C++ calls to Fortran routines.
# This works only with GNU compiler so it is not executed if the
# compiler id is not GNU.
# This file will check for Fortran/C++ headers in the dir of inclusion,
# and in all dirs below.
# F90SRC is the list of Fortran 90 sources
# CXXINCLUDES is the list of CXX headers
# Both must be defined before calling this file. 
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
# The following variables needs to be defined before including
# this file
#
if ( (NOT F90SRC) OR (NOT CXXINCLUDES) )
   message (FATAL_ERROR "F90SRC and CXXINCLUDES must be defined before \
including MFIX_Typecheck.cmake ")
endif ()

#
# Directory for typecheck 
#
set ( TYPECHECK_DIR  ${CMAKE_BINARY_DIR}/TypeCheckTemp )


#
# Object library to generate before typecheck
# 
add_library ( typecheckobjs OBJECT EXCLUDE_FROM_ALL ${F90SRC} )
set_target_properties ( typecheckobjs
   PROPERTIES
   Fortran_MODULE_DIRECTORY ${TYPECHECK_DIR} ) 


   
if ( ENABLE_SUPERBUILD )
   # AMReX must be installed before performing typecheck
   add_dependencies ( typecheckobjs amrex )
endif ()


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
# Find defines needed for typecheck
#
set (DEFINES)
string (STRIP ${MFIX_DEFINES} DEFINES) 
string (REPLACE " " ";" DEFINES ${DEFINES})

#
# Find C headers needed for type check
#

# Filter out %_f.H and %_F.H from list of includes
# Note that this can be achieved natively in CMake >= 3.6
# by using " list (FILTER ...) "
foreach ( item ${CXXINCLUDES}} )
   string ( REGEX MATCH "_f.H" TMP1 ${item})
   string ( REGEX MATCH "_F.H" TMP2 ${item})
   if ( (NOT TMP1) AND (NOT TMP2) )
      list (REMOVE_ITEM CXXINCLUDES ${item})	
   endif ()
endforeach ()

#
# Generate cppd files
#
if (ENABLE_DP)
   set (AMREX_REAL double)
else ()
   set (AMREX_REAL  float)
endif ()

if (ENABLE_DP_PARTICLES)
   set (AMREX_PARTICLE_REAL double)
else ()
   set (AMREX_PARTICLE_REAL  float)
endif ()

set (CXXHEADERS)
foreach ( file ${CXXINCLUDES} )
   get_filename_component ( fname ${file} NAME ) # This strips away the path
   set ( CPPD_FILE ${TYPECHECK_DIR}/${fname}-cppd.h )
   add_custom_command ( OUTPUT  ${CPPD_FILE} COMMAND ${CMAKE_C_COMPILER}
      ARGS ${DEFINES} ${INCLUDES} -E -P -x c -std=c99 ${file} > ${CPPD_FILE}
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
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating ${CPPD_FILE} " )
   list (APPEND CXXHEADERS ${CPPD_FILE})
endforeach ()


#
# Generate Origin Files
#
set (F90ORIG)
foreach ( file ${F90SRC} )
   get_filename_component ( fname ${file} NAME ) # This strips away the path
   set ( ORIG_FILE ${TYPECHECK_DIR}/${fname}.orig )
   add_custom_command ( OUTPUT  ${ORIG_FILE} COMMAND ${CMAKE_Fortran_COMPILER}
      ARGS ${DEFINES} ${INCLUDES} -fsyntax-only -fdump-fortran-original ${file} > ${ORIG_FILE}
      DEPENDS ${file} typecheckobjs
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating ${ORIG_FILE} " )
   list (APPEND F90ORIG ${ORIG_FILE}) 
endforeach ()

# 
# Add typecheck target
#
add_custom_target ( typecheck
   COMMAND python3  ${AMREX_TYPECHECKER}
   --workdir ${TYPECHECK_DIR} --output ${TYPECHECK_DIR}/amrex_typecheck.ou
   DEPENDS ${F90ORIG} ${CXXHEADERS}
   WORKING_DIRECTORY ${TYPECHECK_DIR}
   COMMENT "Type-checking") 

#
# Undefine variables so this file may be included elsewhere
# ( this must be tested )
#
set (F90SRC "")
set (CXXINCLUDES "")
