#
# This file provides the following variables
#
#  MFIX_FFLAGS_DEBUG 
#  MFIX_FFLAGS_RELEASE
#  MFIX_FFLAGS_REQUIRED
#  MFIX_FFLAGS_FPE
#  MFIX_CXXFLAGS_DEBUG
#  MFIX_CXXFLAGS_RELEASE
#  MFIX_CXXFLAGS_REQUIRED
#  MFIX_CXXFLAGS_FPE
#

#
# Check wether the compiler ID has been defined
# 
if (  NOT (DEFINED CMAKE_Fortran_COMPILER_ID) OR
      NOT (DEFINED CMAKE_C_COMPILER_ID) OR 
      NOT (DEFINED CMAKE_CXX_COMPILER_ID) )
   message ( FATAL_ERROR "Compiler ID is UNDEFINED" )
endif ()

#
# Check the same compiler suite is used for all languages
# 
set ( COMPILER ${CMAKE_C_COMPILER_ID} )

if (  NOT (${COMPILER} STREQUAL ${CMAKE_Fortran_COMPILER_ID}) OR
      NOT (${COMPILER} STREQUAL ${CMAKE_CXX_COMPILER_ID}) )
   message ( FATAL_ERROR "C compiler ID does not match Fortran/C++ compiler ID" )
endif ()

# 
# Select Compiler flags based on compiler ID
# 
if ( ${COMPILER} STREQUAL "GNU" )

   # GNU compiler specific flags
   set (MFIX_FFLAGS_DEBUG "-g -O0 -ggdb -fbounds-check -fbacktrace\
 -Wuninitialized -Wunused -finit-real=snan  -finit-integer=2147483647")
   set (MFIX_FFLAGS_RELEASE "-O3 -DNDEBUG")
   set (MFIX_FFLAGS_REQUIRED "-ffixed-line-length-none -ffree-line-length-none\
 -fno-range-check -fno-second-underscore")
   set (MFIX_FFLAGS_FPE "-ffpe-trap=invalid,zero -ftrapv" )

   set (MFIX_CXXFLAGS_DEBUG "-g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare")
   set (MFIX_CXXFLAGS_RELEASE "-O3 -DNDEBUG")
   set (MFIX_CXXFLAGS_REQUIRED "") #-ftemplate-depth-64 -Wno-deprecated")
   set (MFIX_CXXFLAGS_FPE "-ftrapv")

elseif ( ${COMPILER} STREQUAL "Intel" )
   
   # Intel compiler specific flags
   set (MFIX_FFLAGS_DEBUG "-g -O0 -traceback -check bounds,uninit,pointers")
   set (MFIX_FFLAGS_RELEASE "-O3 -ip -qopt-report=5 -qopt-report-phase=vec")
   set (MFIX_FFLAGS_REQUIRED "-extend_source")
   set (MFIX_FFLAGS_FPE "")
   
   set (MFIX_CXXFLAGS_DEBUG "-g -O0 -traceback -Wcheck")
   set (MFIX_CXXFLAGS_RELEASE "-O3 -ip -qopt-report=5 -qopt-report-phase=vec")
   set (MFIX_CXXFLAGS_REQUIRED  "-std=c++11")#-ftemplate-depth-64 -Wno-deprecated")
   set (MFIX_CXXFLAGS_FPE "")

elseif (${COMPILER} STREQUAL "PGI")

   # PGI compiler specific flags
   set (MFIX_FFLAGS_DEBUG "-g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr")
   set (MFIX_FFLAGS_RELEASE "-gopt -fast")
   set (MFIX_FFLAGS_REQUIRED "-extend")
   set (MFIX_FFLAGS_FPE "")
   
   set (MFIX_CXXFLAGS_DEBUG "-O0 -Mbounds")
   set (MFIX_CXXFLAGS_RELEASE "-gopt -fast")
   set (MFIX_CXXFLAGS_REQUIRED "")#-ftemplate-depth-64 -Wno-deprecated")
   set (MFIX_CXXFLAGS_FPE "")

elseif ( ${COMPILER} STREQUAL "Cray" )

   # Cray compiler specific flags
   set (MFIX_FFLAGS_DEBUG "-g -O0 -e i")
   set (MFIX_FFLAGS_RELEASE "-O2")
   set (MFIX_FFLAGS_REQUIRED "-N 255 -h list=a")
   set (MFIX_FFLAGS_FPE "")
   
   set (MFIX_CXXFLAGS_DEBUG "-g -O0")
   set (MFIX_CXXFLAGS_RELEASE "-O2")
   set (MFIX_CXXFLAGS_REQUIRED "-h std=c++11 -h list=a")#-ftemplate-depth-64 -Wno-deprecated")
   set (MFIX_CXXFLAGS_FPE "")

elseif ()

   message ( FATAL_ERROR "Compiler NOT recognized: ID is ${COMPILER}" )

endif ()



