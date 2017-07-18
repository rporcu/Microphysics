###############################################

# Here we define all the global variable that #
# will be used in the CMake files

###############################################

#
# Check if project() has been called; abort if not
# 
if ( NOT PROJECT_NAME )
   message ( FATAL_ERROR "MFIX_CMakeVariables.cmake cannot be included\
before calling project()" )
endif ()

# 
# Check if this file has been loaded already
#
if (DEFINED __MFIX_CMAKEVARIABLES__)
   return ()
endif ()
set (__MFIX_CMAKEVARIABLES__ "")

# Set paths for build system
set ( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)

# The type of build ( will need to be uppercase )
set ( MFIX_BUILD_TYPE )

# Flags to enable/disable superbuild
set ( ENABLE_SUPERBUILD )

# Variable to hold amrex install path
set ( AMREX_INSTALL_PATH )

# Flags to accumulate preprocessor directives
set ( MFIX_DEFINES ) 

# Shorter-name variables for the compilers id
set ( FC_ID ${CMAKE_Fortran_COMPILER_ID} )
set ( CC_ID ${CMAKE_C_COMPILER_ID} )
set ( CXX_ID ${CMAKE_CXX_COMPILER_ID} )

# Compiler flags
set ( MFIX_Fortran_FLAGS )
set ( MFIX_C_FLAGS )
set ( MFIX_CXX_FLAGS )

# GNU compiler specific flags
set (MFIX_GNU_FFLAGS_DEBUG "-g -O0 -ggdb -fbounds-check -fbacktrace\
 -Wuninitialized -Wunused -finit-real=snan  -finit-integer=2147483647")
set (MFIX_GNU_FFLAGS_RELEASE "-O3 -DNDEBUG")
set (MFIX_GNU_FFLAGS_REQUIRED "-ffixed-line-length-none -ffree-line-length-none\
 -fno-range-check -fno-second-underscore")
set (MFIX_GNU_FFLAGS_FPE "-ffpe-trap=invalid,zero -ftrapv" )

set (MFIX_GNU_CXXFLAGS_DEBUG "-g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare")
set (MFIX_GNU_CXXFLAGS_RELEASE "-O3 -DNDEBUG")
set (MFIX_GNU_CXXFLAGS_REQUIRED "") #-ftemplate-depth-64 -Wno-deprecated")
set (MFIX_GNU_CXXFLAGS_FPE "-ftrapv")

# Intel compiler specific flags
set (MFIX_Intel_FFLAGS_DEBUG "-g -O0 -traceback -check bounds,uninit,pointers")
set (MFIX_Intel_FFLAGS_RELEASE "-O2 -ip -qopt-report=5 -qopt-report-phase=vec")
set (MFIX_Intel_FFLAGS_REQUIRED "-extend_source")
set (MFIX_Intel_FFLAGS_FPE "")

set (MFIX_Intel_CXXFLAGS_DEBUG "-g -O0 -traceback -Wcheck")
set (MFIX_Intel_CXXFLAGS_RELEASE "-O2 -ip -qopt-report=5 -qopt-report-phase=vec")
set (MFIX_Intel_CXXFLAGS_REQUIRED "")#-ftemplate-depth-64 -Wno-deprecated")
set (MFIX_Intel_CXXFLAGS_FPE "")

# PGI compiler specific flags
set (MFIX_PGI_FFLAGS_DEBUG "-O0 -Mbounds -Ktrap=divz,inv -Mchkptr")
set (MFIX_PGI_FFLAGS_RELEASE "-gopt -fast")
set (MFIX_PGI_FFLAGS_REQUIRED "-extend")
set (MFIX_PGI_FFLAGS_FPE "")

set (MFIX_PGI_CXXFLAGS_DEBUG "-O0 -Mbounds")
set (MFIX_PGI_CXXFLAGS_RELEASE "-gopt -fast")
set (MFIX_PGI_CXXFLAGS_REQUIRED "")#-ftemplate-depth-64 -Wno-deprecated")
set (MFIX_PGI_CXXFLAGS_FPE "")

# Cray compiler specific flags
set (MFIX_Cray_FFLAGS_DEBUG "-O0 -e -i")
set (MFIX_Cray_FFLAGS_RELEASE "-02")
set (MFIX_Cray_FFLAGS_REQUIRED "-extend")
set (MFIX_Cray_FFLAGS_FPE "")

set (MFIX_Cray_CXXFLAGS_DEBUG "-O0")
set (MFIX_Cray_CXXFLAGS_RELEASE "-02")
set (MFIX_Cray_CXXFLAGS_REQUIRED "")#-ftemplate-depth-64 -Wno-deprecated")
set (MFIX_Cray_CXXFLAGS_FPE "")

#
# Compile- and link-time variables 
#
set (MFIX_EXTRA_C_INCLUDE_PATH)
set (MFIX_EXTRA_CXX_INCLUDE_PATH)
set (MFIX_EXTRA_Fortran_INCLUDE_PATH)
set (MFIX_EXTRA_C_FLAGS)
set (MFIX_EXTRA_CXX_FLAGS)
set (MFIX_EXTRA_Fortran_FLAGS)
set (MFIX_EXTRA_C_LIBRARIES)
set (MFIX_EXTRA_CXX_LIBRARIES)
set (MFIX_EXTRA_Fortran_LIBRARIES)
set (MFIX_EXTRA_C_LIBS_DIR)
set (MFIX_EXTRA_CXX_LIBS_DIR)
set (MFIX_EXTRA_C_LINK_FLAGS)
set (MFIX_EXTRA_CXX_LINK_FLAGS)
set (MFIX_EXTRA_Fortran_LINK_FLAGS)
set (MFIX_EXTRA_C_LINK_LINE)
set (MFIX_EXTRA_CXX_LINK_LINE)
set (MFIX_EXTRA_Fortran_LINK_LINE)
set (MFIX_EXTRA_LINK_LINE)
set (MFIX_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/src
   ${PROJECT_SOURCE_DIR}/src/include)
set (MFIX_EXTRA_INCLUDE_PATH)
set (MFIX_EXTRA_LIBRARIES_PATH)


# AMReX Git variables
set (AMREX_GIT_REPO "https://github.com/AMReX-Codes/amrex.git" )
set (AMREX_GIT_COMMIT_MASTER  3506f5aea50d27237dda43df3ba4611fd4eda638 )
set (AMREX_GIT_COMMIT_DEVELOP b4d550904abff4bd269963681f39e5af7d12be32 )
set (AMREX_GIT_TAG)  # The commit id or branch to download 

# AMReX Superbuild variables
set (AMREX_SUPERBUILD_DIR   ${PROJECT_BINARY_DIR}/ThirdParty)

# MFIX git variables
set (MFIX_GIT_COMMIT)
set (MFIX_GIT_BRANCH)
