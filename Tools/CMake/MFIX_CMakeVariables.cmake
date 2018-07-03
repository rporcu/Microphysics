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

# The type of build ( will need to be uppercase )
set ( MFIX_BUILD_TYPE )

# Variable to hold amrex install path
set ( AMREX_INSTALL_PATH )

# The name of the core mfix library target
set ( MFIX_LIBNAME mfixcore )

# MFIX git variables
#set (MFIX_GIT_COMMIT)#
