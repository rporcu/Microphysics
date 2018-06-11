cmake_minimum_required (VERSION 3.5)

#
# Set the search path for cmake modules
#
set ( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Tools/CMake )

#
# Here we set the only options that decide if superbuild
# is enabled or not
# 
set ( AMREX_INSTALL_DIR "" CACHE PATH
   "Path to installation directory (leave empty for superbuild)")

set ( SUPERBUILD  ON )

if ( AMREX_INSTALL_DIR )
   set ( SUPERBUILD OFF )
endif ()

#
# Now decide if invoking superbuild. If so, this listfile will be
# re-read by the MFIX_Superbuild.cmake 
# 
if ( SUPERBUILD )
   include ( MFIX_Superbuild )
   return ()
endif ()


#
# This is where the MFIX-standalone config starts
# 
project (MFIX-Exa)

#
# Define the languages used by the project
#
enable_language (C)
enable_language (CXX)
enable_language (Fortran)

#
# Require C++11 standard
#
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON) 
set (CMAKE_CXX_EXTENSIONS OFF)

#
# Always verbose output during make
# 
# set ( CMAKE_VERBOSE_MAKEFILE on )


#
# Include module with custom macros and functions
#
include ( MFIX_Utils )

get_git_info ( var1 var2 )

#
# Stop if building in source
# 
check_build_tree_path ()

#
# Initialize all the global variables
#
include ( MFIX_CMakeVariables )

# 
# Load MFIX options
# 
include ( MFIX_Options )

#
# Config MFIX build
# 
include ( MFIX_Config )


# This is the main file
set ( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp )

# Add the external include paths and link directories
link_directories (${MFIX_EXTRA_LIBRARIES_PATH})

# Main
add_subdirectory (src)

# We should check what this line does
file (GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
add_executable (mfix ${MAIN} ${USR_OVERRIDES} )
target_link_libraries ( mfix ${MFIX_LIBNAME} ${MFIX_EXTRA_LINK_LINE} ) 


# Configure type-checking target
include ( MFIX_Typecheck )

# Add rule to generate TAGS
# on macOS ctags-exuberant is just ctags
find_program(CTAGS_EXU "ctags-exuberant")
find_program(CTAGS_CMD "ctags")
if(CTAGS_EXU)
    add_custom_target ( tags
        COMMAND ctags-exuberant -R    --fortran-kinds=+i  --exclude="build/amrex/installdir/*" . 
        COMMAND ctags-exuberant -R -e --fortran-kinds=+i  --exclude="build/amrex/installdir/*" . 
        COMMENT "Generating tags files"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    # Some file systems are case-insensitive so the above will fail:
    #  the first command would overwrite the second!
    add_custom_target ( ctags
        COMMAND ctags-exuberant -R    --fortran-kinds=+i --exclude="build/amrex/installdir/*" .
        COMMENT "Generating only ctags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    add_custom_target( etags
        COMMAND ctags-exuberant -R -e --fortran-kinds=+i --exclude="build/amrex/installdir/*" .
        COMMENT "Generating only etags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    # Avoid "collisions" with AMReX => focus entirely on MFiX sources:
    add_custom_target ( mfix_ctags
        COMMAND ctags-exuberant -R    --fortran-kinds=+i --exclude="build/*" .
        COMMENT "Generating only ctags file for MFiX sources exclusively"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
elseif(CTAGS_CMD)
    add_custom_target ( tags
        COMMAND ctags -R    --fortran-kinds=+i --exclude="build/amrex/installdir/*" .
        COMMAND ctags -R -e --fortran-kinds=+i --exclude="build/amrex/installdir/*" .
        COMMENT "Generating tags files"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    # Some file systems are case-insensitive so the above will fail:
    #  the first command would overwrite the second!
    add_custom_target ( ctags
        COMMAND ctags -R    --fortran-kinds=+i --exclude="build/amrex/installdir/*" .
        COMMENT "Generating only ctags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    add_custom_target( etags
        COMMAND ctags -R -e --fortran-kinds=+i --exclude="build/amrex/installdir/*" .
        COMMENT "Generating only etags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    # Avoid "collisions" with AMReX => focus entirely on MFiX sources:
    add_custom_target ( mfix_ctags
        COMMAND ctags -R    --fortran-kinds=+i --exclude="build/*" .
        COMMENT "Generating only ctags file for MFiX sources exclusively"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
endif()

# Tools
set (MFIX_TOOLS_EXE ${PROJECT_BINARY_DIR}/Tools/) 
add_subdirectory (Tools)

# Tests
enable_testing ()
add_subdirectory ( tests )


# Add debug runs
add_subdirectory ( debug )

# Add a target to generate API documentation with Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND)
   set(DOC_DIR ${CMAKE_SOURCE_DIR}/doc)
   configure_file(${DOC_DIR}/Doxyfile.in ${DOC_DIR}/Doxyfile)
   add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOC_DIR}/Doxyfile
      WORKING_DIRECTORY ${DOC_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
endif (DOXYGEN_FOUND)

# Copy compile_commands.json from the build directory to the project root directory
if ( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
endif ()
