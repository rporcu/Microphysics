function ( build_amrex )

   message( STATUS "Configuring build to compile AMReX as part of MFIX")

   # First make sure the submodule has been pulled and updated
   set( GIT_OUTPUT  "" )

   execute_process( COMMAND git submodule init  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE  GIT_OUTPUT OUTPUT_STRIP_TRAILING_WHITESPACE )

   if ( NOT "${GIT_OUTPUT}" STREQUAL "" )
      message( STATUS ${GIT_OUTPUT} )
   endif ( NOT "${GIT_OUTPUT}" STREQUAL "" )

   execute_process( COMMAND git submodule update  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE  GIT_OUTPUT OUTPUT_STRIP_TRAILING_WHITESPACE )

   if ( NOT "${GIT_OUTPUT}" STREQUAL "" )
      message( STATUS "${GIT_OUTPUT}" )
   endif ( NOT "${GIT_OUTPUT}" STREQUAL "" )


   # Include cmake config files to build external projects
   include(ExternalProject)

   # Set AMReX paths to use in mfix config
   set(AMREX_SOURCE_DIR  ${CMAKE_SOURCE_DIR}/ThirdParty/amrex )
   set(AMREX_ROOT        ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty PARENT_SCOPE )
   set(AMREX_LIB_DIR     ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty/lib PARENT_SCOPE )
   set(AMREX_BIN_DIR     ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty/bin PARENT_SCOPE )
   set(AMREX_GIT         ${CMAKE_SOURCE_DIR}/ThirdParty/amrex PARENT_SCOPE )
   set(AMREX_TOOLS       ${CMAKE_SOURCE_DIR}/ThirdParty/amrex/Tools PARENT_SCOPE )
   set(AMREX_INCLUDES    ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty/include PARENT_SCOPE )
   
   # Set mpi-specific flags ( defined by CCSEFind )
   if (ENABLE_MPI)
      find_package(MPI REQUIRED)
      set( MPI_INCLUDES ${MPI_CXX_INCLUDE_PATH} PARENT_SCOPE )
      set( MPI_LFLAGS   ${MPI_CXX_LINK_FLAGS}   PARENT_SCOPE )
      set( MPI_LIBS     ${MPI_CXX_LIBRARIES}    PARENT_SCOPE )
   endif (ENABLE_MPI)

   # Set manually list of libraries
   # ( not using cmake find utilities in this case)
   set( AMREX_LIBRARIES fboxlib;cboxlib;fboxlib;cfboxlib;box_camrdata PARENT_SCOPE )
   
   # Add Cmake Tools from AMReX
   set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/ThirdParty/amrex/Tools/CMake/)

   ExternalProject_Add(
      amrex
      PREFIX          ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty
      SOURCE_DIR      ${AMREX_SOURCE_DIR}
      INSTALL_DIR     ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty
      CMAKE_ARGS  -DENABLE_MPI=${ENABLE_MPI} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBL_USE_PARTICLES=1 -DBL_SPACEDIM=${BL_SPACEDIM} -DENABLE_OpenMP=${ENABLE_OpenMP} -DBL_PRECISION=${BL_PRECISION}
        -DENABLE_PROFILING=${ENABLE_PROFILING} -DENABLE_BACKTRACE=${ENABLE_BACKTRACE}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DMPI_CXX_INCLUDE_PATH=${MPI_CXX_INCLUDE_PATH}
  -DMPI_CXX_LIBRARIES=${MPI_CXX_LIBRARIES}
      LOG_CONFIGURE 1
      LOG_BUILD 1
      LOG_INSTALL 1
      )

   # AMReX defines for compilation
   include(CCSEOptions)

endfunction (build_amrex)

function (find_amrex)

   message( STATUS "Configuring build to use external AMReX install")

   # Exit if AMREX_HOME is not set
   if ( "$ENV{AMREX_HOME}" STREQUAL "" )
      message( FATAL_ERROR "AMREX_HOME is not set" )
   endif ( "$ENV{AMREX_HOME}" STREQUAL "" )

   message( STATUS "AMREX_HOME is set to $ENV{AMREX_HOME}")

   # Exit if path given by AMREX_HOME is not set
   check_path( $ENV{AMREX_HOME} FATAL_ERROR )

   # Check if cmake folder is in the installation path
   check_path( $ENV{AMREX_HOME}/cmake FATAL_ERROR )

   # Check if FindCCSE.cmake is in path given by AMREX_TOOLS
   check_path( $ENV{AMREX_HOME}/cmake/FindCCSE.cmake FATAL_ERROR )

   # Issue warning
   if (ENABLE_MPI)
      message( STATUS "")
      message( STATUS "==================== WARNING ====================" )
      message( STATUS " MFIX is being compiled with MPI enabled" )
      message( STATUS " Make sure that AMReX is built with ENABLE_MPI=1" )
      message( STATUS " to avoid linking problems" )
      message( STATUS "=================================================" )
      message( STATUS "")
   endif (ENABLE_MPI)

   set(CCSE_DIR $ENV{AMREX_HOME})
   list(APPEND CMAKE_MODULE_PATH $ENV{AMREX_HOME}/cmake)
   find_package(CCSE REQUIRED)

   # Set AMReX paths to use in mfix config
   set(AMREX_ROOT        $ENV{AMREX_HOME}       PARENT_SCOPE )
   set(AMREX_BIN_DIR     $ENV{AMREX_HOME}/bin   PARENT_SCOPE )
   set(AMREX_TOOLS       $ENV{AMREX_HOME}/Tools PARENT_SCOPE )
   set(AMREX_INCLUDES    ${CCSE_INCLUDE_DIR}    PARENT_SCOPE )
   set(AMREX_LIBRARIES   ${CCSE_LIBRARIES}      PARENT_SCOPE )
   set(AMREX_LIB_DIR     ${CCSE_LIBRARY_DIR}    PARENT_SCOPE )

   # Set mpi-specific flags ( defined by CCSEFind )
   if (ENABLE_MPI)
      set( MPI_INCLUDES ${MPI_CXX_INCLUDE_PATH} PARENT_SCOPE )
      set( MPI_LFLAGS   ${MPI_CXX_LINK_FLAGS}   PARENT_SCOPE )
      set( MPI_LIBS     ${MPI_CXX_LIBRARIES}    PARENT_SCOPE )
   endif (ENABLE_MPI)

endfunction (find_amrex)
