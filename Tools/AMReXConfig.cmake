function ( build_amrex )

   message( STATUS "Configuring build to compile AMReX as part of MFIX")

   # Include cmake config files to build external projects
   include(ExternalProject)

   # Set AMReX paths to use in mfix config
   set(AMREX_SOURCE_DIR  ${PROJECT_SOURCE_DIR}/ThirdParty/amrex)
   set(AMREX_ROOT        ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty )
   set(AMREX_LIB_DIR     ${AMREX_ROOT}/lib)
   set(AMREX_GIT         ${AMREX_SOURCE_DIR})
   set(AMREX_TOOLS       ${AMREX_SOURCE_DIR}/Tools)
   set(AMREX_INCLUDES    ${AMREX_ROOT}/include)
   set(AMREX_LIBRARIES    )

   # Set manually list of libraries 
   # ( not using cmake find utilities in this case)
   if (ENABLE_MPI)
      set( AMREX_LIBRARIES fboxlib;cboxlib;fboxlib;cfboxlib;box_camrdata)
   else ()
      set( AMREX_LIBRARIES cboxlib;fboxlib;cfboxlib;box_camrdata)
   endif ()

   # Add Cmake Tools from AMReX
   set(CMAKE_MODULE_PATH ${AMREX_SOURCE_DIR}/Tools/CMake/)

   ExternalProject_Add(
      amrex
      PREFIX          ${AMREX_ROOT}
      SOURCE_DIR      ${AMREX_SOURCE_DIR}
      INSTALL_DIR     ${AMREX_ROOT}
      CMAKE_ARGS      -DENABLE_MPI=${ENABLE_MPI} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DBL_USE_PARTICLES=1 -DBL_SPACEDIM=3 -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      LOG_DOWNLOAD 1
      LOG_UPDATE 1
      )

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
   set(AMREX_LIB_DIR     $ENV{AMREX_HOME}/lib   PARENT_SCOPE )
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
      