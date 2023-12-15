cmake_minimum_required(VERSION 3.20)

project(Microphysics
    VERSION 1.0.0
    DESCRIPTION "building primordial_chemistry network in Microphysics with CMake"
    LANGUAGES CXX C)

#----------------------------------------------------------------------------------------------------------------------
# general settings and options
#----------------------------------------------------------------------------------------------------------------------

find_package (Python COMPONENTS Interpreter Development)
find_package(Threads)

include(CTest)
include(CMakePackageConfigHelpers)

#use this function to build Microphysics using different EOS for different test problems in external codes
function(setup_target_for_microphysics_compilation network_name output_dir)

  set(VODEparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/VODE/_parameters")
  set(integrationparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/_parameters")

endfunction()


#set a cache variable that controls whether
#we want to fetch and link amrex or not
set(BUILD_AMReX false CACHE BOOL "Do you want to build and link amrex? (true/false)")
message("Building and linking amrex -- ${BUILD_AMReX}")

if(BUILD_AMReX)
   #fetching amrex
   set (AMReX_TAG "development")
   include(FetchContent)
   FetchContent_Declare(
     amrex
     GIT_REPOSITORY https://github.com/AMReX-Codes/amrex
     GIT_TAG "${AMReX_TAG}"
   )

   # CMake 3.14+
   FetchContent_MakeAvailable(amrex)

   include_directories(${amrex_INCLUDE_DIRS_RET})
   link_libraries(AMReX::amrex)
   message("Built and fetched amrex in Microphysics")
else()
   message("Not building and fetching amrex in Microphysics")
endif()

#set a cache variable that controls whether
#we want to build the unit test or not
add_compile_options(-Werror -Wall -Wextra)

#setting sourcefiles and directories needed to make the test here
#so that they are accessible to codes using
#Microphysics as a submodule
#include parent_scope for sourcefiles and directories
#when built as a submodule in an external code
#otherwise these are not visible to the external code
#this will generate a warning but it will build successfully
#do not need unit_test paramfiles when unit_test is not built
