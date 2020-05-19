################################################################################
# CSG/EB Comparison Tests
################################################################################

add_executable(functional_tests EXCLUDE_FROM_ALL
  benchmarks/05-cyl-fluidbed.cpp
  benchmarks/07-hopper.cpp
  benchmarks/08-cyclone.cpp
  benchmarks/02-settling.cpp
  benchmarks/tutorial-clr-prototype.cpp
  benchmarks/main.cpp
  )

target_include_directories(functional_tests PRIVATE
  ${CMAKE_SOURCE_DIR}/subprojects/PEGTL/include
  ${CMAKE_SOURCE_DIR}/subprojects/mfix-parser/src/csg/
  amrex
  )

# Special treatment for CUDA in functional_tests
if(ENABLE_CUDA)
   setup_target_for_cuda_compilation(functional_tests)
else()
   target_compile_features(functional_tests PRIVATE cxx_std_17)
endif()

get_filename_component(TUTORIAL_CLR_CSG
                       "${PROJECT_SOURCE_DIR}/tutorials/clr/prototype/clr.csg"
                       ABSOLUTE)
target_compile_definitions(functional_tests
                           PRIVATE TUTORIAL_CLR_CSG="${TUTORIAL_CLR_CSG}")

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endif()

target_link_libraries(functional_tests Catch2::Catch2 AMReX::amrex csg)
