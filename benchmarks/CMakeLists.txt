set( BENCHMARKS   )

set (BUILD_BENCHMARK build_benchmarks)
add_custom_target(${BUILD_BENCHMARK})

foreach( THIS_BENCHMARK ${BENCHMARKS})

   STRING(REPLACE "/" "_" BENCHMARK_NAME ${THIS_BENCHMARK} )

   set(BENCHMARK_DIR ${CMAKE_BINARY_DIR}/benchmarks/${THIS_BENCHMARK})
   set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/benchmarks/${THIS_BENCHMARK})

   # Copy file to binary dir
   FILE( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_BENCHMARK} DESTINATION ${CMAKE_BINARY_DIR}/benchmarks/ )

   # Setup for building a test
   FILE(GLOB USR_OVERRIDES ${SOURCE_DIR}/*.f90)

   set(BENCHMARK_EXE mfix_BM_${BENCHMARK_NAME})

   add_executable(${BENCHMARK_EXE} EXCLUDE_FROM_ALL ${MAIN} ${USR_OVERRIDES})

   target_include_directories(${BENCHMARK_EXE} PUBLIC ${BENCHMARK_DIR} ${CMAKE_Fortran_MODULE_DIRECTORY})

   set_target_properties( ${BENCHMARK_EXE}  PROPERTIES  Fortran_MODULE_DIRECTORY ${BENCHMARK_DIR}
      RUNTIME_OUTPUT_DIRECTORY ${BENCHMARK_DIR} )

   target_link_libraries(${BENCHMARK_EXE} mfixcore ${MFIX_EXT_LINK_STRING} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_BENCHMARK} ${BENCHMARK_EXE} )

   # Setup for running the a single test
   set(BENCHMARK_RUN run_BM_${BENCHMARK_NAME})
   add_custom_target(${BENCHMARK_RUN} DEPENDS ${BENCHMARK_EXE}
      COMMAND ${BENCHMARK_DIR}/runtests.sh ${BENCHMARK_DIR}/${BENCHMARK_EXE} ${CMAKE_BINARY_DIR}/ThirdParty/bin
      WORKING_DIRECTORY  ${BENCHMARK_DIR} )

   # Setup for running the whole suite of tests via ctest
   add_test(NAME BM_${BENCHMARK_NAME}
      COMMAND make  ${BENCHMARK_RUN}
      WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})

endforeach( THIS_BENCHMARK ${BENCHMARKS} )


if ( ENABLE_MPI )
  set( MPI_BENCHMARKS FLD01-x FLD01-y FLD01-z)

  foreach( THIS_BENCHMARK ${MPI_BENCHMARKS})

    STRING(REPLACE "/" "_" BENCHMARK_NAME ${THIS_BENCHMARK} )

    set(BENCHMARK_DIR ${CMAKE_BINARY_DIR}/benchmarks/${THIS_BENCHMARK})

    # Copy file to binary dir
    FILE( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_BENCHMARK} DESTINATION ${CMAKE_BINARY_DIR}/benchmarks/ )

    # Setup for building a test
    FILE(GLOB USR_OVERRIDES ${BENCHMARK_DIR}/*.f90)

    set(BENCHMARK_EXE mfix_MPI_BM_${BENCHMARK_NAME})

    add_executable(${BENCHMARK_EXE} EXCLUDE_FROM_ALL ${MAIN} ${USR_OVERRIDES})

    target_include_directories(${BENCHMARK_EXE} PUBLIC ${BENCHMARK_DIR} ${CMAKE_Fortran_MODULE_DIRECTORY})

    set_target_properties( ${BENCHMARK_EXE}  PROPERTIES  Fortran_MODULE_DIRECTORY ${BENCHMARK_DIR}
      RUNTIME_OUTPUT_DIRECTORY ${BENCHMARK_DIR} )

    target_link_libraries(${BENCHMARK_EXE} mfixcore ${MFIX_EXT_LINK_STRING} )

    # Accumulate the build for each test into target "build_tests"
    add_dependencies ( ${BUILD_BENCHMARK} ${BENCHMARK_EXE} )

    # Setup for running the a single test
    set(BENCHMARK_RUN run_MPI_BM_${BENCHMARK_NAME})
    add_custom_target(${BENCHMARK_RUN} DEPENDS ${BENCHMARK_EXE}
      COMMAND ${BENCHMARK_DIR}/runtests-dmp.sh ${BENCHMARK_DIR}/${BENCHMARK_EXE} ${CMAKE_BINARY_DIR}/ThirdParty/bin
      WORKING_DIRECTORY  ${BENCHMARK_DIR} )

    # Setup for running the whole suite of tests via ctest
    add_test(NAME MPI_BM_${BENCHMARK_NAME}
      COMMAND make  ${BENCHMARK_RUN}
      WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
  endforeach( THIS_BENCHMARK ${MPI_BM_BENCHMARKS} )

endif (ENABLE_MPI)
