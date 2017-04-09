set( TESTS  FLD01-x FLD01-y FLD01-z
            FLD02-x FLD02-y FLD02-z
            FLD03-x FLD03-y FLD03-z
            DEM01-x DEM01-y DEM01-z
            DEM02-x DEM02-y DEM02-z
            DEM03-x DEM03-y DEM03-z
            DEM04-x DEM04-y DEM04-z
            DEM05-x DEM05-y DEM05-z
            DEM06-x DEM06-y DEM06-z
   )

set (BUILD_TEST build_tests)
add_custom_target(${BUILD_TEST})

foreach( THIS_TEST ${TESTS})

   STRING(REPLACE "/" "_" TEST_NAME ${THIS_TEST} )

   set(TEST_DIR ${CMAKE_BINARY_DIR}/tests/${THIS_TEST})

   # Copy file to binary dir
   FILE( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )

   # Setup for building a test
   FILE(GLOB USR_OVERRIDES ${TEST_DIR}/*.f90)

   set(TEST_EXE mfix_${TEST_NAME})

   add_executable(${TEST_EXE} EXCLUDE_FROM_ALL ${MAIN} ${USR_OVERRIDES})

   target_include_directories(${TEST_EXE} PUBLIC ${TEST_DIR} ${CMAKE_Fortran_MODULE_DIRECTORY})

   set_target_properties( ${TEST_EXE}  PROPERTIES  Fortran_MODULE_DIRECTORY ${TEST_DIR}
      RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR} )

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXT_LINK_STRING} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME})
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${CMAKE_BINARY_DIR}/ThirdParty/bin
      WORKING_DIRECTORY  ${TEST_DIR} )

   # Setup for running the whole suite of tests via ctest
   add_test(NAME ${TEST_NAME}
      COMMAND make  ${TEST_RUN}
      WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})

endforeach( THIS_TEST ${TESTS} )




if ( ENABLE_MPI )
   set( MPI_TESTS FLD01-x FLD01-y FLD01-z
      FLD02-x FLD02-y FLD02-z
      FLD03-x FLD03-y FLD03-z
      )

   foreach( THIS_TEST ${MPI_TESTS})

      STRING(REPLACE "/" "_" TEST_NAME ${THIS_TEST} )

      set(TEST_DIR ${CMAKE_BINARY_DIR}/tests/${THIS_TEST})

      # Copy file to binary dir
      FILE( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )

      # Setup for building a test
      FILE(GLOB USR_OVERRIDES ${TEST_DIR}/*.f90)

      set(TEST_EXE mfix_MPI_${TEST_NAME})

      add_executable(${TEST_EXE} EXCLUDE_FROM_ALL ${MAIN} ${USR_OVERRIDES})

      target_include_directories(${TEST_EXE} PUBLIC ${TEST_DIR} ${CMAKE_Fortran_MODULE_DIRECTORY})

      set_target_properties( ${TEST_EXE}  PROPERTIES  Fortran_MODULE_DIRECTORY ${TEST_DIR}
	 RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR} )

      target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXT_LINK_STRING} )

      # Accumulate the build for each test into target "build_tests"
      add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

      # Setup for running the a single test
      set(TEST_RUN run_MPI_${TEST_NAME})
      add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
	 COMMAND ${TEST_DIR}/runtests-dmp.sh ${TEST_DIR}/${TEST_EXE} ${CMAKE_BINARY_DIR}/ThirdParty/bin
	 WORKING_DIRECTORY  ${TEST_DIR} )

      # Setup for running the whole suite of tests via ctest
      add_test(NAME MPI_${TEST_NAME}
	 COMMAND make  ${TEST_RUN}
	 WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
   endforeach( THIS_TEST ${MPI_TESTS} )

endif (ENABLE_MPI)

