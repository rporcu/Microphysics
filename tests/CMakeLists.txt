set( TESTS  FLD01-x FLD01-y FLD01-z FLD02-x FLD02-y FLD02-z 
            DEM01-x DEM01-y DEM01-z DEM02-x DEM02-y DEM02-z
            DEM03-x DEM03-y DEM03-z DEM04-x DEM04-y DEM04-z
            DEM05-x DEM05-y DEM05-z DEM06-y DEM06-z DEM06-x
   )

set (BUILD_TEST build_test)
add_custom_target(${BUILD_TEST})

foreach( rel_testdir ${TESTS})

   STRING(REPLACE "/" "_" TEST_NAME ${rel_testdir} )

   set(testdir ${CMAKE_CURRENT_SOURCE_DIR}/${rel_testdir})

   FILE(GLOB USR_OVERRIDES ${testdir}/*.f90)

   set(TEST_EXE mfix_${TEST_NAME})

   add_executable(${TEST_EXE} EXCLUDE_FROM_ALL ${MAIN} ${USR_OVERRIDES})

   target_include_directories(${TEST_EXE} PUBLIC ${testdir} ${CMAKE_Fortran_MODULE_DIRECTORY})
   set_target_properties(${TEST_EXE} PROPERTIES Fortran_MODULE_DIRECTORY ${testdir})

   add_test(NAME ${TEST_NAME}
      COMMAND ${CMAKE_COMMAND}
      -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
      -DCMAKE_COMMAND=${CMAKE_COMMAND}
      -DTEST_EXE=${TEST_EXE}
      -DTEST_NAME=${TEST_NAME}
      -Dtestdir=${testdir}
      -P ${PROJECT_SOURCE_DIR}/runtests.cmake)


   add_dependencies ( ${BUILD_TEST} ${TEST_EXE})

   set_tests_properties ( ${TEST_NAME}
      PROPERTIES ENVIRONMENT "MFIX=${CMAKE_BINARY_DIR}/${TEST_EXE}"
      )

   target_link_libraries(${TEST_EXE} mfixcore ${AMREX_LIBRARIES} ${MPI_LIBRARIES})

endforeach( rel_testdir ${TESTS})