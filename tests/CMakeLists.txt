set( SGS_TESTS
     FLD01-x FLD01-y FLD01-z
     FLD02-x FLD02-y FLD02-z
     FLD03-x FLD03-y FLD03-z
     FLD04-x FLD04-y FLD04-z
     DEM01-x DEM01-y DEM01-z
     DEM02-x DEM02-y DEM02-z
     DEM03-x DEM03-y DEM03-z
     DEM04-x DEM04-y DEM04-z
     DEM05-x DEM05-y DEM05-z
     DEM06-x DEM06-y DEM06-z
     DEM06p-x DEM06p-y DEM06p-z
     DEM07
   )

set (BUILD_TEST build_tests)
add_custom_target(${BUILD_TEST})

foreach( THIS_TEST ${SGS_TESTS})

   STRING(REPLACE "/" "_" TEST_NAME ${THIS_TEST} )

   set(TEST_DIR ${CMAKE_BINARY_DIR}/tests/${THIS_TEST})
   set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/tests/${THIS_TEST})

   # Copy file to binary dir
   FILE( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )

   # Setup for building a test
   FILE(GLOB USR_OVERRIDES ${SOURCE_DIR}/*.f90)

   set(TEST_EXE mfix_${TEST_NAME})

   add_executable(${TEST_EXE} EXCLUDE_FROM_ALL ${MAIN} ${USR_OVERRIDES})

   target_include_directories(${TEST_EXE} PUBLIC ${TEST_DIR} ${CMAKE_Fortran_MODULE_DIRECTORY})

   set_target_properties( ${TEST_EXE}  PROPERTIES  Fortran_MODULE_DIRECTORY ${TEST_DIR}
      RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR} )

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXT_LINK_STRING} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Create test dir in build path if needed
   # Update test files if needed
   add_custom_command(TARGET ${TEST_EXE} PRE_BUILD
      COMMAND [ -d ${TEST_DIR} ] || mkdir ${TEST_DIR}
      COMMAND cp -R --update ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST}/* .
      WORKING_DIRECTORY  ${TEST_DIR}
      COMMENT "Updating test folder")

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_SGS)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=${ENABLE_MPI}
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_single
         WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
   # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

   # Setup for running the whole suite of tests via ctest
   add_test(NAME ${TEST_NAME}_SGS
      COMMAND make  ${TEST_RUN}
      WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})

endforeach( THIS_TEST ${SGS_TESTS} )




set( MGS_TESTS
     FLD01-x FLD01-y FLD01-z
     FLD02-x FLD02-y FLD02-z
     FLD03-x FLD03-y FLD03-z
     FLD04-x FLD04-y FLD04-z
     DEM06-x DEM06-y DEM06-z
     DEM07
   )

foreach( THIS_TEST ${MGS_TESTS})

   STRING(REPLACE "/" "_" TEST_NAME ${THIS_TEST} )

   set(TEST_DIR ${CMAKE_BINARY_DIR}/tests/${THIS_TEST})
   set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/tests/${THIS_TEST})

   # Copy file to binary dir
   FILE( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )

   # Setup for building a test
   FILE(GLOB USR_OVERRIDES ${SOURCE_DIR}/*.f90)

   set(TEST_EXE mfix_${TEST_NAME})


   target_include_directories(${TEST_EXE} PUBLIC ${TEST_DIR} ${CMAKE_Fortran_MODULE_DIRECTORY})

   set_target_properties( ${TEST_EXE}  PROPERTIES  Fortran_MODULE_DIRECTORY ${TEST_DIR}
      RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR} )

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXT_LINK_STRING} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Create test dir in build path if needed
   # Update test files if needed
   add_custom_command(TARGET ${TEST_EXE} PRE_BUILD
      COMMAND [ -d ${TEST_DIR} ] || mkdir ${TEST_DIR}
      COMMAND cp -R --update ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST}/* .
      WORKING_DIRECTORY  ${TEST_DIR}
      COMMENT "Updating test folder")

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_MGS)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=${ENABLE_MPI}
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_multiple
         WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
   # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

   # Setup for running the whole suite of tests via ctest
   add_test(NAME ${TEST_NAME}_MGS
      COMMAND make  ${TEST_RUN}
      WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})

endforeach( THIS_TEST ${MGS_TESTS} )



set( TGS_TESTS
     FLD01-x FLD01-y FLD01-z
     FLD02-x FLD02-y FLD02-z
     FLD03-x FLD03-y FLD03-z
     FLD04-x FLD04-y FLD04-z
     DEM06-x DEM06-y DEM06-z
     DEM07
   )


foreach( THIS_TEST ${TGS_TESTS})

   STRING(REPLACE "/" "_" TEST_NAME ${THIS_TEST} )

   set(TEST_DIR ${CMAKE_BINARY_DIR}/tests/${THIS_TEST})
   set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/tests/${THIS_TEST})

   # Copy file to binary dir
   FILE( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )

   # Setup for building a test
   FILE(GLOB USR_OVERRIDES ${SOURCE_DIR}/*.f90)

   set(TEST_EXE mfix_${TEST_NAME})


   target_include_directories(${TEST_EXE} PUBLIC ${TEST_DIR} ${CMAKE_Fortran_MODULE_DIRECTORY})

   set_target_properties( ${TEST_EXE}  PROPERTIES  Fortran_MODULE_DIRECTORY ${TEST_DIR}
      RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR} )

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXT_LINK_STRING} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Create test dir in build path if needed
   # Update test files if needed
   add_custom_command(TARGET ${TEST_EXE} PRE_BUILD
      COMMAND [ -d ${TEST_DIR} ] || mkdir ${TEST_DIR}
      COMMAND cp -R --update ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST}/* .
      WORKING_DIRECTORY  ${TEST_DIR}
      COMMENT "Updating test folder")

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_TGS)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=${ENABLE_MPI}
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_tiled
         WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
   # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

   # Setup for running the whole suite of tests via ctest
   add_test(NAME ${TEST_NAME}_TGS
      COMMAND make  ${TEST_RUN}
      WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})

endforeach( THIS_TEST ${TGS_TESTS} )


if ( ENABLE_MPI )
set( MGP_TESTS
     FLD01-x FLD01-y FLD01-z
     FLD02-x FLD02-y FLD02-z
     FLD03-x FLD03-y FLD03-z
     FLD04-x FLD04-y FLD04-z
   )


foreach( THIS_TEST ${MGP_TESTS})

   STRING(REPLACE "/" "_" TEST_NAME ${THIS_TEST} )

   set(TEST_DIR ${CMAKE_BINARY_DIR}/tests/${THIS_TEST})
   set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/tests/${THIS_TEST})

   # Copy file to binary dir
   FILE( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )

   # Setup for building a test
   FILE(GLOB USR_OVERRIDES ${SOURCE_DIR}/*.f90)

   set(TEST_EXE mfix_${TEST_NAME})

   target_include_directories(${TEST_EXE} PUBLIC ${TEST_DIR} ${CMAKE_Fortran_MODULE_DIRECTORY})

   set_target_properties( ${TEST_EXE}  PROPERTIES  Fortran_MODULE_DIRECTORY ${TEST_DIR}
      RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR} )

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXT_LINK_STRING} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Create test dir in build path if needed
   # Update test files if needed
   add_custom_command(TARGET ${TEST_EXE} PRE_BUILD
      COMMAND [ -d ${TEST_DIR} ] || mkdir ${TEST_DIR}
      COMMAND cp -R --update ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST}/* .
      WORKING_DIRECTORY  ${TEST_DIR}
      COMMENT "Updating test folder")

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_MGP)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=${ENABLE_MPI}
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_multiple
         WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
   # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

   # Setup for running the whole suite of tests via ctest
   add_test(NAME ${TEST_NAME}_MGP
      COMMAND make  ${TEST_RUN}
      WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})

endforeach( THIS_TEST ${MGP_TESTS} )
endif (ENABLE_MPI)
