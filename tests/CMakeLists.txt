set( SGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z  
   FLD04-x FLD04-y FLD04-z
   DEM01-x DEM01-y DEM01-z DEM01-eb
   DEM02-x DEM02-y DEM02-z DEM02-eb
   DEM03-x DEM03-y DEM03-z DEM03-eb
   DEM04-x DEM04-y DEM04-z DEM04-eb
   DEM05-x DEM05-y DEM05-z DEM05-eb
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

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Create test dir in build path if needed
   # Update test files if needed
   if(APPLE)
      add_custom_command(TARGET ${TEST_EXE} PRE_BUILD
         COMMAND [ -d ${TEST_DIR} ] || mkdir ${TEST_DIR}
         COMMAND cp -R ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST}/* .
         WORKING_DIRECTORY ${TEST_DIR}
         COMMENT "Updation test folder")
   else()
      add_custom_command(TARGET ${TEST_EXE} PRE_BUILD
         COMMAND [ -d ${TEST_DIR} ] || mkdir ${TEST_DIR}
         COMMAND cp -R --update ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST}/* .
         WORKING_DIRECTORY  ${TEST_DIR}
         COMMENT "Updating test folder")
   endif()

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_SGS)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=0 ENABLE_OMP=0 OMP_NUM_THREADS=1
      ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_single
      WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
   # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

   # Setup for running the whole suite of tests via ctest
   if (ENABLE_STESTS)
      add_test(NAME ${TEST_NAME}_SGS
	 COMMAND make  ${TEST_RUN}
	 WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
   endif ()
endforeach( THIS_TEST ${SGS_TESTS} )


set( MGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   DEM06-x DEM06-y DEM06-z
   DEM06p-x DEM06p-y DEM06p-z
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

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )


   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_MGS)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=0 ENABLE_OMP=0 OMP_NUM_THREADS=1
      ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_multiple
      WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
   # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

   # Setup for running the whole suite of tests via ctest
   if (ENABLE_STESTS)
      add_test(NAME ${TEST_NAME}_MGS
	 COMMAND make  ${TEST_RUN}
	 WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
   endif ()
   
endforeach( THIS_TEST ${MGS_TESTS} )



set( TGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   DEM06-x DEM06-y DEM06-z
   DEM06p-x DEM06p-y DEM06p-z
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

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_TGS)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=0 ENABLE_OMP=0 OMP_NUM_THREADS=1
      ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_tiled
      WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
   # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

   # Setup for running the whole suite of tests via ctest
   if (ENABLE_STESTS)
      add_test(NAME ${TEST_NAME}_TGS
	 COMMAND make  ${TEST_RUN}
	 WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
   endif ()
endforeach( THIS_TEST ${TGS_TESTS} )


if ( ENABLE_OMP )

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

      target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

      # Accumulate the build for each test into target "build_tests"
      add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

      # Setup for running the a single test
      set(TEST_RUN run_${TEST_NAME}_TGSO)
      add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
	 COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=0 ENABLE_OMP=1 OMP_NUM_THREADS=4
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_tiled
         WORKING_DIRECTORY  ${TEST_DIR} )

      # Make sure that fcompare and fextract are always built before a test is run
      add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
      # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

      # Setup for running the whole suite of tests via ctest
      if (ENABLE_STESTS)
	 add_test(NAME ${TEST_NAME}_TGSO
	    COMMAND make  ${TEST_RUN}
	    WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
      endif ()
      
   endforeach( THIS_TEST ${TGS_TESTS} )

endif ()


set( TGM_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   DEM06-x DEM06-y DEM06-z
   DEM06p-x DEM06p-y DEM06p-z
   )

if (ENABLE_OMP AND ENABLE_MPI)
   foreach( THIS_TEST ${TGM_TESTS})

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

      target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

      # Accumulate the build for each test into target "build_tests"
      add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

      # Setup for running the a single test
      set(TEST_RUN run_${TEST_NAME}_TGPO)
      add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
	 COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=1 ENABLE_OMP=1 OMP_NUM_THREADS=2
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_multiple_tiled
         WORKING_DIRECTORY  ${TEST_DIR} )

      # Make sure that fcompare and fextract are always built before a test is run
      add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
      # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

      # Setup for running the whole suite of tests via ctest
      if (ENABLE_STESTS) 
	 add_test(NAME ${TEST_NAME}_TGPO
	    COMMAND make  ${TEST_RUN}
	    WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
      endif() 
   endforeach( THIS_TEST ${TGM_TESTS} )
endif()



if ( ENABLE_MPI )
   set( MGP_TESTS
      FLD01-x FLD01-y FLD01-z
      FLD02-x FLD02-y FLD02-z
      FLD03-x FLD03-y FLD03-z
      FLD04-x FLD04-y FLD04-z
      DEM06-x DEM06-y DEM06-z
      DEM06p-x DEM06p-y DEM06p-z
      DEM07
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

      target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

      # Accumulate the build for each test into target "build_tests"
      add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

      # Setup for running the a single test
      set(TEST_RUN run_${TEST_NAME}_MGP)
      add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
	 COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=1 ENABLE_OMP=0 OMP_NUM_THREADS=1
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_multiple
         WORKING_DIRECTORY  ${TEST_DIR} )

      # Make sure that fcompare and fextract are always built before a test is run
      add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
      # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

      # Setup for running the whole suite of tests via ctest
      if (ENABLE_STESTS) 
	 add_test(NAME ${TEST_NAME}_MGP
	    COMMAND make  ${TEST_RUN}
	    WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
      endif ()
   endforeach( THIS_TEST ${MGP_TESTS} )
endif (ENABLE_MPI)

####################################################################################################
#                                                                                                  #
#                                                                                                  #
#                                    PROJECTION METHOD TESTS                                       #
#                                                                                                  #
#                                                                                                  #
####################################################################################################
set( SGS_PTESTS
   FLD01-x_P FLD01-y_P FLD01-z_P
   FLD02-x_P FLD02-y_P FLD02-z_P
   FLD03-x_P FLD03-y_P FLD03-z_P
   FLD04-x_P FLD04-y_P FLD04-z_P 
   )

# set (BUILD_TEST build_tests)
# add_custom_target(${BUILD_TEST})

foreach( THIS_TEST ${SGS_PTESTS})

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

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Create test dir in build path if needed
   # Update test files if needed
   if(APPLE)
      add_custom_command(TARGET ${TEST_EXE} PRE_BUILD
         COMMAND [ -d ${TEST_DIR} ] || mkdir ${TEST_DIR}
         COMMAND cp -R ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST}/* .
         WORKING_DIRECTORY ${TEST_DIR}
         COMMENT "Updation test folder")
   else()
      add_custom_command(TARGET ${TEST_EXE} PRE_BUILD
         COMMAND [ -d ${TEST_DIR} ] || mkdir ${TEST_DIR}
         COMMAND cp -R --update ${CMAKE_CURRENT_SOURCE_DIR}/${THIS_TEST}/* .
         WORKING_DIRECTORY  ${TEST_DIR}
         COMMENT "Updating test folder")
   endif()

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_SGS )
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=0 ENABLE_OMP=0 OMP_NUM_THREADS=1
      ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_single
      WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )

   if (ENABLE_PTESTS)
      # Setup for running the whole suite of tests via ctest
      add_test(NAME ${TEST_NAME}_SGS
	 COMMAND make  ${TEST_RUN}
	 WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
   endif ()

endforeach ()





set( MGS_PTESTS
   FLD01-x_P FLD01-y_P FLD01-z_P
   FLD02-x_P FLD02-y_P FLD02-z_P
   FLD03-x_P FLD03-y_P FLD03-z_P
   FLD04-x_P FLD04-y_P FLD04-z_P
   #  DEM06-x_P DEM06-y DEM06-z
   )

foreach( THIS_TEST ${MGS_PTESTS})

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

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_MGS)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=0 ENABLE_OMP=0 OMP_NUM_THREADS=1
      ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_multiple
      WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )

   if (ENABLE_PTESTS) 
      # Setup for running the whole suite of tests via ctest
      add_test(NAME ${TEST_NAME}_MGS
         COMMAND make  ${TEST_RUN}
         WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
   endif ()

endforeach()



set( TGS_PTESTS
   FLD01-x_P FLD01-y_P FLD01-z_P
   FLD02-x_P FLD02-y_P FLD02-z_P
   FLD03-x_P FLD03-y_P FLD03-z_P
   FLD04-x_P FLD04-y_P FLD04-z_P
   )


foreach( THIS_TEST ${TGS_PTESTS})

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

   target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

   # Accumulate the build for each test into target "build_tests"
   add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

   # Setup for running the a single test
   set(TEST_RUN run_${TEST_NAME}_TGS)
   add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
      COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=0 ENABLE_OMP=0 OMP_NUM_THREADS=1
      ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_tiled
      WORKING_DIRECTORY  ${TEST_DIR} )

   # Make sure that fcompare and fextract are always built before a test is run
   add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
   # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

   # Setup for running the whole suite of tests via ctest
   if (ENABLE_PTESTS)
      add_test(NAME ${TEST_NAME}_TGS
	 COMMAND make  ${TEST_RUN}
	 WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
   endif ()
endforeach( )


set( TGM_PTESTS
   FLD01-x_P FLD01-y_P FLD01-z_P
   FLD02-x_P FLD02-y_P FLD02-z_P
   FLD03-x_P FLD03-y_P FLD03-z_P
   FLD04-x_P FLD04-y_P FLD04-z_P
   )

if (ENABLE_OMP AND ENABLE_MPI)
   foreach( THIS_TEST ${TGM_PTESTS})

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

      target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

      # Accumulate the build for each test into target "build_tests"
      add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

      # Setup for running the a single test
      set(TEST_RUN run_${TEST_NAME}_TGPO)
      add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
	 COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=1 ENABLE_OMP=1 OMP_NUM_THREADS=2
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_multiple_tiled
         WORKING_DIRECTORY  ${TEST_DIR} )

      # Make sure that fcompare and fextract are always built before a test is run
      add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
      # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

      # Setup for running the whole suite of tests via ctest
      if (ENABLE_PTESTS) 
	 add_test(NAME ${TEST_NAME}_TGPO
	    COMMAND make  ${TEST_RUN}
	    WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
      endif() 
   endforeach()
endif()


set( MGP_PTESTS
   FLD01-x_P FLD01-y_P FLD01-z_P
   FLD02-x_P FLD02-y_P FLD02-z_P
   FLD03-x_P FLD03-y_P FLD03-z_P
   FLD04-x_P FLD04-y_P FLD04-z_P
   )


if ( ENABLE_MPI )

   foreach( THIS_TEST ${MGP_PTESTS})

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

      target_link_libraries(${TEST_EXE} mfixcore ${MFIX_EXTRA_LINK_LINE} )

      # Accumulate the build for each test into target "build_tests"
      add_dependencies ( ${BUILD_TEST} ${TEST_EXE} )

      # Setup for running the a single test
      set(TEST_RUN run_${TEST_NAME}_MGP)
      add_custom_target(${TEST_RUN} DEPENDS ${TEST_EXE}
	 COMMAND ${CMAKE_COMMAND} -E env ENABLE_MPI=1 ENABLE_OMP=0 OMP_NUM_THREADS=1
         ${TEST_DIR}/runtests.sh ${TEST_DIR}/${TEST_EXE} ${MFIX_TOOLS_EXE} inputs_multiple
         WORKING_DIRECTORY  ${TEST_DIR} )

      # Make sure that fcompare and fextract are always built before a test is run
      add_dependencies ( ${TEST_RUN} ${MFIX_F90_TOOLS} )
      # add_dependencies ( ${TEST_RUN} ${FCOMPARE_PAR_EXE} ${FEXTRACT_EXE})

      # Setup for running the whole suite of tests via ctest
      if (ENABLE_PTESTS) 
	 add_test(NAME ${TEST_NAME}_MGP
	    COMMAND make  ${TEST_RUN}
	    WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})
      endif ()
   endforeach ( )
   
endif (ENABLE_MPI)
