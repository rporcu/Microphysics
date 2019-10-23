#
# Define a general build target to build all tests
#
set (BUILD_TEST build_tests)
add_custom_target(${BUILD_TEST})

#
# Define a macro to setup a test
#
macro (setup_test _test _suffix _env _input_file )

    string(REPLACE "/" "_" _test_name ${_test} )

    if (NOT TEST ${_test_name}_${_suffix})

        set(_test_dir  ${CMAKE_BINARY_DIR}/tests/${_test})
        set(_src_dir   ${PROJECT_SOURCE_DIR}/tests/${_test})

        # Copy file to binary dir
        file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${_test} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )

        # Setup for building a test
        FILE(GLOB _udf  ${_src_dir}/*.f90 ${_src_dir}/*.cpp)

        set(_test_exe mfix_${_test})

        # Create target if needed
        if (NOT TARGET ${_test_exe})
            add_executable(${_test_exe} EXCLUDE_FROM_ALL ${MAIN} ${_udf})

            target_include_directories(${_test_exe} PUBLIC ${_test_dir} ${CMAKE_Fortran_MODULE_DIRECTORY})

            set_target_properties( ${_test_exe}
                PROPERTIES
                Fortran_MODULE_DIRECTORY ${_test_dir}
                RUNTIME_OUTPUT_DIRECTORY ${_test_dir} )

            target_link_libraries(${_test_exe} ${MFIX_LIBNAME} )

            # Accumulate the build for each test into target "build_tests"
            add_dependencies ( ${BUILD_TEST} ${_test_exe} )

            # Create test dir in build path if needed
            # Update test files if needed

            add_custom_command(TARGET ${_test_exe} PRE_BUILD
                COMMAND mkdir -p ${_test_dir}
                COMMAND cp -R ${CMAKE_CURRENT_SOURCE_DIR}/${_test}/* .
                WORKING_DIRECTORY ${_test_dir}
                COMMENT "Copying files to test run directory")
        endif ()

        # Setup for running the a single test
        set(_test_run run_${_test_name}_${_suffix})

        add_custom_target(${_test_run} DEPENDS ${_test_exe}
            COMMAND ${CMAKE_COMMAND} -E env ${_env}
            ${_test_dir}/runtests.sh ${_test_dir}/${_test_exe} ${MFIX_TOOLS_EXE} ${_input_file}
            WORKING_DIRECTORY  ${_test_dir}
            )

        # Make sure that fcompare and fextract are always built before a test is run
        add_dependencies ( ${_test_run} ${MFIX_F90_TOOLS} )

        add_test(NAME ${_test_name}_${_suffix}
            COMMAND make  ${_test_run}
            WORKING_DIRECTORY ${PROJECT_BINDARY_DIR})

    endif ()

endmacro()



############################################################
#                                                          #
#                      SGS TESTS                           #
#                                                          #
############################################################
set( SGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   DEM01-x DEM01-y DEM01-z
   DEM03-x DEM03-y DEM03-z
   DEM04-x DEM04-y DEM04-z
   DEM05-x DEM05-y DEM05-z
   DEM06-x DEM06-y DEM06-z
   DEM07
   )


foreach( THIS_TEST ${SGS_TESTS})

    setup_test(
        ${THIS_TEST}
        "SGS"
        "ENABLE_MPI=0;ENABLE_OMP=0;OMP_NUM_THREADS=1"
        "inputs_single")

endforeach( THIS_TEST ${SGS_TESTS} )


############################################################
#                                                          #
#                      MGS & MGP TESTS                     #
#                                                          #
############################################################
set( MGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   DEM06-x DEM06-y DEM06-z
   DEM07
   )

foreach( THIS_TEST ${MGS_TESTS})

    setup_test(
        ${THIS_TEST}
        "MGS"
        "ENABLE_MPI=0;ENABLE_OMP=0;OMP_NUM_THREADS=1"
        "inputs_multiple")

    if (ENABLE_MPI)
        setup_test(
            ${THIS_TEST}
            "MGP"
            "ENABLE_MPI=1;ENABLE_OMP=0;OMP_NUM_THREADS=1"
            "inputs_multiple")
    endif ()

endforeach( THIS_TEST ${MGS_TESTS} )


############################################################
#                                                          #
#                      TGS & TGSO TESTS                    #
#                                                          #
############################################################
set( TGS_TESTS
   FLD01-x FLD01-y FLD01-z
   FLD02-x FLD02-y FLD02-z
   FLD03-x FLD03-y FLD03-z
   FLD04-x FLD04-y FLD04-z
   DEM06-x DEM06-y DEM06-z
   DEM07
   )


foreach( THIS_TEST ${TGS_TESTS})

    setup_test(
        ${THIS_TEST}
        "TGS"
        "ENABLE_MPI=0;ENABLE_OMP=0;OMP_NUM_THREADS=1"
        "inputs_tiled")

    if (ENABLE_OMP)
        setup_test(
            ${THIS_TEST}
            "TGSO"
            "ENABLE_MPI=0;ENABLE_OMP=1;OMP_NUM_THREADS=4"
            "inputs_tiled")
    endif ()

endforeach( THIS_TEST ${TGS_TESTS} )


############################################################
#                                                          #
#                      TGM TESTS                           #
#                                                          #
############################################################
if ( ENABLE_OMP AND ENABLE_MPI )

    set( TGM_TESTS
        FLD01-x FLD01-y FLD01-z
        FLD02-x FLD02-y FLD02-z
        FLD03-x FLD03-y FLD03-z
        FLD04-x FLD04-y FLD04-z
        DEM06-x DEM06-y DEM06-z
        )

    foreach( THIS_TEST ${TGM_TESTS})

        setup_test(
            ${THIS_TEST}
            "TGPO"
            "ENABLE_MPI=1;ENABLE_OMP=1;OMP_NUM_THREADS=2"
            "inputs_multiple_tiled")

    endforeach( )

endif()
