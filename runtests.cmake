macro(EXEC_CHECK)
  execute_process(COMMAND ${ARGV}
    WORKING_DIRECTORY ${testdir}
    RESULT_VARIABLE
    CMD_RESULT)
  if(CMD_RESULT)
    message(FATAL_ERROR "Error running ${ARGV} in ${testdir} got result ${CMD_RESULT}")
  endif()
endmacro()
exec_check(${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${TEST_EXE})
exec_check(${testdir}/runtests.sh ${CMAKE_CURRENT_BINARY_DIR}/${TEST_EXE})
