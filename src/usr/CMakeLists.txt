#
# We need to treat separately the files in this folder to avoid duplication
# of symbols when building the tests
# If we do not do this and include these sources in the mfixcore library,
# CUDA compilation of the tests will fail
#
set(MFIX_UDFS_STUBS
   ${CMAKE_CURRENT_LIST_DIR}/mfix_usr.cpp
   ${CMAKE_CURRENT_LIST_DIR}/mfix_usr_drag.cpp
   ${CMAKE_CURRENT_LIST_DIR}/mfix_pc_usr.cpp
   CACHE INTERNAL "MFIX default UDFS")
