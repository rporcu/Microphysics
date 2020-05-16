target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_usr.cpp
   MFIX_PC_usr.cpp
   usr_drag.cpp
   )
