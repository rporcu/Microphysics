target_include_directories(mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources( mfixcore
   PRIVATE
   desnamelist.inc
   initial_conditions.inc
   usr_hooks.inc
   )
