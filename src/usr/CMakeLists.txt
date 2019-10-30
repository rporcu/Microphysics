target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   usr0.f90
   usr0_des.f90
   usr1.f90
   usr1_des.f90
   usr2.f90
   usr2_des.f90
   usr3.f90
   usr3_des.f90
   usr_mod.f90
   write_usr0.f90
   write_usr1.f90
   )
