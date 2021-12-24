target_include_directories( mfixcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_io.cpp
   mfix_build_info.cpp
   mfix_catalyst.cpp
   mfix_chk.cpp
   mfix_plt.cpp
   mfix_restart.cpp
   mfix_ascent.cpp
   mfix_pic2pic.cpp
   mfix_pic2dem.cpp
   mfix_rw.H
   mfix_rw.cpp
   )
