set( FCOMPARE_PAR_EXE fcompare_par )
set( FEXTRACT_EXE     fextract )
set( PLT_COMPARE_EXE  plt_compare_diff_grids )

add_executable( ${FCOMPARE_PAR_EXE} EXCLUDE_FROM_ALL fcompare_par.f90 )

add_executable( ${FEXTRACT_EXE} EXCLUDE_FROM_ALL fextract.f90 )
target_link_libraries( ${FEXTRACT_EXE} ${MFIX_EXT_LINK_STRING} )
set_target_properties(${FEXTRACT_EXE} PROPERTIES LINKER_LANGUAGE CXX)


add_executable( ${PLT_COMPARE_EXE} EXCLUDE_FROM_ALL plt_compare_diff_grids.f90 )
target_link_libraries( ${PLT_COMPARE_EXE} ${MFIX_EXT_LINK_STRING} )
set_target_properties(${PLT_COMPARE_EXE} PROPERTIES LINKER_LANGUAGE CXX)


# Communicate to main Cmakelists file
set( MFIX_F90_TOOLS ${FEXTRACT_EXE} ${PLT_COMPARE_EXE} ${FEXTRACT_EXE} PARENT_SCOPE  )
set( PLT_COMPARE_EXE ${PLT_COMPARE_EXE} PARENT_SCOPE )
set( FEXTRACT_EXE ${FEXTRACT_EXE} PARENT_SCOPE )
set( FCOMPARE_PAR_EXE ${FCOMPARE_PAR_EXE} PARENT_SCOPE )
