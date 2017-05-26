set( FCOMPARE_PAR_EXE fcompare_par.exe )
set( FEXTRACT_EXE     fextract.exe )

add_executable( ${FCOMPARE_PAR_EXE} fcompare_par.f90 )

add_executable( ${FEXTRACT_EXE} fextract.f90 )
target_link_libraries( ${FEXTRACT_EXE} ${MFIX_EXT_LINK_STRING} )
SET_TARGET_PROPERTIES(${FEXTRACT_EXE} PROPERTIES LINKER_LANGUAGE CXX)

# Communicate to main Cmakelists file
set( FCOMPARE_PAR_EXE ${FCOMPARE_PAR_EXE} PARENT_SCOPE )
set( FEXTRACT_EXE ${FEXTRACT_EXE} PARENT_SCOPE )
