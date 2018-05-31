set( FCOMPARE_PAR_EXE fcompare_par )
set( FJOIN_PAR_EXE    fjoin_par )
set( FEXTRACT_EXE     fextract )
set( PLT_COMPARE_EXE  plt_compare_diff_grids )

add_executable( ${FCOMPARE_PAR_EXE} EXCLUDE_FROM_ALL F90/fcompare_par.f90 )

add_executable( ${FJOIN_PAR_EXE} EXCLUDE_FROM_ALL F90/fjoin_par.f90 )

add_executable( ${FEXTRACT_EXE} EXCLUDE_FROM_ALL C++/fextract.cpp )
target_link_libraries( ${FEXTRACT_EXE} ${MFIX_EXTRA_LINK_LINE} )

# fextract relies on some AMReX Fortran internals which can be linked against
# only if Fortran compiler libraries (like libgfortran.a) are included
# during the link phase. The following line takes care of that (not sure
# it's the cleanest way of doing it )
set_target_properties(${FEXTRACT_EXE} PROPERTIES LINKER_LANGUAGE Fortran )

# Communicate to main Cmakelists file
set( MFIX_F90_TOOLS ${FEXTRACT_EXE} ${FEXTRACT_EXE} ${FJOIN_PAR_EXE} PARENT_SCOPE  )
set( PLT_COMPARE_EXE ${PLT_COMPARE_EXE} PARENT_SCOPE )
set( FEXTRACT_EXE ${FEXTRACT_EXE} PARENT_SCOPE )
set( FCOMPARE_PAR_EXE ${FCOMPARE_PAR_EXE} PARENT_SCOPE )
set( FJOIN_PAR_EXE ${FJOIN_PAR_EXE} PARENT_SCOPE )
