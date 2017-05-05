set( F90_SRC fcompare_par.f90 )

set( FCOMPARE_PAR_EXE fcompare_par.exe ) 
add_executable( ${FCOMPARE_PAR_EXE} ${F90_SRC} )

# Communicate to main Cmakelists file
set( FCOMPARE_PAR_EXE ${FCOMPARE_PAR_EXE}  PARENT_SCOPE )
