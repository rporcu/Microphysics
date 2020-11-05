set( FJOIN_PAR_EXE    fjoin_par )
set( FEXTRACT_EXE     fextract )

add_executable( ${FJOIN_PAR_EXE} EXCLUDE_FROM_ALL C++/fjoin_par.cpp )
target_link_libraries( ${FJOIN_PAR_EXE} AMReX::amrex )

add_executable( ${FEXTRACT_EXE} EXCLUDE_FROM_ALL C++/fextract.cpp )
target_link_libraries( ${FEXTRACT_EXE} AMReX::amrex )

if (MFIX_CUDA)
   setup_target_for_cuda_compilation(${FEXTRACT_EXE})
endif ()

# Communicate to main Cmakelists file
set( MFIX_TOOLS ${FEXTRACT_EXE} ${FJOIN_PAR_EXE} PARENT_SCOPE  )
set( FEXTRACT_EXE ${FEXTRACT_EXE} PARENT_SCOPE )
set( FJOIN_PAR_EXE ${FJOIN_PAR_EXE} PARENT_SCOPE )
