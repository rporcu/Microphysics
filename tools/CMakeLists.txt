set ( FJOIN_PAR_EXE    fjoin_par )
set ( FEXTRACT_EXE     fextract )

add_executable( ${FJOIN_PAR_EXE} EXCLUDE_FROM_ALL F90/fjoin_par.f90 )

add_executable( ${FEXTRACT_EXE} EXCLUDE_FROM_ALL C++/fextract.cpp )
target_link_libraries( ${FEXTRACT_EXE} AMReX::amrex )

if (ENABLE_CUDA)
   setup_target_for_cuda_compilation(${FEXTRACT_EXE})
endif ()

# Communicate to main Cmakelists file
set ( MFIX_F90_TOOLS ${FEXTRACT_EXE} ${FEXTRACT_EXE} ${FJOIN_PAR_EXE} PARENT_SCOPE  )
set ( FEXTRACT_EXE ${FEXTRACT_EXE} PARENT_SCOPE )
set ( FJOIN_PAR_EXE ${FJOIN_PAR_EXE} PARENT_SCOPE )
