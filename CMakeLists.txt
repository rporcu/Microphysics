cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
enable_language(Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if (CMAKE_Fortran_COMPILER MATCHES ".*gfortran.*")
  set (CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -g -Wall -O2" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS_DEBUG   "-fmax-errors=2 -std=c++11 -g -Wall -Wextra -O0" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-g -std=f2008 -fcheck=all -Wall -O2" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fmax-errors=2 -g -std=f2008 -fcheck=all -Wall -Wextra -Wimplicit-interface" CACHE STRING "" FORCE)

endif (CMAKE_Fortran_COMPILER MATCHES ".*gfortran.*")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set some defaults
set(BL_SPACEDIM 3 CACHE INT "Dimension of AMReX build")
set(ENABLE_MPI 0 CACHE INT "Enable build with MPI")
set(ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
set(BL_PRECISION "DOUBLE" CACHE INT "Precision of AMReX build")
set(BL_USE_PARTICLES 1 CACHE INT "Include Particles classes in AMReX build")
set(ENABLE_PROFILING 0 CACHE INT "Include profiling information in AMReX build")
set(ENABLE_BACKTRACE 1 CACHE INT "Include backtrace information in AMReX build")


# This is where the CCSE lib is installed
set(CCSE_DIR $ENV{AMREX_HOME})

set(CMAKE_MODULE_PATH ${CCSE_DIR}/cmake)
find_package(CCSE REQUIRED)

# Find CCSE check for mpi libraries if ENABLE_MPI is on.
# Thus all the variables defined by FindMPI are available
# Notice that only need c++-related variables are needed, since
# mfix itself does not call mpi function explicitly 
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(src)
include_directories(src/include)
link_directories(${CCSE_LIBRARY_DIR})


set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)

project(MFIXEXA)

# Generate buildInfo.cpp
if ( $ENV{AMREX_HOME} )
  message( FATAL_ERROR "AMREX_HOME is not set" )
endif()

add_custom_command(
  COMMAND  $(AMREX_HOME)/Tools/C_scripts/makebuildinfo_C.py
  --amrex_home "$ENV{AMREX_HOME}"
  --COMP ${CMAKE_C_COMPILER_ID} --COMP_VERSION ${CMAKE_C_COMPILER_VERSION}
  --FCOMP ${CMAKE_Fortran_COMPILER_ID} --FCOMP_VERSION ${CMAKE_C_COMPILER_VERSION}
  --GIT ". $ENV{AMREX_HOME}"
  OUTPUT AMReX_buildInfo.cpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating buildInfo.cpp" )


set(CXX_sources
  src/main.cpp
  src/mfix_level.cpp
  src/solve_bicgstab.cpp
  src/ParticleContainer.cpp
  src/io/mfixIO.cpp )
set(CXX_includes
  src/mfix_level.H
  src/ParticleContainer.H
  src/include/mfix_F.H)
set(F90_sources
  src/adjust_a_g.f90
  src/adjust_dt_mod.f90
  src/adjust_leq.f90
  src/bc_mod.f90
  src/calc_cell.f90
  src/calc_coeff.f90
  src/calc_d.f90
  src/calc_mflux.f90
  src/calc_mu_g.f90
  src/calc_outflow.f90
  src/calc_ro_g.f90
  src/calc_trd_g.f90
  src/check_bqend.f90
  src/check_convergence.f90
  src/check_data/check_axis.f90
  src/check_data/check_bc_dem.f90
  src/check_data/check_bc_geometry.f90
  src/check_data/check_bc_inflow.f90
  src/check_data/check_bc_outflow.f90
  src/check_data/check_bc_walls.f90
  src/check_data/check_boundary_conditions.f90
  src/check_data/check_gas_phase.f90
  src/check_data/check_geometry_prereqs.f90
  src/check_data/check_ic_common_discrete.f90
  src/check_data/check_initial_conditions.f90
  src/check_data/check_numerics.f90
  src/check_data/check_output_control.f90
  src/check_data/check_point_sources.f90
  src/check_data/check_run_control.f90
  src/check_data/check_solids_common_all.f90
  src/check_data/check_solids_common_discrete.f90
  src/check_data/check_solids_dem.f90
  src/check_data/check_solids_model_prereqs.f90
  src/check_data/check_solids_phases.f90
  src/check_plane.f90
  src/constant_mod.f90
  src/conv_dif_u_g.f90
  src/conv_dif_v_g.f90
  src/conv_dif_w_g.f90
  src/conv_pp_g.f90
  src/conv_rop.f90
  src/correct_0.f90
  src/deprecated.f90
  src/des/calc_collision_wall_mod.f90
  src/des/calc_drag_des.f90
  src/des/calc_force_dem.f90
  src/des/calc_grad_des.f90
  src/des/calc_pg_grad.f90
  src/des/cffctowall.f90
  src/des/cfnewvalues.f90
  src/des/cfrelvel.f90
  src/des/cfslide.f90
  src/des/comp_mean_fields.f90
  src/des/des_bc_mod.f90
  src/des/des_init_namelist.f90
  src/des/des_time_march.f90
  src/des/desmpi_mod.f90
  src/des/discretelement_mod.f90
  src/des/drag_gp_des.f90
  src/des/drag_gs_des.f90
  src/des/gas_drag.f90
  src/des/make_arrays_des.f90
  src/des/randomno_mod.f90
  src/des/read_particle_input.f90
  src/des/set_phase_index.f90
  src/des/stl_functions_des_mod.f90
  src/des/vtp_mod.f90
  src/des/write_des_data.f90
  src/discretization_mod.f90
  src/display_resid.f90
  src/dmp_modules/compar_mod.f90
  src/drag_mod.f90
  src/eos_mod.f90
  src/error_manager_mod.f90
  src/exit.f90
  src/fill_bc0.f90
  src/fld_constants_mod.f90
  src/flow_to_vel.f90
  src/functions_mod.f90
  src/geometry_mod.f90
  src/get_bc_area.f90
  src/get_data.f90
  src/get_ps.f90
  src/get_tunit.f90
  src/goal_seek_mflux.f90
  src/ic_mod.f90
  src/init_fluid.f90
  src/init_namelist.f90
  src/io/funits_mod.f90
  src/io/open_file.f90
  src/io/open_files.f90
  src/io/out_bin_512.f90
  src/io/out_bin_512i.f90
  src/io/out_bin_512r.f90
  src/io/output_manager.f90
  src/io/write_error.f90
  src/io/write_header.f90
  src/io/write_out0.f90
  src/io/write_table.f90
  src/location_check.f90
  src/machine_mod.f90
  src/matrix_mod.f90
  src/mfix.f90
  src/output_mod.f90
  src/param1_mod.f90
  src/param_mod.f90
  src/parse_line.f90
  src/parse_resid_string.f90
  src/physical_prop.f90
  src/ps_mod.f90
  src/read_namelist.f90
  src/remove_comment.f90
  src/residual_mod.f90
  src/run_mod.f90
  src/scales_mod.f90
  src/set_bc0.f90
  src/set_bc_flow.f90
  src/set_bc_type.f90
  src/set_domain.f90
  src/set_param.f90
  src/set_ps.f90
  src/solve_pp_g.f90
  src/solve_vel_star.f90
  src/solver_params.f90
  src/solvers/leqsol_mod.f90
  src/solvers/matvec.f90
  src/source_pp_g.f90
  src/source_u_g.f90
  src/source_v_g.f90
  src/source_w_g.f90
  src/tau_u_g.f90
  src/tau_v_g.f90
  src/tau_w_g.f90
  src/time_cpu_mod.f90
  src/toleranc_mod.f90
  src/ur_facs_mod.f90
  src/usr/usr0.f90
  src/usr/usr0_des.f90
  src/usr/usr1.f90
  src/usr/usr1_des.f90
  src/usr/usr2.f90
  src/usr/usr2_des.f90
  src/usr/usr3.f90
  src/usr/usr3_des.f90
  src/usr/usr_drag.f90
  src/usr/usr_init_namelist.f90
  src/usr/usr_mod.f90
  src/usr/write_usr0.f90
  src/usr/write_usr1.f90
  src/utilities_mod.f90
  src/vavg_g.f90
  src/xsi_mod.f90
  src/zero_norm_vel.f90
  src/boxlib_to_mfix.F90)
set(F77_sources  )
set(FPP_includes )
set(FPP_sources  )

include(PreprocessAMReXFortran)
preprocess_boxlib_fortran(FPP_out ${FPP_sources})

FILE(GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.f90)
set(local_includes ${CXX_includes} ${FPP_includes})
set(local_sources  ${CXX_sources} ${F90_sources} ${F77_sources} ${FPP_out} AMReX_buildInfo.cpp)

add_executable(mfix src/main.cpp ${USR_OVERRIDES})

add_library(model STATIC ${local_sources} ${local_includes})

if (MPI_CXX_LINK_FLAGS) 
  string( STRIP ${MPI_CXX_LINK_FLAGS} MPI_CXX_LINK_FLAGS)
endif()

target_link_libraries( mfix model ${CCSE_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES} )

# Remove buildInfo.cpp
add_custom_command(TARGET mfix POST_BUILD COMMAND rm AMReX_buildInfo.cpp
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   COMMENT "Removing AMReX_buildInfo.cpp" )

enable_testing()

set(TESTS
  tests/FLD01-x
  tests/FLD01-y
  tests/FLD01-z
  tests/FLD02-x
  tests/FLD02-y
  tests/FLD02-z
  tests/DEM01-x
  tests/DEM01-y
  tests/DEM01-z
  tests/DEM02-x
  tests/DEM02-y
  tests/DEM02-z
  tests/DEM03-x
  tests/DEM03-y
  tests/DEM03-z
  tests/DEM04-x
  tests/DEM04-y
  tests/DEM04-z
  tests/DEM05-x
  tests/DEM05-y
  tests/DEM05-z
  tests/DEM06-y
  tests/DEM06-z
  tests/DEM06-x
  #tests/DEM07
  #benchmarks/DEM01
  #benchmarks/DEM02
  #benchmarks/DEM03
  #benchmarks/DEM04
  #benchmarks/FLD01
  )

set (BUILD_TEST build_test)
add_custom_target(${BUILD_TEST})

foreach( rel_testdir ${TESTS})

  STRING(REPLACE "/" "_" TEST_NAME ${rel_testdir} )

  set(testdir ${CMAKE_CURRENT_SOURCE_DIR}/${rel_testdir})

  FILE(GLOB USR_OVERRIDES ${testdir}/*.f90)

  set(TEST_EXE mfix_${TEST_NAME})

  add_executable(${TEST_EXE} EXCLUDE_FROM_ALL src/main.cpp ${USR_OVERRIDES})

  target_include_directories(${TEST_EXE} PUBLIC ${testdir} ${CMAKE_Fortran_MODULE_DIRECTORY})
  set_target_properties(${TEST_EXE} PROPERTIES Fortran_MODULE_DIRECTORY ${testdir})

  add_test(NAME ${TEST_NAME}
    COMMAND ${CMAKE_COMMAND}
            -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
            -DCMAKE_COMMAND=${CMAKE_COMMAND}
            -DTEST_EXE=${TEST_EXE}
            -DTEST_NAME=${TEST_NAME}
            -Dtestdir=${testdir}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtests.cmake)


  add_dependencies ( ${BUILD_TEST} ${TEST_EXE})

  set_tests_properties ( ${TEST_NAME}
    PROPERTIES ENVIRONMENT "MFIX=${CMAKE_BINARY_DIR}/${TEST_EXE}"
    )

  target_link_libraries(${TEST_EXE} model ${CCSE_LIBRARIES} ${MPI_LIBRARIES})

endforeach( rel_testdir ${TESTS})
