MODULE GET_DATA_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GET_DATA                                                C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_DATA

      USE check_gas_phase_module, only: check_gas_phase
      USE check_geometry_module, only: check_geometry, check_geometry_des
      USE check_numerics_module, only: check_numerics
      USE check_output_control_module, only: check_output_control
      USE check_run_control_module, only: check_run_control
      USE check_solids_phases_module, only: check_solids_phases
      USE desgrid        , only: desgrid_init
      USE error_manager  , only: init_error_manager
      USE get_bc_area_module, only: get_bc_area
      USE gridmap        , only: gridmap_init
      USE init_namelist_module, only: init_namelist
      USE mpi_init_des   , only: desmpi_init
      USE open_files_mod, only: open_files
      USE read_namelist_module, only: read_namelist
      USE run            , only: run_type, run_name
      USE run, only: dem_solids
      USE set_bc_flow_module, only: set_bc_flow
      USE set_flags_module, only: set_flags
      USE set_geometry_des_module, only: set_geometry_des
      USE set_icbc_flags_module, only: set_icbc_flag
      USE set_max2_module, only: set_max2
      USE stl_preproc_des, only: DES_STL_PREPROCESSING
      USE write_header_module, only: write_header

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! shift DX, DY and DZ values
      LOGICAL, PARAMETER :: SHIFT = .TRUE.

! This module call routines to initialize the namelist variables.
      CALL INIT_NAMELIST
! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(0)

! Initialize the error manager. This call occurs after the mfix.dat
! is read so that message verbosity can be set and the .LOG file
! can be opened.
      CALL INIT_ERROR_MANAGER

! Write header in the .LOG file and to screen.
! Not sure if the values are correct or useful
      CALL WRITE_HEADER

! Open files
      CALL OPEN_FILES(RUN_NAME, RUN_TYPE)

! These checks verify that sufficient information was provided
! to setup the domain indices and DMP gridmap.
      CALL CHECK_GEOMETRY_PREREQS
      CALL CHECK_DMP_PREREQS

! Set up the physical domain indices (cell index max/min values).
      CALL SET_MAX2

! Partition the domain and set indices
      CALL GRIDMAP_INIT

! Check the minimum solids phase requirements.
      CALL CHECK_SOLIDS_MODEL_PREREQS

      CALL CHECK_RUN_CONTROL
      CALL CHECK_NUMERICS
      CALL CHECK_OUTPUT_CONTROL

      CALL CHECK_GAS_PHASE
      CALL CHECK_SOLIDS_PHASES
      CALL SET_PARAMETERS

! Basic geometry checks.
      CALL CHECK_GEOMETRY(SHIFT)
      IF(DEM_SOLIDS) CALL CHECK_GEOMETRY_DES

      END SUBROUTINE GET_DATA
END MODULE GET_DATA_MODULE
