MODULE GET_DATA_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

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
      SUBROUTINE GET_DATA(time, dt)

      USE check_gas_phase_module, only: check_gas_phase
      USE check_geometry_prereqs_module, only: check_geometry_prereqs
      USE check_numerics_module, only: check_numerics
      USE check_output_control_module, only: check_output_control
      USE check_run_control_module, only: check_run_control
      USE check_solids_model_prereqs_module, only: check_solids_model_prereqs
      USE check_solids_phases_module, only: check_solids_phases
      USE error_manager  , only: init_error_manager
      USE get_bc_area_module, only: get_bc_area
      USE init_namelist_module, only: init_namelist
      USE open_files_mod, only: open_files
      USE read_namelist_module, only: read_namelist
      USE run            , only: run_type, run_name
      USE set_bc_flow_module, only: set_bc_flow
      USE set_parameters_module, only: set_parameters
      USE write_header_module, only: write_header

      use geometry, only: domlo, domhi
      use geometry, only: imax, jmax, kmax

      IMPLICIT NONE

      real(c_real), intent(  out) :: time, dt

! This module call routines to initialize the namelist variables.
      CALL INIT_NAMELIST
! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(time, dt)

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

      domlo(1) = 0
      domlo(2) = 0
      domlo(3) = 0
      domhi(1) = imax-1
      domhi(2) = jmax-1
      domhi(3) = kmax-1

      ! Check the minimum solids phase requirements.
      CALL CHECK_SOLIDS_MODEL_PREREQS

      CALL CHECK_RUN_CONTROL(time, dt)
      CALL CHECK_NUMERICS
      CALL CHECK_OUTPUT_CONTROL

      CALL CHECK_GAS_PHASE
      CALL CHECK_SOLIDS_PHASES
      CALL SET_PARAMETERS

      END SUBROUTINE GET_DATA
END MODULE GET_DATA_MODULE
