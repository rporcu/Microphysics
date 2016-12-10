!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GET_DATA                                                C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_DATA(ro_sol)

      USE desgrid        , only: desgrid_init
      USE error_manager  , only: init_error_manager
      USE geometry, only: flag
      USE gridmap        , only: gridmap_init
      USE mpi_init_des   , only: desmpi_init
      USE open_files_mod, only: open_files
      USE read_namelist_module, only: read_namelist
      USE run            , only: run_type, run_name
      USE run, only: dem_solids
      USE set_bc_flow_module, only: set_bc_flow
      USE set_geometry_des_module, only: set_geometry_des
      USE set_icbc_flags_module, only: set_icbc_flag
      USE stl_preproc_des, only: DES_STL_PREPROCESSING

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ro_sol

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

! Set grid spacing variables.
      CALL SET_GEOMETRY
      IF(DEM_SOLIDS) CALL SET_GEOMETRY_DES

      CALL CHECK_INITIAL_CONDITIONS
      CALL CHECK_BOUNDARY_CONDITIONS
      CALL CHECK_POINT_SOURCES


!----------------------  DOMAIN SPECIFIC CHECKS  --------------------!


! This call needs to occur before any of the IC/BC checks.
      CALL SET_ICBC_FLAG(flag)

! Compute area of boundary surfaces.
      CALL GET_BC_AREA

! Convert (mass, volume) flows to velocities.
      CALL SET_BC_FLOW

! Set the flags for identifying computational cells
      CALL SET_FLAGS

      IF(DEM_SOLIDS) THEN
         CALL DESGRID_INIT
         CALL DESMPI_INIT(ro_sol)
         CALL DES_STL_PREPROCESSING
      ENDIF

!--------------------------  ARRAY ALLOCATION -----------------------!

      END SUBROUTINE GET_DATA
