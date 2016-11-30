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

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE des_allocate   , only: des_allocate_arrays
      USE desgrid        , only: desgrid_init
      USE discretelement , only: discrete_element
      USE error_manager  , only: init_error_manager
      USE fldvar         , only: ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g
      USE gridmap        , only: gridmap_init
      USE mpi_init_des   , only: desmpi_init
      USE open_files_mod, only: open_files
      USE param1         , only: undefined
      USE run            , only: run_type, run_name
      USE stl_preproc_des, only: DES_STL_PREPROCESSING

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

! Set constants
      CALL SET_CONSTANTS

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
      IF(DISCRETE_ELEMENT) CALL CHECK_GEOMETRY_DES

! Set grid spacing variables.
      CALL SET_GEOMETRY
      IF(DISCRETE_ELEMENT) CALL SET_GEOMETRY_DES

      CALL CHECK_INITIAL_CONDITIONS
      CALL CHECK_BOUNDARY_CONDITIONS
      CALL CHECK_POINT_SOURCES


!----------------------  DOMAIN SPECIFIC CHECKS  --------------------!


! This call needs to occur before any of the IC/BC checks.
      CALL SET_ICBC_FLAG

! Compute area of boundary surfaces.
      CALL GET_BC_AREA

! Convert (mass, volume) flows to velocities.
      CALL SET_BC_FLOW

! Set the flags for identifying computational cells
      CALL SET_FLAGS

      IF(DISCRETE_ELEMENT) THEN
         CALL DESGRID_INIT
         CALL DESMPI_INIT
         CALL DES_STL_PREPROCESSING
      ENDIF

!--------------------------  ARRAY ALLOCATION -----------------------!

! Allocate array storage.
      CALL ALLOCATE_ARRAYS
      IF(DISCRETE_ELEMENT) CALL DES_ALLOCATE_ARRAYS

      ! Initialize arrays.
      IF(allocated(  EP_G)) EP_G = UNDEFINED
      IF(allocated(  P_G))   P_G = UNDEFINED
      IF(allocated( RO_G))  RO_G = UNDEFINED
      IF(allocated(ROP_G)) ROP_G = UNDEFINED
      IF(allocated(  U_G))   U_G = UNDEFINED
      IF(allocated(  V_G))   V_G = UNDEFINED
      IF(allocated(  W_G))   W_G = UNDEFINED

      IF (DISCRETE_ELEMENT) CALL DES_INIT_ARRAYS

      END SUBROUTINE GET_DATA
