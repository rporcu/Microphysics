!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: run                                                    C
!  Purpose: Common block containing run control data                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE run

! Modules
!---------------------------------------------------------------------//

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use param, only: dim_M, dim_eqs
      use param1, only: UNDEFINED_I

!---------------------------------------------------------------------//


! Main filename to be used for output files  Name must
! still be legal after extensions are added to it.
      CHARACTER(LEN=60) :: RUN_NAME

! Brief description of the problem.
      CHARACTER(LEN=60) :: DESCRIPTION

! Units for data input and output: SI.
      CHARACTER(LEN=16) :: UNITS

! Type of run: NEW, RESTART
      CHARACTER(LEN=16) :: RUN_TYPE

! Variable which triggers automatic restart
      LOGICAL :: AUTOMATIC_RESTART

! counter to keep track of how many auto_retart were performed
      INTEGER :: ITER_RESTART

! version.release of software
      CHARACTER(LEN=10) :: ID_VERSION

! Stop-time of the run.
      real(c_real) :: TSTOP

! Discretization scheme for different equations
      INTEGER :: DISCRETIZE(DIM_EQS)

! If .TRUE. call user-defined subroutines
      LOGICAL :: CALL_USR

! Drag model options (see drag_gs for full details)
! default is syam_obrien (may enforce a corrected Umf by defining
! drag_c1 and drag_d1 accordingly)
      CHARACTER(64) :: DRAG_TYPE
      INTEGER :: DRAG_TYPE_ENUM
      INTEGER,PARAMETER :: SYAM_OBRIEN=0
      INTEGER,PARAMETER :: GIDASPOW=1
      INTEGER,PARAMETER :: GIDASPOW_PCF=2
      INTEGER,PARAMETER :: GIDASPOW_BLEND=3
      INTEGER,PARAMETER :: GIDASPOW_BLEND_PCF=4
      INTEGER,PARAMETER :: WEN_YU=5
      INTEGER,PARAMETER :: WEN_YU_PCF=6
      INTEGER,PARAMETER :: KOCH_HILL=7
      INTEGER,PARAMETER :: KOCH_HILL_PCF=8
      INTEGER,PARAMETER :: BVK=9
      INTEGER,PARAMETER :: HYS=10
      INTEGER,PARAMETER :: USER_DRAG=11

! Single particle drag correlation
      CHARACTER(64) :: CD_FUNCTION


! STOP Trigger mechanism to terminate MFIX normally before batch
! queue terminates flag variable to check for end of batch queue when
! set to TRUE check performed at the beginning of each time step and
! termination of mfix triggered after saving all files if condition
! is met
      LOGICAL :: CHK_BATCHQ_END
! variable to store the total wall clock duration of the batch queue
! session wall clock time specified in seconds
! for jaguarcnl@NCCS max wall clock limit is 2.5 hr limit up to 512
! processors
      real(c_real) :: BATCH_WALLCLOCK
! variable to set a buffer time before the batch queue session ends to
! make sure once MFIX is triggered to shutdown, there is sufficient
! time to save files, make copies to HPSS storage before batch queue
! time runs out. Current logic in MFIX checks for:
!    if CPU_TIME > (BATCH_WALLCLOCK - TERM_BUFFER) then
!    save all .RES .SP files and trigger shutdown
      real(c_real) :: TERM_BUFFER

! parameters for dynamically adjusting time step
! +1 -> increase dt; -1 decrease dt
      INTEGER :: DT_dir = -1

! Maximum Time step.
      real(c_real) :: DT_MAX

! Minimum Time step.
      real(c_real) :: DT_MIN

! Time step adjustment factor (<1.0)
      real(c_real) :: DT_FAC

! If .TRUE. reduce time step when residuals do not decrease
      LOGICAL :: DETECT_STALL

! String which controls reduction of global sums for residual
! calculations
      LOGICAL :: DEBUG_RESID

! Generate log files when negative gas density is detected.
      LOGICAL :: REPORT_NEG_DENSITY


! Flags indicating variable solids density.
      LOGICAL :: SOLVE_ROs(DIM_M), ANY_SOLVE_ROs

! Specifies the type of solids: TFM, DEM, MPPIC
      CHARACTER(len=3), DIMENSION(DIM_M) :: SOLIDS_MODEL

! Flags for various solids phase models.
      LOGICAL :: TFM_SOLIDS
      LOGICAL :: DEM_SOLIDS
      LOGICAL :: PIC_SOLIDS
! The number of the various solids phases.
      INTEGER :: TFM_COUNT = 0
      INTEGER :: DEM_COUNT = 0
      INTEGER :: PIC_COUNT = 0
      logical :: bDist_IO

      END MODULE RUN
