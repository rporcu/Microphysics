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

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use param, only: dim_M, dim_eqs
      use param1, only: UNDEFINED_I

!---------------------------------------------------------------------//


! Main filename to be used for output files  Name must
! still be legal after extensions are added to it.
      CHARACTER(LEN=60) :: RUN_NAME

! Brief description of the problem.
      CHARACTER(LEN=60) :: DESCRIPTION

! Type of run: NEW, RESTART
      CHARACTER(LEN=16) :: RUN_TYPE

! Stop-time of the run.
      real(c_real) :: TSTOP

! Discretization scheme for different equations
      integer :: DISCRETIZE(DIM_EQS)

! If .TRUE. call user-defined subroutines
      logical :: CALL_USR

! Maximum Time step.
      real(c_real) :: DT_MAX

! Minimum Time step.
      real(c_real) :: DT_MIN

! Time step adjustment factor (<1.0)
      real(c_real) :: DT_FAC

! If .TRUE. reduce time step when residuals do not decrease
      logical :: DETECT_STALL

! Specifies the type of solids: TFM, DEM, MPPIC
      CHARACTER(len=3), DIMENSION(DIM_M) :: SOLIDS_MODEL

! Flags for various solids phase models.
      logical :: TFM_SOLIDS
      logical :: DEM_SOLIDS
      logical :: PIC_SOLIDS
! The number of the various solids phases.
      integer :: TFM_COUNT = 0
      integer :: DEM_COUNT = 0
      integer :: PIC_COUNT = 0

      END MODULE RUN
