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
      use param, only: UNDEFINED_I

!---------------------------------------------------------------------//


! Main filename to be used for output files  Name must
! still be legal after extensions are added to it.
      CHARACTER(LEN=60) :: RUN_NAME

! Brief description of the problem.
      CHARACTER(LEN=60) :: DESCRIPTION

! Stop-time of the run.
      real(c_real) :: tstop

! Start-time of des-loop
! (updated by des_init_time_loop)
      real(c_real) :: des_tstart

! Discretization scheme for different equations
      integer :: DISCRETIZE(DIM_EQS)

! If .TRUE. call user-defined subroutines
      logical :: CALL_USR

! Maximum Time step.
      real(c_real) :: dt_max

! Minimum Time step.
      real(c_real) :: dt_min

! Time step adjustment factor (<1.0)
      real(c_real) :: dt_fac

! Current time step (super-step)
! (updated by des_init_time_loop)
      real(c_real) :: des_dt

  integer :: nlog
  ! Flag to display messages and residuals on the screen
  logical :: full_log


! Flags for various solids phase models.
      logical :: TFM_SOLIDS
      logical :: DEM_SOLIDS
      logical :: PIC_SOLIDS
! The number of the various solids phases.
      integer :: TFM_COUNT = 0
      integer :: DEM_COUNT = 0
      integer :: PIC_COUNT = 0

! Enable output (via output_manger) on the subdt time scale
! (updated by des_init_time_loop)
      logical :: subdt_io = .false.

      END MODULE RUN
