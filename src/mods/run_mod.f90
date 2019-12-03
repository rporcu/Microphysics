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

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use param, only: dim_M, dim_eqs
      use param, only: UNDEFINED_I

!---------------------------------------------------------------------//

! Brief description of the problem.


! Start-time of des-loop
! (updated by des_init_time_loop)
      real(rt) :: des_tstart

! If .TRUE. call user-defined subroutines


! Current time step (super-step)
! (updated by des_init_time_loop)
      real(rt) :: des_dt

! Flags for various solids phase models.
      logical :: TFM_SOLIDS
      logical :: DEM_SOLIDS
      logical :: PIC_SOLIDS
! The number of the various solids phases.
      integer :: TFM_COUNT = 0
      integer :: DEM_COUNT = 0
      integer :: PIC_COUNT = 0

      END MODULE RUN
