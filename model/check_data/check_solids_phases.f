!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  SUBROUTINE: CHECK_SOLIDS_PHASES                                     !
!                                                                      !
!  Purpose: Driver routine for calls to solids phase checks.           !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_PHASES


! Global Variables:
!---------------------------------------------------------------------//
! Runtime flag specifying DEM solids
      use run, only: DEM_SOLIDS

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//

!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_PHASES")

! Impose the various model limitations.
      CALL CHECK_SOLIDS_MODEL_LIMITATIONS

! Checks common to all solids models.
      CALL CHECK_SOLIDS_COMMON_ALL

! Checks common to discrete solids phases (DEM).
      IF(DEM_SOLIDS) CALL CHECK_SOLIDS_COMMON_DISCRETE

! Checks specific to the particular solids phase.
      IF(DEM_SOLIDS) CALL CHECK_SOLIDS_DEM

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_PHASES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  SUBROUTINE: CHECK_SOLIDS_MODEL_LIMITATIONS                          !
!                                                                      !
!  Purpose: Impose the limitations of the various solids models. These !
!  checks should be 'high level' in that they only ensure that models  !
!  are only used with phases that support them.                        !
!                                                                      !
!  Author: J.Musser                                   Date: 28-Feb-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_MODEL_LIMITATIONS

! Global Variables:
!---------------------------------------------------------------------//
! Number of solid phases specified by the user/TFM model.
      use physprop, only: SMAX
! Number of discrete solids.
      use discretelement, only: DES_MMAX

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Total number of solids phases.
      INTEGER :: MMAX_TOT
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_MODEL_LIMITATIONS")

! Set up the total number of solids.
      MMAX_TOT = SMAX + DES_MMAX

! The cohesion model is only implemented for DEM simulations

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_SOLIDS_MODEL_LIMITATIONS
