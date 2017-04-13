MODULE CHECK_SOLIDS_PHASES_MODULE

   use check_solids_common_all_module, only: check_solids_common_all
   use check_des_solids_module, only: check_solids_dem

   CONTAINS
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
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival

      implicit none

! Local Variables:
!---------------------------------------------------------------------//

!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_PHASES")

! Checks common to all solids models.
      CALL CHECK_SOLIDS_COMMON_ALL

! Checks specific to the particular solids phase.
      IF(DEM_SOLIDS) CALL CHECK_SOLIDS_DEM

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_PHASES
END MODULE CHECK_SOLIDS_PHASES_MODULE
