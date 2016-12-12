MODULE SET_FILTER_DES_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_FILTER_DES                                          !
!  Author: J.Musser                                   Date: 25-Nov-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_FILTER_DES(flag)

      use cfassign_module, only: compute_volume_of_nodes
      use error_manager, only: init_err_msg, finl_err_msg
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_SCHEME_ENUM

      IMPLICIT NONE

      INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: FLAG

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_FILTER_DES")

! Calculate reused quanties
      SELECT CASE(DES_INTERP_SCHEME_ENUM)

      CASE(DES_INTERP_GARG)
! Compute the volume of nodes needed in drag_fgs_des0.f
         CALL COMPUTE_VOLUME_OF_NODES(flag)
      END SELECT

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_FILTER_DES
END MODULE SET_FILTER_DES_MODULE
