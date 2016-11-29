!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_FILTER_DES                                          !
!  Author: J.Musser                                   Date: 25-Nov-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_FILTER_DES

! Runtime Flag: Invoke gas/solids coupled simulation.
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG

      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION :: DXYZ_MIN, DG_DXYZ_MIN


!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_FILTER_DES")

! Calculate reused quanties
      SELECT CASE(DES_INTERP_SCHEME_ENUM)

      CASE(DES_INTERP_GARG)
! Compute the volume of nodes needed in drag_fgs_des0.f
         CALL COMPUTE_VOLUME_OF_NODES
! Setup MPI exchange arrys for nodes
         ! CALL DES_SETNODEINDICES
      END SELECT

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_FILTER_DES
