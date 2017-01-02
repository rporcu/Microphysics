!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR2_DES(des_pos_new, des_vel_new, omega_new)

      Use discretelement, only: MAX_PIP

      IMPLICIT NONE

      INTEGER :: LL
      double precision, intent(inout) :: des_pos_new(max_pip,3)
      double precision, intent(inout) :: des_vel_new(max_pip,3)
      double precision, intent(inout) :: omega_new(max_pip,3)

! Move particles 63-93 below particles 32-62 to fake a wall.
      DO LL=63,MAX_PIP
         DES_VEL_NEW(LL,:) = 0.0d0
         DES_POS_NEW(LL,1) = DES_POS_NEW(LL-31,1)
         DES_POS_NEW(LL,3) = DES_POS_NEW(LL-31,3)
         DES_POS_NEW(LL,2) = 0.0475d0
         OMEGA_NEW(LL,:) = 0.0d0
      ENDDO

      RETURN
      END SUBROUTINE USR2_DES
