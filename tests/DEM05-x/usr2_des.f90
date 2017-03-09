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
      SUBROUTINE USR2_DES(max_pip, des_pos_new, des_vel_new, omega_new)

      use amrex_fort_module, only : c_real => amrex_real

      implicit none

      integer     , intent(in   ) :: max_pip
      real(c_real), intent(inout) :: des_pos_new(max_pip,3)
      real(c_real), intent(inout) :: des_vel_new(max_pip,3)
      real(c_real), intent(inout) :: omega_new(max_pip,3)

      integer :: ll

! Move particles 63-93 below particles 32-62 to fake a wall.
      DO LL=63,MAX_PIP
         DES_VEL_NEW(LL,:) = 0.0d0
         DES_POS_NEW(LL,1) = 0.0475d0
         DES_POS_NEW(LL,2) = DES_POS_NEW(LL-31,2)
         DES_POS_NEW(LL,3) = DES_POS_NEW(LL-31,3)
         OMEGA_NEW(LL,:) = 0.0d0
      ENDDO

      END SUBROUTINE USR2_DES
