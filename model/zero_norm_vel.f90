module zero_norm_vel_module
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: zero_norm_vel                                           C
!  Purpose: Set the velocity component normal to a wall to zero        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      subroutine zero_norm_vel(slo,shi,lo,hi,u_g,v_g,w_g,flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1   , only: zero
      USE geometry , only: imax2, jmax2, kmax2
      use iso_c_binding, only: c_double, c_int
      USE functions, only: iminus, jminus, kminus

      IMPLICIT NONE

      integer(c_int), intent(in) :: slo(3), shi(3), lo(3), hi(3)

      real(c_double), intent(inout) ::  u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_double), intent(inout) ::  v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_double), intent(inout) ::  w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      integer :: i,j,k

      do k = slo(3),shi(3)
        do j = slo(2),shi(2)
          do i = slo(1),shi(1)

            if (flag(i,j,k,1)<100) then

               if (flag(i,j,k,2) < 1000) u_g(i,j,k) = zero
               if (flag(i,j,k,3) < 1000) v_g(i,j,k) = zero
               if (flag(i,j,k,4) < 1000) w_g(i,j,k) = zero
            else

               U_G(I,J,K) = ZERO
               V_G(I,J,K) = ZERO
               W_G(I,J,K) = ZERO

               if(flag(i,j,k,1) /= 106 .and. flag(i,j,k,1) /= 107) then
                  if(i < imax2) u_g(iminus(i,j,k),j,k) = zero
                  if(j < jmax2) v_g(i,jminus(i,j,k),k) = zero
                  if(k < kmax2) w_g(i,j,kminus(i,j,k)) = zero
               endif
            endif

          end do
        end do
      end do

      end subroutine zero_norm_vel
END module zero_norm_vel_module
