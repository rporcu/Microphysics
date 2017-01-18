module zero_norm_vel_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

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

      subroutine zero_norm_vel(slo,shi,lo,hi,u_g,ulo,uhi,v_g,vlo,vhi,w_g,wlo,whi,flag)

      use ic       , only: CYCL_, CYCP_
      use param1   , only: zero
      use functions, only: iminus, jminus, kminus

      implicit none

      integer(c_int), intent(in) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

      real(c_real), intent(inout) ::  u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) ::  v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) ::  w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      integer(c_int), intent(in) ::  flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      integer :: i,j,k

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            ! Not a wall
            if (flag(i,j,k,1) < 100) then

               ! Test and set the hi-side of the cell
               if (flag(i,j,k,2) < 1000) u_g(i,j,k) = zero
               if (flag(i,j,k,3) < 1000) v_g(i,j,k) = zero
               if (flag(i,j,k,4) < 1000) w_g(i,j,k) = zero

            ! Periodic
            else if (flag(i,j,k,1) /= CYCL_ .and. flag(i,j,k,1) /= CYCP_) then

            ! Wall
            else 

               u_g(i            ,j,k) = zero
               u_g(iminus(i,j,k),j,k) = zero
               v_g(i,j            ,k) = zero
               v_g(i,jminus(i,j,k),k) = zero
               w_g(i,j,k            ) = zero
               w_g(i,j,kminus(i,j,k)) = zero

            endif

          end do
        end do
      end do

      end subroutine zero_norm_vel

end module zero_norm_vel_module
