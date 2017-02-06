module calc_trd_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_trD_g                                             !
!  Purpose: Calculate the trace of gas phase rate of strain tensor     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_trd_g(slo, shi, lo, hi, trd_g, &
      u_g, v_g, w_g, dx, dy, dz)

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      real(c_real),   intent(in   ) :: dx, dy, dz

! Dummy Arguments
!-----------------------------------------------
      real(c_real), intent(inout) :: trd_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables
!-----------------------------------------------
      integer :: i, j, k
      real(c_real) :: odx, ody, odz

      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               trd_g(i,j,k) = &
                  (u_g(i,j,k)-u_g(i-1,j,k))*odx + &
                  (v_g(i,j,k)-v_g(i,j-1,k))*ody + &
                  (w_g(i,j,k)-w_g(i,j,k-1))*odz
            end do
         end do
      end do

   end subroutine calc_trd_g
end module calc_trd_g_module
