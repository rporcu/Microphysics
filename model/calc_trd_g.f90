module calc_trd_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_trD_g                                             C
!  Purpose: Calculate the trace of gas phase rate of strain tensor     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine calc_trd_g(slo,shi,lo,hi,trd_g,u_g, &
         ulo,uhi,v_g,vlo,vhi,w_g,wlo,whi,flag,dx,dy,dz) &
         bind(C, name="calc_trd_g")

      IMPLICIT NONE

      integer(c_int), intent(in ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(c_real), intent(inout) :: trd_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dx, dy, dz

      integer(c_int) :: i, j, k
      real(c_real)   :: odx, ody, odz

      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      ! We only loop over interior "valid" cells here
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            if (flag(i,j,k,1) < 100) then

              TRD_G(i,j,k) = &
                 (U_G(I,J,K)-U_G(i-1,j,k))*odx + &
                 (V_G(I,J,K)-V_G(i,j-1,k))*ody + &
                 (W_G(I,J,K)-W_G(i,j,k-1))*odz

             end if
          end do
        end do
      end do

      end subroutine calc_trd_g
end module calc_trd_g_module
