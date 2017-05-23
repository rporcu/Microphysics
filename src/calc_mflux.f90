module calc_mflux_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MFLUX                                              C
!  Purpose: Calculate the convection fluxes. Master routine.           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine calc_mflux (ulo, uhi, vlo, vhi, wlo, whi, &
         u_g, v_g, w_g, ropX, ropY, ropZ, &
         fluxX, fluxY, fluxZ, dx, dy, dz) bind(C, name="calc_mflux")

      implicit none

      integer(c_int), intent(in   ) :: ulo(3),uhi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      integer(c_int), intent(in   ) :: wlo(3),whi(3)
      real(c_real),   intent(in   ) :: dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: ropX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: ropY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: ropZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(  out) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: fluxY&
        (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: fluxZ&
        (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))


! Local variables
!---------------------------------------------------------------------//
! Indices
      integer      :: i,j,k
      real(c_real) :: ayz, axz, axy
!---------------------------------------------------------------------//

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      do k = ulo(3), uhi(3)
         do j = ulo(2), uhi(2)
            do i = ulo(1), uhi(1)
               fluxX(i,j,k) = ropX(i,j,k)*ayz*u_g(i,j,k)
            enddo
         enddo
      enddo

      do k = vlo(3), vhi(3)
         do j = vlo(2), vhi(2)
            do i = vlo(1), vhi(1)
               fluxY(i,j,k) = ropY(i,j,k)*axz*v_g(i,j,k)
            enddo
         enddo
      enddo

      do k = wlo(3), whi(3)
         do j = wlo(2), whi(2)
            do i = wlo(1), whi(1)
               fluxZ(i,j,k) = ropZ(i,j,k)*axy*w_g(i,j,k)
            enddo
         enddo
      enddo

      end subroutine calc_mflux

end module calc_mflux_module
