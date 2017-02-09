module calc_mflux_module

   use bl_fort_module, only : c_real
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
      subroutine calc_mflux (slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
         u_g, v_g, w_g, rop_e, rop_n, rop_t, &
         flux_e, flux_n, flux_t, dx, dy, dz) bind(C, name="calc_mflux")

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      real(c_real),   intent(in   ) :: dx, dy, dz

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: rop_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: rop_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: rop_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(  out) :: flux_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: flux_n&
        (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: flux_t&
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
      print *,"ULO IN CMF ",ulo(:), uhi(:) 

! East face (i+1/2, j, k)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = slo(1), hi(1)
               flux_e(i,j,k) = rop_e(i,j,k)*ayz*u_g(i,j,k)
              if (j.eq.0 .and. k.eq.0) print *,'FLUX COMP ', i, flux_e(i,j,k), rop_e(i,j,k), u_g(i,j,k)
            enddo
         enddo
      enddo

! North face (i, j+1/2, k)
      do k = lo(3), hi(3)
         do j = slo(2), hi(2)
            do i = lo(1), hi(1)
               flux_n(i,j,k) = rop_n(i,j,k)*axz*v_g(i,j,k)
            enddo
         enddo
      enddo

! Top face (i, j, k+1/2)
      do k = slo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               flux_t(i,j,k) = rop_t(i,j,k)*axy*w_g(i,j,k)
            enddo
         enddo
      enddo

      end subroutine calc_mflux
end module calc_mflux_module
