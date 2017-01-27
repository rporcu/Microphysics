MODULE CALC_MFLUX_MODULE

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
      SUBROUTINE CALC_MFLUX (slo, shi, lo, hi, u, v, w, rop_e, rop_n, rop_t, &
         flux_e, flux_n, flux_t, flag, dx, dy, dz) bind(C, name="calc_mflux")

      USE functions, only: iminus, jminus, kminus

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      real(c_real), intent(in) :: dx, dy, dz

      real(c_real), intent(in   ) :: u&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(  out) :: flux_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: flux_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: flux_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))


! Local variables
!---------------------------------------------------------------------//
! Indices
      integer      :: i,j,k
      integer      :: im1, jm1, km1
      real(c_real) :: ayz, axz, axy
!---------------------------------------------------------------------//

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               im1 = i-1
               jm1 = j-1
               km1 = k-1

               ! East face (i+1/2, j, k)
               flux_e(i,j,k) = rop_e(i,j,k)*ayz*u(i,j,k)

               ! West face (i-1/2, j, k)
               if (flag(im1,j,k,1) /= 1) &
                  flux_e(im1,j,k) = rop_e(im1,j,k)*ayz*u(im1,j,k)

               ! North face (i, j+1/2, k)
               flux_n(i,j,k) = rop_n(i,j,k)*axz*v(i,j,k)

               ! South face (i, j-1/2, k)
               if (flag(i,jm1,k,1) /= 1) &
                  flux_n(i,jm1,k) = rop_n(i,jm1,k)*axz*v(i,jm1,k)

               ! Top face (i, j, k+1/2)
               flux_t(i,j,k) = rop_t(i,j,k)*axy*w(i,j,k)

               ! Bottom face (i, j, k-1/2)
               if (flag(i,j,kminus(i,j,k),1) /= 1) &
                  flux_t(i,j,km1) = rop_t(i,j,km1)*axy*w(i,j,km1)

            enddo
         enddo
      enddo

      END SUBROUTINE CALC_MFLUX
END MODULE CALC_MFLUX_MODULE
