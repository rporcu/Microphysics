module v_g_conv_dif

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding    , only: c_int
   use param        , only: half, one, zero

   implicit none

   private
   public :: conv_dif_v_g

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: conv_div_v_g(A_m)                                      !
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  !
!  The off-diagonal coefficients calculated here must be positive. The !
!  center coefficient and the source vector are negative;              !
!  See source_v_g                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine conv_dif_v_g(lo, hi, &
      slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
      A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)


      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      integer(c_int), intent(in   ) :: wlo(3),whi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      real(c_real)  , intent(in   ) :: dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

!---------------------------------------------------------------------//

      call store_a_v_g0(lo, hi, &
           slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
           A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)

      end subroutine conv_dif_v_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: store_a_v_g0                                            !
!  Purpose: Determine convection diffusion terms for V_g momentum eqs. !
!  The off-diagonal coefficients calculated here must be positive.     !
!  The center coefficient and the source vector are negative. See      !
!  source_v_g.                                                         !
!                                                                      !
!  Implement FOUP discretization                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine store_a_v_g0(lo, hi, &
         slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
         A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)

      use functions, only: avg_h
      use matrix, only: e, w, n, s, t, b

      implicit none

      integer     , intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer(c_int) :: i, j, k
      real(c_real)   :: lflux_lo, lflux_hi
      real(c_real)   :: ayz_x,axz_y,axy_z

      ayz_x = dy*dz / dx
      axz_y = dx*dz / dy
      axy_z = dx*dy / dz

!---------------------------------------------------------------------//

      ! Diffusion terms
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               A_m(i,j,k,e) = avg_h(avg_h(mu_g(i  ,j-1,k),mu_g(i+1,j-1,k)),&
                                    avg_h(mu_g(i  ,j  ,k),mu_g(i+1,j  ,k))) * ayz_x
               A_m(i,j,k,w) = avg_h(avg_h(mu_g(i-1,j-1,k),mu_g(i  ,j-1,k)),&
                                    avg_h(mu_g(i-1,j  ,k),mu_g(i  ,j  ,k))) * ayz_x

               A_m(i,j,k,n) = mu_g(i,j  ,k) * axz_y
               A_m(i,j,k,s) = mu_g(i,j-1,k) * axz_y
                              
               A_m(i,j,k,t) = avg_h(avg_h(mu_g(i,j-1,k  ),mu_g(i,j-1,k+1)),&
                                    avg_h(mu_g(i,j  ,k  ),mu_g(i,j  ,k+1))) * axy_z
               A_m(i,j,k,b) = avg_h(avg_h(mu_g(i,j-1,k-1),mu_g(i,j-1,k  )),&
                                    avg_h(mu_g(i,j  ,k-1),mu_g(i,j  ,k  ))) * axy_z
                             
            enddo
         enddo
      enddo

      ! Calculate convection-diffusion fluxes through each of the faces
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               lflux_hi = HALF * (fluxX(i+1,j-1,k) + fluxX(i+1,j,k))
               lflux_lo = HALF * (fluxX(i  ,j-1,k) + fluxX(i  ,j,k))

               if (lflux_hi .lt. zero) &
                  A_m(i,j,k,e) = A_m(i,j,k,e) - lflux_hi

               if (lflux_lo .ge. zero) &
                  A_m(i,j,k,w) = A_m(i,j,k,w) + lflux_lo

               ! ******************************************************* 

               lflux_hi = HALF * (fluxY(i,j  ,k) + fluxY(i,j+1,k))
               lflux_lo = HALF * (fluxY(i,j-1,k) + fluxY(i,j  ,k))

               if (lflux_hi .lt. zero) &
                  A_m(i,j,k,n) = A_m(i,j,k,n) - lflux_hi

               if (lflux_lo .ge. zero) &
                  A_m(i,j,k,s) = A_m(i,j,k,s) + lflux_lo

               ! ******************************************************* 

               lflux_hi = HALF * (fluxZ(i,j-1,k+1) + fluxZ(i,j  ,k+1))
               lflux_lo = HALF * (fluxZ(i,j-1,k  ) + fluxZ(i,j  ,k  ))

               if (lflux_hi .lt. zero) &
                  A_m(i,j,k,t) = A_m(i,j,k,t) - lflux_hi
 
               if (lflux_lo .ge. zero) &
                  A_m(i,j,k,b) = A_m(i,j,k,b) + lflux_lo

            enddo
         enddo
      enddo

      end subroutine store_a_v_g0

end module v_g_conv_dif
