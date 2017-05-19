module conv_rop_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: conv_rop                                                !
!  Purpose: Calculate the face value of density used for calculating   !
!  convection fluxes. Master routine.                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine conv_rop(lo, hi, rop, slo, shi, &
                          u, ulo, uhi, v, vlo, vhi, w, wlo, whi, &
                          ropX, rxlo, rxhi, ropY, rylo, ryhi, ropZ, rzlo, rzhi) & 
                          bind(C, name="conv_rop")

      use param, only: zero

      integer(c_int), intent(in   ) ::   lo(3),  hi(3)
      integer(c_int), intent(in   ) ::  slo(3),  shi(3)
      integer(c_int), intent(in   ) ::  ulo(3),  uhi(3)
      integer(c_int), intent(in   ) ::  vlo(3),  vhi(3)
      integer(c_int), intent(in   ) ::  wlo(3),  whi(3)
      integer(c_int), intent(in   ) :: rxlo(3), rxhi(3)
      integer(c_int), intent(in   ) :: rylo(3), ryhi(3)
      integer(c_int), intent(in   ) :: rzlo(3), rzhi(3)

      real(c_real), intent(in   ) :: u&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      ! macroscopic density (rho_prime)
      real(c_real), intent(in   ) :: rop&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Face value of density (for calculating convective fluxes)
      real(c_real), intent(  out) :: ropX&
         (rxlo(1):rxhi(1),rxlo(2):rxhi(2),rxlo(3):rxhi(3))
      real(c_real), intent(  out) :: ropY&
         (rylo(1):ryhi(1),rylo(2):ryhi(2),rylo(3):ryhi(3))
      real(c_real), intent(  out) :: ropZ&
         (rzlo(1):rzhi(1),rzlo(2):rzhi(2),rzlo(3):rzhi(3))

      integer :: i,j,k

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1

            if (u(i,j,k) >= zero) then
               ropX(i,j,k) = rop(i-1,j,k)
            else
               ropX(i,j,k) = rop(i  ,j,k)
            end if

          end do
        end do
      end do

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)

            if (v(i,j,k) >= zero) then
               ropY(i,j,k) = rop(i,j-1,k)
            else
               ropY(i,j,k) = rop(i,j  ,k)
            end if

          end do
        end do
      end do

      do k = lo(3),hi(3)+1
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            if (w(i,j,k) >= zero) then
               ropZ(i,j,k) = rop(i,j,k-1)
            else
               ropZ(i,j,k) = rop(i,j,k  )
            end if

          end do
        end do
      end do

      end subroutine conv_rop

end module conv_rop_module
