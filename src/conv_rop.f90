module conv_rop_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CONV_rop                                                !
!  Purpose: Calculate the face value of density used for calculating   !
!  convection fluxes. Master routine.                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine conv_rop(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                          u_g, v_g, w_g, rop_g, ropX, ropY, ropZ, &
                          dt, dx, dy, dz, domlo, domhi) &
                          bind(C, name="conv_rop")

      use run, only: discretize

      integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(  out) :: ropX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: ropY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: ropZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz
!---------------------------------------------------------------------//


      if (discretize(1) == 0) then       ! 0 & 1 => first order upwinding
         call conv_rop0 (slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                         rop_g, U_g, V_g, W_g, ropX, ropY, ropZ)
      else
         call conv_rop1 (discretize(1), &
                         slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                         rop_g, u_g, v_g, w_g, &
                         ropX, ropY, ropZ, dt, dx, dy, dz, domlo, domhi)
      end if

!  contains

   end subroutine conv_rop

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CONV_rop                                                !
!  Purpose: Calculate the face value of density used for calculating   !
!  convection fluxes. FOU routine.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine conv_rop0(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                           rop, U, V, W, ropX, ropY, ropZ)

! Modules
!---------------------------------------------------------------------//
      use param1, only: zero


      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

      ! macroscopic density (rho_prime)
      real(c_real), intent(in   ) :: rop&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Velocity components
      real(c_real), intent(in   ) :: u&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      ! Face value of density (for calculating convective fluxes)
      real(c_real), intent(  out) :: ropX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: ropY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: ropZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer :: I, J, K

      do K = lo(3),hi(3)
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)+1

            if (u(i,j,k) >= ZERO) THEN
               ropX(i,j,k) = rop(i-1,j,k)
            else
               ropX(i,j,k) = rop(i  ,j,k)
            end if

          end do
        end do
      end do

      do K = lo(3),hi(3)
        do J = lo(2),hi(2)+1
          do I = lo(1),hi(1)

            if (V(i,j,k) >= ZERO) THEN
               ropY(i,j,k) = rop(i,j-1,k)
            else
               ropY(i,j,k) = rop(i,j  ,k)
            end if

          end do
        end do
      end do

      do K = lo(3),hi(3)+1
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)

            if (W(i,j,k) >= ZERO) THEN
               ropZ(i,j,k) = rop(i,j,k-1)
            else
               ropZ(i,j,k) = rop(i,j,k  )
            end if

          end do
        end do
      end do

      END SUBROUTINE CONV_rop0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: conv_rop1                                               !
!  Purpose: Calculate the face value of density used for calculating   !
!  convection fluxes. HR routine.  Here interpolate the face value of  !
!  density.                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine conv_rop1(DISC, &
                           slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                           rop, u_g, v_g, w_g, &
                           ropX, ropY, ropZ, &
                           dt, dx, dy, dz, domlo, domhi)

! Modules
!---------------------------------------------------------------------//
      use param1, only: one
      use xsi   , only: calc_xsi_x, calc_xsi_y, calc_xsi_z
      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

! Discretization scheme
      integer, INTENT(IN) :: DISC

! macroscopic density (rho_prime)
      real(c_real), INTENT(in) :: rop&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Velocity components
      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Face value of density (for calculating convective fluxes)
      real(c_real), intent(  out) :: ropX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: ropY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: ropZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in) :: dt, dx, dy, dz

      integer :: xhi(3)
      integer :: i,j,k

      real(c_real), allocatable :: xsi_(:,:,:)

      ! Calculate factors

      xhi(1) = hi(1)+1
      xhi(2) = hi(2)
      xhi(3) = hi(3)

      allocate( xsi_(lo(1):xhi(1), lo(2):xhi(2), lo(3):xhi(3)) )
      call calc_xsi_x (DISC, rop, slo, shi, u_g, ulo, uhi, &
         xsi_,  lo, xhi, dt, dx, dy, dz, domlo, domhi)

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
            ropX(i,j,k) = ((one - xsi_(i,j,k))*rop(i-1,j,k) + &
                                  xsi_(i,j,k) *rop(i  ,j,k) )
          end do
        end do
      end do
      deallocate(xsi_)

!---------------------------------------------------------------------//

      xhi(1) = hi(1)
      xhi(2) = hi(2)+1
      xhi(3) = hi(3)

      allocate( xsi_(lo(1):xhi(1), lo(2):xhi(2), lo(3):xhi(3)) )
      call calc_xsi_y (DISC, rop, slo, shi, v_g, vlo, vhi, &
         xsi_,  lo, xhi, dt, dx, dy, dz, domlo, domhi)

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
            ropY(i,j,k) = ((one - xsi_(i,j,k))*rop(i,j-1,k)+&
                                  xsi_(i,j,k) *rop(i,j  ,k))
          end do
        end do
      end do
      deallocate(xsi_)

!---------------------------------------------------------------------//

      xhi(1) = hi(1)
      xhi(2) = hi(2)
      xhi(3) = hi(3)+1

      allocate( xsi_(lo(1):xhi(1), lo(2):xhi(2), lo(3):xhi(3)) )
      call calc_xsi_z (DISC, rop, slo, shi, w_g, wlo, whi, &
         xsi_,  lo, xhi, dt, dx, dy, dz, domlo, domhi)

      do k = lo(3),hi(3)+1
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            ropZ(i,j,k) = ((one - xsi_(i,j,k))*rop(i,j,k-1) + &
                                  xsi_(i,j,k) *rop(i,j,k  ) )
          end do
        end do
      end do
      deallocate(xsi_)

      end subroutine conv_rop1

end module conv_rop_module
