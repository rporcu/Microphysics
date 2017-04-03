module conv_rop_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_rop                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. Master routine.                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine conv_rop(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                          u_g, v_g, w_g, rop_g, rop_ge, rop_gn, rop_gt, &
                          dt, dx, dy, dz) &
                          bind(C, name="conv_rop")

      use run, only: discretize

      integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(  out) :: rop_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: rop_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: rop_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz
!---------------------------------------------------------------------//


      if (discretize(1) == 0) then       ! 0 & 1 => first order upwinding
         call conv_rop0 (slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                         rop_g, U_g, V_g, W_g, rop_gE, rop_gN, rop_gT)
      else
         call conv_rop1 (discretize(1), &
                         slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                         rop_g, u_g, v_g, w_g, &
                         rop_ge, rop_gn, rop_gt, dt, dx, dy, dz)
      end if

!  contains

   end subroutine conv_rop

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_rop                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. FOU routine.                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine conv_rop0(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                           rop, U, V, W, rop_E, rop_N, rop_T)

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
      real(c_real), intent(  out) :: rop_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: rop_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: rop_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer :: I, J, K

      do K = lo(3),hi(3)
        do J = lo(2),hi(2)
          do I = slo(1),hi(1)

            ! East face (i+1/2, j, k)
            if (U(i,j,k) >= ZERO) THEN
               rop_E(i,j,k) = rop(i,j,k)
            else
               rop_E(i,j,k) = rop(i+1,j,k)
            end if

          end do
        end do
      end do

      do K = lo(3),hi(3)
        do J = slo(2),hi(2)
          do I = lo(1),hi(1)

            ! North face (i, j+1/2, k)
            if (V(i,j,k) >= ZERO) THEN
               rop_N(i,j,k) = rop(i,j,k)
            else
               rop_N(i,j,k) = rop(i,j+1,k)
            end if

          end do
        end do
      end do

      do K = slo(3),hi(3)
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)

            ! Top face (i, j, k+1/2)
            if (W(i,j,k) >= ZERO) THEN
               rop_T(i,j,k) = rop(i,j,k)
            else
               rop_T(i,j,k) = rop(i,j,k+1)
            end if

          end do
        end do
      end do

      END SUBROUTINE CONV_rop0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: conv_rop1                                               C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. HR routine.  Here interpolate the face value of  C
!  density.                                                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine conv_rop1(DISC, &
                           slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
                           rop, u_g, v_g, w_g, &
                           rop_e, rop_n, rop_t, &
                           dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use param1, only: one
      use xsi   , only: calc_xsi_e, calc_xsi_n, calc_xsi_t
      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

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
      real(c_real), intent(  out) :: rop_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(  out) :: rop_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(  out) :: rop_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in) :: dt, dx, dy, dz

      integer :: xlo(3)
      integer :: i,j,k

      real(c_real), allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)

      ! Calculate factors

      xlo(1) = lo(1)-1
      xlo(2) = lo(2)
      xlo(3) = lo(3)
      allocate( xsi_e(xlo(1): hi(1),xlo(2): hi(2),xlo(3): hi(3)) )
      call calc_xsi_e (DISC, rop, slo, shi, u_g, ulo, uhi, xsi_e, xlo,  hi, &
                       dt, dx, dy, dz)

      xlo(1) = lo(1)
      xlo(2) = lo(2)-1
      xlo(3) = lo(3)
      allocate( xsi_n(xlo(1): hi(1),xlo(2): hi(2),xlo(3): hi(3)) )
      call calc_xsi_n (DISC, rop, slo, shi, v_g, vlo, vhi, xsi_n, xlo,  hi, &
                       dt, dx, dy, dz)

      xlo(1) = lo(1)
      xlo(2) = lo(2)
      xlo(3) = lo(3)-1
      allocate( xsi_t(xlo(1): hi(1),xlo(2): hi(2),xlo(3): hi(3)) )
      call calc_xsi_t (DISC, rop, slo, shi, w_g, wlo, whi, xsi_t, xlo,  hi, &
                       dt, dx, dy, dz)

      ! East face (i+1/2, j, k)
      do K = lo(3),hi(3)
        do J = lo(2),hi(2)
          do I = lo(1)-1,hi(1)
            rop_e(i,j,k) = ((ONE - XSI_E(i,j,k))*rop(i,j,k) + &
                                   XSI_E(i,j,k) *rop(i+1,j,k) )
          end do
        end do
      end do

      ! North face (i, j+1/2, k)
      do K = lo(3),hi(3)
        do J = lo(2)-1,hi(2)
          do I = lo(1),hi(1)
            rop_n(i,j,k) = ((ONE - XSI_N(i,j,k))*rop(i,j,k)+&
                                   XSI_N(i,j,k) *rop(i,j+1,k))
          end do
        end do
      end do

      ! Top face (i, j, k+1/2)
      do K = lo(3)-1,hi(3)
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)
            rop_t(i,j,k) = ((ONE - XSI_T(i,j,k))*rop(i,j,k  ) + &
                                   XSI_T(i,j,k) *rop(i,j,k+1) )
          end do
        end do
      end do

      deallocate( xsi_e, xsi_n, xsi_t)

      end subroutine conv_rop1

end module conv_rop_module
