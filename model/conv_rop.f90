module conv_rop_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. Master routine.                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine conv_rop(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
         u_g, v_g, w_g, rop_g, rop_ge, rop_gn, rop_gt, &
         dt, dx, dy, dz) bind(C, name="conv_rop")

! Modules
!---------------------------------------------------------------------//
      USE run, only: discretize

      IMPLICIT NONE

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
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: rop_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: rop_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz
!---------------------------------------------------------------------//


      if (discretize(1) == 0) then       ! 0 & 1 => first order upwinding
         call conv_rop0 (slo, shi, lo, hi, ROP_g, U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT)
      else
         call conv_rop1 (discretize(1), &
                         slo, shi, lo, hi, &
                         rop_g, u_g, v_g, w_g, &
                         rop_ge, rop_gn, rop_gt, dt, dx, dy, dz)
      end if

!  contains

   end subroutine conv_rop 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. FOU routine.                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine conv_rop0(slo, shi, lo, hi, ROP, U, V, W, ROP_E, ROP_N, ROP_T)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! macroscopic density (rho_prime)
      real(c_real), intent(in   ) :: rop&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Velocity components
      real(c_real), intent(in   ) :: u&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Face value of density (for calculating convective fluxes)
      real(c_real), intent(  out) :: rop_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: rop_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: rop_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer :: I, J, K

      do K = lo(3),hi(3)
        do J = lo(2),hi(2)
          do I = slo(1),hi(1)

            ! East face (i+1/2, j, k)
            if (U(i,j,k) >= ZERO) THEN
               ROP_E(i,j,k) = ROP(i,j,k)
            else
               ROP_E(i,j,k) = ROP(i+1,j,k)
            end if

          end do
        end do
      end do

      do K = lo(3),hi(3)
        do J = slo(2),hi(2)
          do I = lo(1),hi(1)

            ! North face (i, j+1/2, k)
            if (V(i,j,k) >= ZERO) THEN
               ROP_N(i,j,k) = ROP(i,j,k)
            else
               ROP_N(i,j,k) = ROP(i,j+1,k)
            end if

          end do
        end do
      end do

      do K = slo(3),hi(3)
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)

            ! Top face (i, j, k+1/2)
            if (W(i,j,k) >= ZERO) THEN
               ROP_T(i,j,k) = ROP(i,j,k)
            else
               ROP_T(i,j,k) = ROP(i,j,k+1)
            end if

          end do
        end do
      end do

      END SUBROUTINE CONV_ROP0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP1                                               C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. HR routine.  Here interpolate the face value of  C
!  density.                                                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP1(DISC, &
                           slo, shi, lo, hi, &
                           rop, u, v, w, &
                           rop_e, rop_n, rop_t, &
                           dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: one
      USE xsi, only: calc_xsi
      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

! Discretization scheme
      INTEGER, INTENT(IN) :: DISC

! macroscopic density (rho_prime)
      real(c_real), INTENT(in) :: rop&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Velocity components
      real(c_real), INTENT(IN   ) :: u&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: v&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: w&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Face value of density (for calculating convective fluxes)
      real(c_real), intent(  out) :: rop_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: rop_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: rop_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in) :: dt, dx, dy, dz
!
! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I,J,K
      Integer :: incr
      real(c_real), allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)

!---------------------------------------------------------------------//

      allocate( xsi_e(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate( xsi_n(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate( xsi_t(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )


! Calculate factors
      incr=0
      CALL CALC_XSI (DISC, slo, shi, hi, &
                     ROP, U, V, W, XSI_E, XSI_N, XSI_T, dt, dx, dy, dz)

      do K = lo(3),hi(3)
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)

            ! East face (i+1/2, j, k)
            ROP_E(i,j,k) = ((ONE - XSI_E(i,j,k))*ROP(i,j,k) + &
                                   XSI_E(i,j,k) *ROP(i+1,j,k) )

            ! West face (i-1/2, j, k)
            ROP_E(i-1,j,k) = ((ONE - XSI_E(i-1,j,k))*ROP(i-1,j,k) + &
                                     XSI_E(i-1,j,k) *ROP(i,j,k) )

            ! North face (i, j+1/2, k)
            ROP_N(i,j,k) = ((ONE - XSI_N(i,j,k))*ROP(i,j,k)+&
                                   XSI_N(i,j,k) *ROP(i,j+1,k))

            ! South face (i, j-1/2, k)
            ROP_N(i,j-1,k) = ((ONE - XSI_N(i,j-1,k))*ROP(i,j-1,k) + &
                                     XSI_N(i,j-1,k) *ROP(i,j,k) )

            ! Top face (i, j, k+1/2)
            ROP_T(i,j,k) = ((ONE - XSI_T(i,j,k))*ROP(i,j,k  ) + &
                                   XSI_T(i,j,k) *ROP(i,j,k+1) )

            ! Bottom face (i, j, k-1/2)
            ROP_T(i,j,k-1) = ((ONE - XSI_T(i,j,k-1))*ROP(i,j,k-1) + &
                                     XSI_T(i,j,k-1) *ROP(i,j,k) )

      end do
      end do
      end do

      deallocate( xsi_e, xsi_n, xsi_t)

      end subroutine conv_rop1

end module conv_rop_module
