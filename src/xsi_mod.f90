!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: XSI                                                    !
!                                                                      !
!  Purpose: Determine convection weighting factors for higher order    !
!  discretization.                                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE XSI

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use discretization, only: phi_c_of
   use discretization, only: superbee
   use discretization, only: smart
   use discretization, only: ultra_quick
   use discretization, only: quickest
   use discretization, only: muscl
   use discretization, only: vanleer
   use discretization, only: minmod
   use discretization, only: central_scheme

   use param1, only: zero

   implicit none

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Procedure: calc_xsi_x                                               !
!                                                                      !
!  Purpose: Determine convection weighting factors for higher order    !
!  discretization in x-axial direction.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_xsi_x(phi, philo, phihi, vel, vello, velhi, &
                         xsi_e, xlo, xhi, is_centered)

      integer     , intent(in   ) :: philo(3),phihi(3)
      integer     , intent(in   ) :: vello(3),velhi(3)
      integer     , intent(in   ) :: xlo(3),xhi(3)
      logical     , intent(in   ) :: is_centered

      ! convected quantity
      real(c_real), intent(in) :: phi&
         (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))

      ! Velocity components
      real(c_real), intent(in) :: vel&
         (vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))

      ! Convection weighting factors
      real(c_real), intent(out) :: xsi_e&
         (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))

      integer(c_int) :: i, j, k
      real(c_real)   :: phi_c, dwf

      do k = xlo(3),xhi(3)
        do j = xlo(2),xhi(2)
          do i = xlo(1),xhi(1)

            if (.not. is_centered) then
               if (vel(i,j,k) >= zero) then
                  phi_c = phi_c_of(phi(i-1,j,k),phi(i,j,k),phi(i+1,j,k))
               else
                  phi_c = phi_c_of(phi(i+2,j,k),phi(i+1,j,k),phi(i,j,k))
               endif
            else
               if (vel(i,j,k) >= zero) then
                  phi_c = phi_c_of(phi(i-2,j,k),phi(i-1,j,k),phi(i,j,k))
               else
                  phi_c = phi_c_of(phi(i+1,j,k),phi(i,j,k),phi(i-1,j,k))
               endif
            endif 

            dwf = superbee(phi_c)
            xsi_e(i,j,k) = xsi_func(vel(i,j,k),dwf)

          end do
        end do
      end do

      end subroutine calc_xsi_x


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Procedure: calc_xsi_y                                               !
!                                                                      !
!  Purpose: Determine convection weighting factors for higher order    !
!  discretization in y-axial direction.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_xsi_y(phi, philo, phihi, vel, vello, velhi, &
                            xsi_n, xlo, xhi, is_centered)

      integer     , intent(in   ) :: philo(3),phihi(3),vello(3),velhi(3),xlo(3),xhi(3)
      logical     , intent(in   ) :: is_centered

      ! convected quantity
      real(c_real), intent(IN) :: phi&
         (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))

      ! Velocity components
      real(c_real), intent(IN) :: vel&
         (vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))

      ! Convection weighting factors
      real(c_real), intent(out) :: xsi_n&
         (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))

      integer(c_int) :: i, j, k
      real(c_real)   :: phi_c, dwf

      do k = xlo(3),xhi(3)
        do j = xlo(2),xhi(2)
          do i = xlo(1),xhi(1)

            if (.not. is_centered) then
               if (vel(i,j,k) >= zero) then
                  phi_c = phi_c_of(phi(i,j-1,k),phi(i,j,k),phi(i,j+1,k))
               else
                  phi_c = phi_c_of(phi(i,j+2,k),phi(i,j+1,k),phi(i,j,k))
               endif
            else
               if (vel(i,j,k) >= zero) then
                  phi_c = phi_c_of(phi(i,j-2,k),phi(i,j-1,k),phi(i,j,k))
               else
                  phi_c = phi_c_of(phi(i,j+1,k),phi(i,j,k),phi(i,j-1,k))
               endif
            endif 

            dwf = superbee(phi_c)
            xsi_n(i,j,k) = xsi_func(vel(i,j,k),dwf)

          end do
        end do
      end do

      end subroutine calc_xsi_y

!---------------------------------------------------------------------//

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Procedure: calc_xsi_z                                               !
!                                                                      !
!  Purpose: Determine convection weighting factors for higher order    !
!  discretization in z-axial direction.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_xsi_z(phi, philo, phihi, vel, vello, velhi, &
                            xsi_t, xlo, xhi, is_centered)

      integer     , intent(in   ) :: philo(3),phihi(3),vello(3),velhi(3),xlo(3),xhi(3)
      logical     , intent(in   ) :: is_centered

      ! convected quantity
      real(c_real), intent(IN) :: phi&
         (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))

      ! Velocity components
      real(c_real), intent(IN) :: vel&
         (vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))

      ! Convection weighting factors
      real(c_real), intent(out) :: xsi_t&
         (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))

      integer(c_int) :: i, j, k
      real(c_real)   :: phi_c, dwf

      do k = xlo(3),xhi(3)
        do j = xlo(2),xhi(2)
          do i = xlo(1),xhi(1)

            if (.not. is_centered) then
               if (vel(i,j,k) >= zero) then
                  phi_c = phi_c_of(phi(i,j,k-1),phi(i,j,k),phi(i,j,k+1))
               else
                  phi_c = phi_c_of(phi(i,j,k+2),phi(i,j,k+1),phi(i,j,k))
               endif
            else
               if (vel(i,j,k) >= zero) then
                  phi_c = phi_c_of(phi(i,j,k-2),phi(i,j,k-1),phi(i,j,k))
               else
                  phi_c = phi_c_of(phi(i,j,k+1),phi(i,j,k),phi(i,j,k-1))
               endif
            endif 

            dwf = superbee(phi_C)
            XSI_T(i,j,k) = XSI_func(vel(i,j,k),DWF)

          end do
        end do
      end do

      end subroutine calc_xsi_z

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: Xsi_func                                                  C
!  Purpose: Special function for xsi that should be similar to:        C
!      xsi(v,dw) = merge( v >= 0, dwf, one-dwf)                        C
!                                                                      C
!  Slight difference when v is exactly zero but may not be             C
!  significant.                                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      real(c_real) FUNCTION XSI_func(XXXv,XXXdwf)

      IMPLICIT NONE
      real(c_real), intent(IN) :: XXXv, XXXdwf
      XSI_func = (sign(1d0, -XXXv ) +1.d0) * 0.5d0 + &
                  sign(1d0,  XXXv )*XXXdwf
      END FUNCTION XSI_func

      END MODULE XSI
