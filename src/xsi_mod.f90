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
   subroutine calc_xsi_x(DISCR, phi, philo, phihi, vel, vello, velhi, &
      xsi_e, xlo, xhi, dt, dx, dy, dz, domlo, domhi)

      integer     , intent(in   ) :: philo(3),phihi(3)
      integer     , intent(in   ) :: vello(3),velhi(3)
      integer     , intent(in   ) :: xlo(3),xhi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3)

      ! discretization method
      integer, intent(IN) :: DISCR

      ! convected quantity
      real(c_real), intent(in) :: phi&
         (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))

      ! Velocity components
      real(c_real), intent(in) :: vel&
         (vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))

      ! Convection weighting factors
      real(c_real), intent(out) :: xsi_e&
         (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz

!---------------------------------------------------------------------//
! Indices
      integer :: ic, id, iu
      integer :: i, j, k
!
      real(c_real) :: phi_C

      ! down wind factor
      real(c_real) :: dwf

      ! Courant number
      real(c_real) :: cf

      ! cell widths for QUICKEST
      real(c_real) :: odxc, odxuc
      real(c_real) :: odx, ody, odz

      integer :: jm_shift
!---------------------------------------------------------------------//

       odx = 1.d0 / dx
       ody = 1.d0 / dy
       odz = 1.d0 / dz


      do k = xlo(3),xhi(3)
        do j = xlo(2),xhi(2)
          do i = xlo(1),xhi(1)

            IF (vel(i,j,k) >= ZERO) THEN
               IC = i
               ID = i+1
               IU = i-1
            ELSE
               IC = i+1
               ID = I
               IU = i+2
            ENDIF
            phi_c = phi_c_of(phi(iu,j,k),phi(ic,j,k),phi(id,j,k))
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
      subroutine calc_xsi_y(DISCR, phi, philo, phihi, V, vello, velhi, xsi_n, xlo, xhi, dt, dx, dy, dz, domlo, domhi)

      integer     , intent(in   ) :: philo(3),phihi(3),vello(3),velhi(3),xlo(3),xhi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3)

      ! discretization method
      integer, intent(IN) :: DISCR

      ! convected quantity
      real(c_real), intent(IN) :: phi&
         (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))

      ! Velocity components
      real(c_real), intent(IN) :: V&
         (vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))

      ! Convection weighting factors
      real(c_real), intent(out) :: xsi_n&
         (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz

      integer :: JC, JD, JU
      integer :: i, j, k

      real(c_real) :: phi_C

      ! down wind factor
      real(c_real) :: dwf

      ! Courant number
      real(c_real) :: cf

      ! cell widths for QUICKEST
      real(c_real) :: odyc, odyuc
      real(c_real) :: odx, ody, odz
!---------------------------------------------------------------------//

      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz

      do k = xlo(3),xhi(3)
        do j = xlo(2),xhi(2)
          do i = xlo(1),xhi(1)

            if (v(i,j,k) >= zero) then
               ju = j-1
               jc = j
               jd = j+1
            else
               ju = j+2
               jc = j+1
               jd = j
            endif

            phi_c = phi_c_of(phi(i,ju,k),phi(i,jc,k),phi(i,jd,k))
            dwf = superbee(phi_c)
            xsi_n(i,j,k) = xsi_func(v(i,j,k),dwf)

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
      subroutine calc_xsi_z(DISCR, phi, philo, phihi, W, vello, velhi, xsi_t, xlo, xhi, dt, dx, dy, dz, domlo, domhi)

      integer     , intent(in   ) :: philo(3),phihi(3),vello(3),velhi(3),xlo(3),xhi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3)

      ! discretization method
      integer, intent(IN) :: DISCR

      ! convected quantity
      real(c_real), intent(IN) :: phi&
         (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))

      ! Velocity components
      real(c_real), intent(IN) :: W&
         (vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))

      ! Convection weighting factors
      real(c_real), intent(out) :: xsi_t&
         (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz
! Local variables
!---------------------------------------------------------------------//
! Indices
      integer :: KC, KD, KU
      integer :: i, j, k
!
      real(c_real) :: phi_C
! down wind factor
      real(c_real) :: dwf
! Courant number
      real(c_real) :: cf
! cell widths for QUICKEST
      real(c_real) :: odzc, odzuc
      real(c_real) :: odx, ody, odz
!---------------------------------------------------------------------//

       odx = 1.d0 / dx
       ody = 1.d0 / dy
       odz = 1.d0 / dz

      do k = xlo(3),xhi(3)
        do j = xlo(2),xhi(2)
          do i = xlo(1),xhi(1)

            if (w(i,j,k) >= zero) then
               ku = k-1
               kc = k
               kd = k+1
            else
               ku = k+2
               kc = k+1
               kd = k
            endif

            phi_C = phi_C_OF(phi(i,j,KU),phi(i,j,KC),phi(i,j,KD))
            DWF = SUPERBEE(phi_C)
            XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)

          end do
        end do
      end do

      END subroutine calc_xsi_z

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
