!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: XSI                                                    C
!  Author: M. Syamlal                                 Date: 6-MAR-97   C
!                                                                      C
!  Purpose: Determine convection weighting factors for higher order    C
!  discretization.                                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE XSI

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int
      use geometry      , only: domlo, domhi

      USE discretization, only: phi_c_of
      USE discretization, only: superbee
      USE discretization, only: smart
      USE discretization, only: ultra_quick
      USE discretization, only: quickest
      USE discretization, only: muscl
      USE discretization, only: vanleer
      USE discretization, only: minmod
      USE discretization, only: central_scheme

      USE geometry , only: domlo, domhi
      USE param1, only: zero
      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: ival, flush_err_msg

      IMPLICIT NONE

      CONTAINS

      subroutine calc_xsi_e(DISCR, phi, philo, phihi, U, vello, velhi, xsi_e, xlo, xhi, dt, dx, dy, dz)

      integer     , intent(in   ) :: philo(3),phihi(3),vello(3),velhi(3),xlo(3),xhi(3)

      ! discretization method
      integer, intent(IN) :: DISCR

      ! convected quantity
      real(c_real), intent(IN) :: phi&
         (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))

      ! Velocity components
      real(c_real), intent(IN) :: U&
         (vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))

      ! Convection weighting factors
      real(c_real), intent(out) :: xsi_e&
         (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz

!---------------------------------------------------------------------//
! Indices
      integer :: IC, ID, IU
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
!---------------------------------------------------------------------//

       odx = 1.d0 / dx
       ody = 1.d0 / dy
       odz = 1.d0 / dz

       SELECT CASE (DISCR)                    !first order upwinding
       CASE (:1)

       do k = xlo(3),xhi(3)
         do j = xlo(2),xhi(2)
           do i = xlo(1),xhi(1)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),ZERO)
           end do
         end do
       end do

       CASE (2)                               !Superbee

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

                IF (U(i,j,k) >= ZERO) THEN
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
                xsi_e(i,j,k) = xsi_func(u(i,j,k),dwf)

              end do
            end do
          end do

       CASE (3)                               !SMART

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

                IF (U(i,j,k) >= ZERO) THEN
                   IC = I
                   ID = i
                   IU = i-1
                ELSE
                   IC = i+1
                   ID = I
                   IU = min(i+2,domhi(1)+1)
                ENDIF
                phi_C = phi_C_OF(phi(IU,j,k),phi(IC,j,k),phi(ID,j,k))
                DWF = SMART(phi_C)
                XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

              end do
            end do
          end do

       CASE (4)                               !ULTRA-QUICK

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = i
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = i
                IU = min(i+2,domhi(1)+1)
             ENDIF
             phi_C = phi_C_OF(phi(iu,j,k),phi(ic,j,k),phi(id,j,k))
             CF = ABS(U(i,j,k))*DT*ODX
             DWF = ULTRA_QUICK(phi_C,CF)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

              end do
            end do
          end do


       CASE (5)                               !QUICKEST

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = i
                ID = i+1
                IU = i-1
                ODXC = ODX
                ODXUC = ODX
             ELSE
                IU = min(i+1,domhi(1)+1)
                ID = I
                IU = min(i+2,domhi(1)+1)
                ODXC = ODX
                ODXUC = ODX
             ENDIF
             phi_C = phi_C_OF(phi(iu,j,k),phi(ic,j,k),phi(id,j,k))
             CF = ABS(U(i,j,k))*DT*ODX
             DWF = QUICKEST(phi_C,CF,ODXC,ODXUC,ODX)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

              end do
            end do
          end do


       CASE (6)                               !MUSCL

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = i
                IU = min(i+2,domhi(1)+1)
             ENDIF
             phi_C = phi_C_OF(phi(iu,j,k),phi(ic,j,k),phi(id,j,k))
             DWF = MUSCL(phi_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

              end do
            end do
          end do

       CASE (7)                               !Van Leer

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = i
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = i
                IU = min(i+2,domhi(1)+1)
             ENDIF
             phi_C = phi_C_OF(phi(iu,j,k),phi(ic,j,k),phi(id,j,k))
             DWF = VANLEER(phi_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

              end do
            end do
          end do

       CASE (8)                               !Minmod

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = I
                IU = min(i+2,domhi(1)+1)
             ENDIF
             phi_C = phi_C_OF(phi(iu,j,k),phi(ic,j,k),phi(id,j,k))
             DWF = MINMOD(phi_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

              end do
            end do
          end do

       CASE (9)                               ! Central

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = I
                IU = min(i+2,domhi(1)+1)
             ENDIF
             phi_C = phi_C_OF(phi(iu,j,k),phi(ic,j,k),phi(id,j,k))
             DWF = CENTRAL_SCHEME()
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

              end do
            end do
          end do

       END SELECT

      end subroutine calc_xsi_e

!---------------------------------------------------------------------//
!---------------------------------------------------------------------//
!---------------------------------------------------------------------//

      subroutine calc_xsi_n(DISCR, phi, philo, phihi, V, vello, velhi, xsi_n, xlo, xhi, dt, dx, dy, dz)

      integer     , intent(in   ) :: philo(3),phihi(3),vello(3),velhi(3),xlo(3),xhi(3)

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

       SELECT CASE (DISCR)                    !first order upwinding
       CASE (:1)

       do k = xlo(3),xhi(3)
         do j = xlo(2),xhi(2)
           do i = xlo(1),xhi(1)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),ZERO)
           end do
         end do
       end do

       CASE (2)                               !Superbee

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


       CASE (3)                               !SMART

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

                IF (V(i,j,k) >= ZERO) THEN
                   JC = J
                   JD = min(j+1,domhi(2)+1)
                   JU = max(j-1,domlo(2)-1)
                ELSE
                   JC = min(j+1,domhi(2)+1)
                   JD = j
                   JU = min(j+2,domhi(2)+1)
                ENDIF
                phi_C = phi_C_OF(phi(i,ju,k),phi(i,jc,k),phi(i,jd,k))
                DWF = SMART(phi_C)
                XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

              end do
            end do
          end do

       CASE (4)                               !ULTRA-QUICK

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = min(j+1,domhi(2)+1)
                JU = max(j-1,domlo(2)-1)
             ELSE
                JC = min(j+1,domhi(2)+1)
                JD = J
                JU = min(j+2,domhi(2)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,ju,k),phi(i,jc,k),phi(i,jd,k))
             CF = ABS(V(i,j,k))*DT*ODY
             DWF = ULTRA_QUICK(phi_C,CF)

             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

              end do
            end do
          end do


       CASE (5)                               !QUICKEST

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = min(j+1,domhi(2)+1)
                JU = max(j-1,domlo(2)-1)
                ODYC = ODY
                ODYUC = ODY
             ELSE
                JC = min(j+1,domhi(2)+1)
                JD = J
                JU = min(j+2,domhi(2)+1)
                ODYC = ODY
                ODYUC = ODY
             ENDIF
             phi_C = phi_C_OF(phi(i,ju,k),phi(i,jc,k),phi(i,jd,k))
             CF = ABS(V(i,j,k))*DT*ODY
             DWF = QUICKEST(phi_C,CF,ODYC,ODYUC,ODY)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

              end do
            end do
          end do


       CASE (6)                               !MUSCL

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = min(j+1,domhi(2)+1)
                JU = max(j-1,domlo(2)-1)
             ELSE
                JC = min(j+1,domhi(2)+1)
                JD = J
                JU = min(j+2,domhi(2)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,ju,k),phi(i,jc,k),phi(i,jd,k))
             DWF = MUSCL(phi_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

              end do
            end do
          end do


       CASE (7)                               !Van Leer

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (V(i,j,k) >= ZERO) THEN
                JC = j
                JD = min(j+1,domhi(2)+1)
                JU = max(j-1,domlo(2)-1)
             ELSE
                JC = min(j+1,domhi(2)+1)
                JD = j
                JU = min(j+2,domhi(2)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,ju,k),phi(i,jc,k),phi(i,jd,k))
             DWF = VANLEER(phi_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

              end do
            end do
          end do


       CASE (8)                               !Minmod

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = min(j+1,domhi(2)+1)
                JU = max(j-1,domlo(2)-1)
             ELSE
                JC = min(j+1,domhi(2)+1)
                JD = J
                JU = min(j+2,domhi(2)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,ju,k),phi(i,jc,k),phi(i,jd,k))
             DWF = MINMOD(phi_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

              end do
            end do
          end do

       CASE (9)                               ! Central

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = min(j+1,domhi(2)+1)
                JU = max(j-1,domlo(2)-1)
             ELSE
                JC = min(j+1,domhi(2)+1)
                JD = J
                JU = min(j+2,domhi(2)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,ju,k),phi(i,jc,k),phi(i,jd,k))
             DWF = CENTRAL_SCHEME()
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

              end do
            end do
          end do

       END SELECT

      end subroutine calc_xsi_n

!---------------------------------------------------------------------//

      subroutine calc_xsi_t(DISCR, phi, philo, phihi, W, vello, velhi, xsi_t, xlo, xhi, dt, dx, dy, dz)

      integer     , intent(in   ) :: philo(3),phihi(3),vello(3),velhi(3),xlo(3),xhi(3)

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

       SELECT CASE (DISCR)                    !first order upwinding
       CASE (:1)

       do k = xlo(3),xhi(3)
         do j = xlo(2),xhi(2)
           do i = xlo(1),xhi(1)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),ZERO)
           end do
         end do
       end do

       CASE (2)                               !Superbee

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

       CASE (3)                               !SMART

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = min(k+1,domhi(3)+1)
                   KU = max(k-1,domlo(3)-1)
                ELSE
                   KC = min(k+1,domhi(3)+1)
                   KD = K
                   KU = min(k+2,domhi(3)+1)
                ENDIF
                phi_C = phi_C_OF(phi(i,j,ku),phi(i,j,kc),phi(i,j,kd))
                DWF = SMART(phi_C)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)

              end do
            end do
          end do

       CASE (4)                               !ULTRA-QUICK

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,j,ku),phi(i,j,kc),phi(i,j,kd))
             CF = ABS(W(i,j,k))*DT*ODZ
             DWF = ULTRA_QUICK(phi_C,CF)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do


       CASE (5)                               !QUICKEST

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
                ODZC = ODZ
                ODZUC = ODZ
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
                ODZC = ODZ
                ODZUC = ODZ
             ENDIF
             phi_C = phi_C_OF(phi(i,j,ku),phi(i,j,kc),phi(i,j,kd))
             CF = ABS(W(i,j,k))*DT*ODZ
             DWF = QUICKEST(phi_C,CF,ODZC,ODZUC,ODZ)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)

              end do
            end do
          end do


       CASE (6)                               !MUSCL

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KU = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,j,ku),phi(i,j,kc),phi(i,j,kd))
             DWF = MUSCL(phi_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do


       CASE (7)                               !Van Leer

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,j,ku),phi(i,j,kc),phi(i,j,kd))
             DWF = VANLEER(phi_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)

              end do
            end do
          end do


       CASE (8)                               !Minmod

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,j,ku),phi(i,j,kc),phi(i,j,kd))
             DWF = MINMOD(phi_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do

       CASE (9)                               ! Central

          do k = xlo(3),xhi(3)
            do j = xlo(2),xhi(2)
              do i = xlo(1),xhi(1)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             phi_C = phi_C_OF(phi(i,j,ku),phi(i,j,kc),phi(i,j,kd))
             DWF = CENTRAL_SCHEME()
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do

       END SELECT

      END subroutine calc_xsi_t

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
