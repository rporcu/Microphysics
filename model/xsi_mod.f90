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

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE CALC_XSI(DISCR, slo, shi, lo, hi, PHI, U, V, W, xsi_e, xsi_n, xsi_t, &
                          dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//

      USE discretization, only: phi_c_of
      USE discretization, only: superbee
      USE discretization, only: smart
      USE discretization, only: ultra_quick
      USE discretization, only: quickest
      USE discretization, only: muscl
      USE discretization, only: vanleer
      USE discretization, only: minmod
      USE discretization, only: central_scheme

      USE functions, only: jsouth, jnorth

      USE geometry , only: domlo, domhi

      USE param1, only: zero

      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: ival, flush_err_msg

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! discretization method
      INTEGER, INTENT(IN) :: DISCR

      ! convected quantity
      real(c_real), INTENT(IN) :: phi&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Velocity components
      real(c_real), INTENT(IN) :: U&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN) :: V&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN) :: W&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Convection weighting factors
      real(c_real), INTENT(out) :: xsi_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(out) :: xsi_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(out) :: xsi_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! shear indicator
      real(c_real), intent(in   ) :: dt, dx, dy, dz
! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IC, ID, IU, JC, JD, JU, KC, KD, KU
      INTEGER :: i, j, k
!
      real(c_real) :: PHI_C
! down wind factor
      real(c_real) :: dwf
! Courant number
      real(c_real) :: cf
! cell widths for QUICKEST
      real(c_real) :: oDXc, oDXuc, oDYc, oDYuc, oDZc, oDZuc
      real(c_real) :: odx, ody, odz
!---------------------------------------------------------------------//

       odx = 1.d0 / dx
       ody = 1.d0 / dy
       odz = 1.d0 / dz

       SELECT CASE (DISCR)                    !first order upwinding
       CASE (:1)

       do k = slo(3),hi(3)
         do j = slo(2),hi(2)
           do i = slo(1),hi(1)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),ZERO)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),ZERO)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),ZERO)
           end do
         end do
       end do

       CASE (2)                               !Superbee

          do k = slo(3),hi(3)
            do j = slo(2),hi(2)
              do i = slo(1),hi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = i
                ID = max(i-1,domlo(1)-1)
                IU = max(i-1,domlo(1)-1)
             ELSE
                IC = i+1
                ID = I
                IU = min(i+2,domhi(1)+1)
                IU = min(i+1,domhi(1)+1)
             ENDIF

             PHI_C = PHI_C_OF(PHI(IU,j,k),PHI(IC,j,k),PHI(ID,j,k))
             DWF = SUPERBEE(PHI_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = jnorth(i,j,k)
                JU = jsouth(i,j,k)
             ELSE
                JC = jnorth(i,j,k)
                JD = J
                JU = jnorth(i,jnorth(i,j,k),k)
             ENDIF

             PHI_C = PHI_C_OF(PHI(i,JU,k),PHI(i,JC,k),PHI(i,JD,k))
             DWF = SUPERBEE(PHI_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF

             PHI_C = PHI_C_OF(PHI(i,j,KU),PHI(i,j,KC),PHI(i,j,KD))
             DWF = SUPERBEE(PHI_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do

       CASE (3)                               !SMART

          do k = slo(3),hi(3)
            do j = slo(2),hi(2)
              do i = slo(1),hi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = i
                IU = i-1
             ELSE
                IC = i+1
                ID = I
                IU = min(i+2,domhi(1)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IU,j,k),PHI(IC,j,k),PHI(ID,j,k))
             DWF = SMART(PHI_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = jnorth(i,j,k)
                JU = jsouth(i,j,k)
             ELSE
                JC = jnorth(i,j,k)
                JD = j
                JU = jnorth(i,jnorth(i,j,k),k)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,ju,k),PHI(i,jc,k),PHI(i,jd,k))
             DWF = SMART(PHI_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
             DWF = SMART(PHI_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do

       CASE (4)                               !ULTRA-QUICK

          do k = slo(3),hi(3)
            do j = slo(2),hi(2)
              do i = slo(1),hi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = i
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = i
                IU = min(i+2,domhi(1)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(iu,j,k),PHI(ic,j,k),PHI(id,j,k))
             CF = ABS(U(i,j,k))*DT*ODX
             DWF = ULTRA_QUICK(PHI_C,CF)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = jnorth(i,j,k)
                JU = jsouth(i,j,k)
             ELSE
                JC = jnorth(i,j,k)
                JD = J
                JU = jnorth(i,jnorth(i,j,k),k)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,ju,k),PHI(i,jc,k),PHI(i,jd,k))
             CF = ABS(V(i,j,k))*DT*ODY
             DWF = ULTRA_QUICK(PHI_C,CF)

             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
             CF = ABS(W(i,j,k))*DT*ODZ
             DWF = ULTRA_QUICK(PHI_C,CF)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do


       CASE (5)                               !QUICKEST

          do k = slo(3),hi(3)
            do j = slo(2),hi(2)
              do i = slo(1),hi(1)

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
             PHI_C = PHI_C_OF(PHI(iu,j,k),PHI(ic,j,k),PHI(id,j,k))
             CF = ABS(U(i,j,k))*DT*ODX
             DWF = QUICKEST(PHI_C,CF,ODXC,ODXUC,ODX)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = jnorth(i,j,k)
                JU = jsouth(i,j,k)
                ODYC = ODY
                ODYUC = ODY
             ELSE
                JC = jnorth(i,j,k)
                JD = J
                JU = jnorth(i,jnorth(i,j,k),k)
                ODYC = ODY
                ODYUC = ODY
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,ju,k),PHI(i,jc,k),PHI(i,jd,k))
             CF = ABS(V(i,j,k))*DT*ODY
             DWF = QUICKEST(PHI_C,CF,ODYC,ODYUC,ODY)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

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
             PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
             CF = ABS(W(i,j,k))*DT*ODZ
             DWF = QUICKEST(PHI_C,CF,ODZC,ODZUC,ODZ)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do


       CASE (6)                               !MUSCL

          do k = slo(3),hi(3)
            do j = slo(2),hi(2)
              do i = slo(1),hi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = i
                IU = min(i+2,domhi(1)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(iu,j,k),PHI(ic,j,k),PHI(id,j,k))
             DWF = MUSCL(PHI_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = jnorth(i,j,k)
                JU = jsouth(i,j,k)
             ELSE
                JC = jnorth(i,j,k)
                JD = J
                JU = jnorth(i,jnorth(i,j,k),k)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,ju,k),PHI(i,jc,k),PHI(i,jd,k))
             DWF = MUSCL(PHI_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KU = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
             DWF = MUSCL(PHI_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do


       CASE (7)                               !Van Leer

          do k = slo(3),hi(3)
            do j = slo(2),hi(2)
              do i = slo(1),hi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = i
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = i
                IU = min(i+2,domhi(1)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(iu,j,k),PHI(ic,j,k),PHI(id,j,k))
             DWF = VANLEER(PHI_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

             IF (V(i,j,k) >= ZERO) THEN
                JC = j
                JD = jnorth(i,j,k)
                JU = jsouth(i,j,k)
             ELSE
                JC = jnorth(i,j,k)
                JD = j
                JU = jnorth(i,jnorth(i,j,k),k)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,ju,k),PHI(i,jc,k),PHI(i,jd,k))
             DWF = VANLEER(PHI_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
             DWF = VANLEER(PHI_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do


       CASE (8)                               !Minmod

          do k = slo(3),hi(3)
            do j = slo(2),hi(2)
              do i = slo(1),hi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = I
                IU = min(i+2,domhi(1)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(iu,j,k),PHI(ic,j,k),PHI(id,j,k))
             DWF = MINMOD(PHI_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = jnorth(i,j,k)
                JU = jsouth(i,j,k)
             ELSE
                JC = jnorth(i,j,k)
                JD = J
                JU = jnorth(i,jnorth(i,j,k),k)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,ju,k),PHI(i,jc,k),PHI(i,jd,k))
             DWF = MINMOD(PHI_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
             DWF = MINMOD(PHI_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do

       CASE (9)                               ! Central

          do k = slo(3),hi(3)
            do j = slo(2),hi(2)
              do i = slo(1),hi(1)

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = i+1
                IU = i-1
             ELSE
                IC = i+1
                ID = I
                IU = min(i+2,domhi(1)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(iu,j,k),PHI(ic,j,k),PHI(id,j,k))
             DWF = CENTRAL_SCHEME(PHI_C)
             XSI_E(i,j,k) = XSI_func(U(i,j,k),DWF)

             IF (V(i,j,k) >= ZERO) THEN
                JC = J
                JD = jnorth(i,j,k)
                JU = jsouth(i,j,k)
             ELSE
                JC = jnorth(i,j,k)
                JD = J
                JU = jnorth(i,jnorth(i,j,k),k)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,ju,k),PHI(i,jc,k),PHI(i,jd,k))
             DWF = CENTRAL_SCHEME(PHI_C)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),DWF)

             IF (W(i,j,k) >= ZERO) THEN
                KC = K
                KD = min(k+1,domhi(3)+1)
                KU = max(k-1,domlo(3)-1)
             ELSE
                KC = min(k+1,domhi(3)+1)
                KD = K
                KU = min(k+2,domhi(3)+1)
             ENDIF
             PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
             DWF = CENTRAL_SCHEME(PHI_C)
             XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
              end do
            end do
          end do

       CASE DEFAULT                           !Error
! should never hit this
          CALL INIT_ERR_MSG("CALC_XSI")
          WRITE(ERR_MSG, 1100) IVAL(DISCR)
 1100 FORMAT('ERROR 1100: Invalid DISCRETIZE= ', A,' The check_data ',&
         'routines should',/, 'have already caught this error and ',&
         'prevented the simulation from ',/,'running. Please notify ',&
         'the MFIX developers.')
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
          CALL FINL_ERR_MSG

       END SELECT

      ! call send_recv(XSI_E,2)
      ! call send_recv(XSI_N,2)
      ! call send_recv(XSI_T,2)

      RETURN
      END SUBROUTINE CALC_XSI

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
      real(c_real), INTENT(IN) :: XXXv, XXXdwf
      XSI_func = (sign(1d0,(-XXXv))+1d0)/(2d0) + &
         sign(1d0,XXXv)*XXXdwf
      END FUNCTION XSI_func

      END MODULE XSI
