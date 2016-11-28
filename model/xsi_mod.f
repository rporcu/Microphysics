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
      IMPLICIT NONE

      CONTAINS


      SUBROUTINE CALC_XSI(DISCR, PHI, U, V, W, xsi_e, xsi_n, xsi_t, incr)

! Modules
!---------------------------------------------------------------------//

      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

      USE discretization, only: phi_c_of
      USE discretization, only: superbee
      USE discretization, only: smart
      USE discretization, only: ultra_quick
      USE discretization, only: quickest
      USE discretization, only: muscl
      USE discretization, only: vanleer
      USE discretization, only: minmod
      USE discretization, only: central_scheme

      USE functions, only: ieast, iwest, jsouth, jnorth, kbot, ktop

      USE geometry, only: do_k
      USE geometry, only: odx, ody, odz

      USE functions, only: im1, ip1, jm1, jp1, km1, kp1

      USE param, only: dimension_3

      USE param1, only: zero

      USE run, only: dt

      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: ival, flush_err_msg
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! discretization method
      INTEGER, INTENT(IN) :: DISCR
! convected quantity
      DOUBLE PRECISION, INTENT(IN) :: phi(istart3:iend3, jstart3:jend3, kstart3:kend3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: V(istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: W(istart3:iend3, jstart3:jend3, kstart3:kend3)
! Convection weighting factors
      DOUBLE PRECISION, INTENT(out) :: xsi_e(istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(out) :: xsi_n(istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(out) :: xsi_t(istart3:iend3, jstart3:jend3, kstart3:kend3)

! shear indicator
      INTEGER, INTENT(IN) :: incr

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IC, ID, IU, JC, JD, JU, KC, KD, KU
      INTEGER :: i, j, k
!
      DOUBLE PRECISION :: PHI_C
! down wind factor
      DOUBLE PRECISION :: dwf
! Courant number
      DOUBLE PRECISION :: cf
! cell widths for QUICKEST
      DOUBLE PRECISION :: oDXc, oDXuc, oDYc, oDYuc, oDZc, oDZuc
!---------------------------------------------------------------------//


       SELECT CASE (DISCR)                    !first order upwinding
       CASE (:1)

       do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3
             XSI_E(i,j,k) = XSI_func(U(i,j,k),ZERO)
             XSI_N(i,j,k) = XSI_func(V(i,j,k),ZERO)
             IF (DO_K) XSI_T(i,j,k) = XSI_func(W(i,j,k),ZERO)
           end do
         end do
       end do

       CASE (2)                               !Superbee

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = ieast(i,j,k)
                IU = iwest(i,j,k)
             ELSE
                IC = ieast(i,j,k)
                ID = I
                IU = ieast(ieast(i,j,k),j,k)
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

             IF (DO_K) THEN
                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = ktop(i,j,k)
                   KU = kbot(i,j,k)
                ELSE
                   KC = ktop(i,j,k)
                   KD = K
                   KU = ktop(i,j,ktop(i,j,k))
                ENDIF

                PHI_C = PHI_C_OF(PHI(i,j,KU),PHI(i,j,KC),PHI(i,j,KD))
                DWF = SUPERBEE(PHI_C)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
             ENDIF
              end do
            end do
          end do

       CASE (3)                               !SMART

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = ieast(i,j,k)
                IU = iwest(i,j,k)
             ELSE
                IC = ieast(i,j,k)
                ID = I
                IU = ieast(ieast(i,j,k),j,k)
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

             IF (DO_K) THEN
                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = ktop(i,j,k)
                   KU = kbot(i,j,k)
                ELSE
                   KC = ktop(i,j,k)
                   KD = K
                   KU = ktop(i,j,ktop(i,j,k))
                ENDIF
                PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
                DWF = SMART(PHI_C)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
             ENDIF
              end do
            end do
          end do

       CASE (4)                               !ULTRA-QUICK

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = ieast(i,j,k)
                IU = iwest(i,j,k)
             ELSE
                IC = ieast(i,j,k)
                ID = i
                IU = ieast(ieast(i,j,k),j,k)
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

             IF (DO_K) THEN
                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = ktop(i,j,k)
                   KU = kbot(i,j,k)
                ELSE
                   KC = ktop(i,j,k)
                   KD = K
                   KU = ktop(i,j,ktop(i,j,k))
                ENDIF
                PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
                CF = ABS(W(i,j,k))*DT*ODZ
                DWF = ULTRA_QUICK(PHI_C,CF)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
             ENDIF
              end do
            end do
          end do


       CASE (5)                               !QUICKEST

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = ieast(i,j,k)
                IU = iwest(i,j,k)
                ODXC = ODX
                ODXUC = ODX
             ELSE
                IC = ieast(i,j,k)
                ID = I
                IU = ieast(ieast(i,j,k),j,k)
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

             IF (DO_K) THEN
                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = ktop(i,j,k)
                   KU = kbot(i,j,k)
                   ODZC = ODZ
                   ODZUC = ODZ
                ELSE
                   KC = ktop(i,j,k)
                   KD = K
                   KU = ktop(i,j,ktop(i,j,k))
                   ODZC = ODZ
                   ODZUC = ODZ
                ENDIF
                PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
                CF = ABS(W(i,j,k))*DT*ODZ
                DWF = QUICKEST(PHI_C,CF,ODZC,ODZUC,ODZ)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
             ENDIF
              end do
            end do
          end do


       CASE (6)                               !MUSCL

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = ieast(i,j,k)
                IU = iwest(i,j,k)
             ELSE
                IC = ieast(i,j,k)
                ID = I
                IU = ieast(ieast(i,j,k),j,k)
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

             IF (DO_K) THEN
                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = ktop(i,j,k)
                   KU = kbot(i,j,k)
                ELSE
                   KC = ktop(i,j,k)
                   KD = K
                   KU = ktop(i,j,ktop(i,j,k))
                ENDIF
                PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
                DWF = MUSCL(PHI_C)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
             ENDIF
              end do
            end do
          end do


       CASE (7)                               !Van Leer

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = ieast(i,j,k)
                IU = iwest(i,j,k)
             ELSE
                IC = ieast(i,j,k)
                ID = I
                IU = ieast(ieast(i,j,k),j,k)
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

             IF (DO_K) THEN
                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = ktop(i,j,k)
                   KU = kbot(i,j,k)
                ELSE
                   KC = ktop(i,j,k)
                   KD = K
                   KU = ktop(i,j,ktop(i,j,k))
                ENDIF
                PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
                DWF = VANLEER(PHI_C)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
             ENDIF
              end do
            end do
          end do


       CASE (8)                               !Minmod

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = ieast(i,j,k)
                IU = iwest(i,j,k)
             ELSE
                IC = ieast(i,j,k)
                ID = I
                IU = ieast(ieast(i,j,k),j,k)
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

             IF (DO_K) THEN
                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = ktop(i,j,k)
                   KU = kbot(i,j,k)
                ELSE
                   KC = ktop(i,j,k)
                   KD = K
                   KU = ktop(i,j,ktop(i,j,k))
                ENDIF

                PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
                DWF = MINMOD(PHI_C)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
             ENDIF
              end do
            end do
          end do

       CASE (9)                               ! Central

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3

             IF (U(i,j,k) >= ZERO) THEN
                IC = I
                ID = ieast(i,j,k)
                IU = iwest(i,j,k)
             ELSE
                IC = ieast(i,j,k)
                ID = I
                IU = ieast(ieast(i,j,k),j,k)
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

             IF (DO_K) THEN
                IF (W(i,j,k) >= ZERO) THEN
                   KC = K
                   KD = ktop(i,j,k)
                   KU = kbot(i,j,k)
                ELSE
                   KC = ktop(i,j,k)
                   KD = K
                   KU = ktop(i,j,ktop(i,j,k))
                ENDIF
                PHI_C = PHI_C_OF(PHI(i,j,ku),PHI(i,j,kc),PHI(i,j,kd))
                DWF = CENTRAL_SCHEME(PHI_C)
                XSI_T(i,j,k) = XSI_func(W(i,j,k),DWF)
             ENDIF
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
      DOUBLE PRECISION FUNCTION XSI_func(XXXv,XXXdwf)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXv, XXXdwf
      XSI_func = (sign(1d0,(-XXXv))+1d0)/(2d0) + &
         sign(1d0,XXXv)*XXXdwf
      END FUNCTION XSI_func

      END MODULE XSI
