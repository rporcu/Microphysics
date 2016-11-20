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


      SUBROUTINE CALC_XSI(DISCR, PHI, U, V, W, XSI_E, XSI_N, XSI_T, &
                          incr)

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

      USE functions, only: east_of, west_of, north_of, south_of
      USE functions, only: top_of, bottom_of
      USE functions, only: ieast, iwest, jsouth, jnorth, kbot, ktop
      USE functions, only: funijk

      USE geometry, only: do_k
      USE geometry, only: odx, odx_e, ody, ody_n, odz, odz_t
      USE geometry, only: ox

      USE indices, only: im1, ip1, jm1, jp1, km1, kp1

      USE param, only: dimension_3

      USE param1, only: zero

      USE run, only: dt

      USE sendrecv, only: send_recv

      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: ival, flush_err_msg
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! discretization method
      INTEGER, INTENT(IN) :: DISCR
! convected quantity
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Convection weighting factors
      DOUBLE PRECISION, INTENT(OUT) :: XSI_e(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: XSI_n(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: XSI_t(DIMENSION_3)
! shear indicator
      INTEGER, INTENT(IN) :: incr

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IJKC, IJKD, IJKU, I, J, K
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
             ijk = funijk(i,j,k)
             XSI_E(IJK) = XSI_func(U(IJK),ZERO)
             XSI_N(IJK) = XSI_func(V(IJK),ZERO)
             IF (DO_K) XSI_T(IJK) = XSI_func(W(IJK),ZERO)
           end do
         end do
       end do

       CASE (2)                               !Superbee

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3
                ijk = funijk(i,j,k)

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = FUNIJK(ieast(i,j,k),j,k)
                IJKU = FUNIJK(iwest(i,j,k),j,k)
             ELSE
                IJKC = FUNIJK(ieast(i,j,k),j,k)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = SUPERBEE(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = FUNIJK(i,jnorth(i,j,k),k)
                IJKU = FUNIJK(i,jsouth(i,j,k),k)
             ELSE
                IJKC = FUNIJK(i,jnorth(i,j,k),k)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = SUPERBEE(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = FUNIJK(i,j,ktop(i,j,k))
                   IJKU = FUNIJK(i,j,kbot(i,j,k))
                ELSE
                   IJKC = FUNIJK(i,j,ktop(i,j,k))
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = SUPERBEE(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
              end do
            end do
          end do

       CASE (3)                               !SMART

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3
                ijk = funijk(i,j,k)

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = FUNIJK(ieast(i,j,k),j,k)
                IJKU = FUNIJK(iwest(i,j,k),j,k)
             ELSE
                IJKC = FUNIJK(ieast(i,j,k),j,k)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = SMART(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = FUNIJK(i,jnorth(i,j,k),k)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = FUNIJK(i,jnorth(i,j,k),k)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = SMART(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = FUNIJK(i,j,ktop(i,j,k))
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = FUNIJK(i,j,ktop(i,j,k))
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = SMART(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
              end do
            end do
          end do

       CASE (4)                               !ULTRA-QUICK

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3
                ijk = funijk(i,j,k)

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             CF = ABS(U(IJK))*DT*ODX_E(I)
             DWF = ULTRA_QUICK(PHI_C,CF)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             CF = ABS(V(IJK))*DT*ODY_N(J)
             DWF = ULTRA_QUICK(PHI_C,CF)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                CF = ABS(W(IJK))*DT*OX(I)*ODZ_T(K)
                DWF = ULTRA_QUICK(PHI_C,CF)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
              end do
            end do
          end do


       CASE (5)                               !QUICKEST

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3
                ijk = funijk(i,j,k)

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
                ODXC = ODX(I)
                ODXUC = ODX_E(IM1(I))
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
                ODXC = ODX(IP1(I))
                ODXUC = ODX_E(IP1(I))
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             CF = ABS(U(IJK))*DT*ODX_E(I)
             DWF = QUICKEST(PHI_C,CF,ODXC,ODXUC,ODX_E(I))
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
                ODYC = ODY(J)
                ODYUC = ODY_N(JM1(J))
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
                ODYC = ODY(JP1(J))
                ODYUC = ODY_N(JP1(J))
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             CF = ABS(V(IJK))*DT*ODY_N(J)
             DWF = QUICKEST(PHI_C,CF,ODYC,ODYUC,ODY_N(J))
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                   ODZC = ODZ(K)
                   ODZUC = ODZ_T(KM1(K))
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                   ODZC = ODZ(KP1(K))
                   ODZUC = ODZ_T(KP1(K))
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                CF = ABS(W(IJK))*DT*OX(I)*ODZ_T(K)
                DWF = QUICKEST(PHI_C,CF,ODZC,ODZUC,ODZ_T(K))
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
              end do
            end do
          end do


       CASE (6)                               !MUSCL

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3
                ijk = funijk(i,j,k)

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = MUSCL(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = MUSCL(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = MUSCL(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
              end do
            end do
          end do


       CASE (7)                               !Van Leer

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3
                ijk = funijk(i,j,k)

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = VANLEER(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = VANLEER(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = VANLEER(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
              end do
            end do
          end do


       CASE (8)                               !Minmod

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3
                ijk = funijk(i,j,k)

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = MINMOD(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = MINMOD(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF

                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = MINMOD(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
              end do
            end do
          end do

       CASE (9)                               ! Central

          do k = kstart3, kend3
            do j = jstart3, jend3
              do i = istart3, iend3
                ijk = funijk(i,j,k)

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = CENTRAL_SCHEME(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = CENTRAL_SCHEME(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = CENTRAL_SCHEME(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
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

      call send_recv(XSI_E,2)
      call send_recv(XSI_N,2)
      call send_recv(XSI_T,2)

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
