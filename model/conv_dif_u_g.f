!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_U_g                                            C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_u_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_DIF_U_G(A_M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3, dimension_m
      USE run, only: momentum_x_eq
      USE run, only: discretize
      use fldvar
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!---------------------------------------------------------------------//

      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN

      IF (DISCRETIZE(3) == 0) THEN
         CALL STORE_A_U_G0(A_M)
      ELSE
         CALL STORE_A_U_G1(A_M)
      ENDIF


      RETURN
      END SUBROUTINE CONV_DIF_U_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a u-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GVTERMS(U, V, WW)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE cutcell, only: cut_u_treatment_at
      USE cutcell, only: theta_ue, theta_ue_bar
      USE cutcell, only: theta_u_ne, theta_u_nw
      USE cutcell, only: theta_u_te, theta_u_tw
      USE cutcell, only: alpha_ue_c, alpha_un_c, alpha_ut_c

      USE fldvar, only: u_g, v_g, w_g

      USE fun_avg, only: avg_x_e, avg_x
      USE functions, only: ip_of
      USE geometry, only: do_k
      USE indices, only: i_of, ip1

      USE param, only: dimension_3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(OUT) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: WW(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: IJK, I, IP, IPJK
! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//

!!!$omp parallel do private(IJK,I,IP,IPJK)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         IP = IP1(I)
         IPJK = IP_OF(IJK)

         IF(CUT_U_TREATMENT_AT(IJK)) THEN

! East face (i+1, j, k)
            U(IJK) = (Theta_Ue_bar(IJK) * U_g(IJK) + &
                      Theta_Ue(IJK) * U_g(IPJK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
               alpha_Ue_c(IJK), AW, HW, VELW)
            U(IJK) = U(IJK) * AW

! North face (i+1/2, j+1/2, k)
            V(IJK) = (Theta_U_nw(IJK) * V_g(IJK) + &
                      Theta_U_ne(IJK) * V_g(IPJK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
               ALPHA_Un_c(IJK), AW, HW, VELW)
            V(IJK) = V(IJK) * AW

! Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN
               WW(IJK) = (Theta_U_tw(IJK) * W_g(IJK) + &
                          Theta_U_te(IJK) * W_g(IPJK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
                  ALPHA_Ut_c(IJK), AW, HW, VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF

         ELSE
            U(IJK) = AVG_X_E(U_G(IJK),U_G(IPJK),IP)
            V(IJK) = AVG_X(V_G(IJK),V_G(IPJK),I)
            IF (DO_K) WW(IJK) = AVG_X(W_G(IJK),W_G(IPJK),I)
         ENDIF
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE GET_UCELL_GVTERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the convective fluxes through the faces of a     C
!  u-momentum cell. Note the fluxes are calculated at all faces of     C
!  regardless of flow_at_e of condition of the west, south, or         C
!  bottom face.                                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GCFLUX_TERMS(FLUX_E, FLUX_W, FLUX_N, &
         FLUX_S, FLUX_T, FLUX_B, IJK)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_u_treatment_at
      USE cutcell, only: theta_ue, theta_ue_bar
      USE cutcell, only: theta_u_ne, theta_u_nw
      USE cutcell, only: theta_u_te, theta_u_tw
      USE cutcell, only: alpha_ue_c, alpha_un_c, alpha_ut_c
      USE functions, only: ip_of, im_of, jm_of, km_of
      USE geometry, only: do_k
      USE mflux, only: flux_ge, flux_gn, flux_gt

      USE param1, only: half
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! fluxes through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: flux_e, flux_w
      DOUBLE PRECISION, INTENT(OUT) :: flux_n, flux_s
      DOUBLE PRECISION, INTENT(OUT) :: flux_t, flux_b
! ijk index
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: ipjk, ipjmk, ipjkm

! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//

! indices
      IPJK = IP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IPJMK = IP_OF(IJMK)
      IPJKM = IP_OF(IJKM)

! First calculate the fluxes at the faces
      IF(CUT_U_TREATMENT_AT(IJK)) THEN
! east face: i+1, j, k
         Flux_e = (Theta_Ue_bar(IJK) * Flux_gE(IJK) + &
                   Theta_Ue(IJK) * Flux_gE(IPJK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
            alpha_Ue_c(IJK), AW, HW, VELW)
         Flux_e = Flux_e * AW
! west face: i, j, k
         Flux_w = (Theta_Ue_bar(IMJK) * Flux_gE(IMJK) + &
                   Theta_Ue(IMJK) * Flux_gE(IJK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
            alpha_Ue_c(IMJK), AW, HW, VELW)
         Flux_w = Flux_w * AW


! north face: i+1/2, j+1/2, k
         Flux_n = (Theta_U_nw(IJK) * Flux_gN(IJK) + &
                   Theta_U_ne(IJK) * Flux_gN(IPJK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM', &
            ALPHA_Un_c(IJK), AW, HW, VELW)
         Flux_n = Flux_n * AW
! south face: i+1/2, j-1/2, k
         Flux_s = (Theta_U_nw(IJMK) * Flux_gN(IJMK) + &
                   Theta_U_ne(IJMK) * Flux_gN(IPJMK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
            ALPHA_Un_c(IJMK), AW, HW, VELW)
         Flux_s = Flux_s * AW

         IF (DO_K) THEN
! top face: i+1/2, j, k+1/2
            Flux_t = (Theta_U_tw(IJK) * Flux_gT(IJK) + &
                     Theta_U_te(IJK) * Flux_gT(IPJK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM', &
               ALPHA_Ut_c(IJK), AW, HW, VELW)
            Flux_t = Flux_t * AW
! bottom face: i+1/2, j, k-1/2
            Flux_b = (Theta_U_tw(IJKM) * Flux_gT(IJKM) + &
                      Theta_U_te(IJKM) * Flux_gT(IPJKM))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
               ALPHA_Ut_c(IJKM), AW, HW, VELW)
            Flux_b = Flux_b * AW
         ENDIF
      ELSE
         Flux_e = HALF * (Flux_gE(IJK) + Flux_gE(IPJK))
         Flux_w = HALF * (Flux_gE(IMJK) + Flux_gE(IJK))
         Flux_n = HALF * (Flux_gN(IJK) + Flux_gN(IPJK))
         Flux_s = HALF * (Flux_gN(IJMK) + Flux_gN(IPJMK))

         IF (DO_K) THEN
            Flux_t = HALF * (Flux_gT(IJK) + Flux_gT(IPJK))
            Flux_b = HALF * (Flux_gT(IJKM) + Flux_gT(IPJKM))
         ENDIF
      ENDIF   ! end if/else cut_u_treatment_at
      RETURN
      END SUBROUTINE GET_UCELL_GCFLUX_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a u-momentum cell. Note the fluxes are calculated at       C
!  all faces regardless of flow_at_e condition of the west, south      C
!  or bottom face.                                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GDIFF_TERMS(D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, IJK)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_u_treatment_at
      USE cutcell, only: oneodx_e_u, oneody_n_u, oneodz_t_u

      USE functions, only: wall_at
      USE functions, only: east_of, north_of, top_of
      USE functions, only: south_of, bottom_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: odx, ody_n, odz_t
      USE geometry, only: do_k
      USE geometry, only: ox_e
      USE geometry, only: ayz_u, axz_u, axy_u

      USE indices, only: i_of, j_of, k_of
      USE indices, only: ip1, jm1, km1

      use matrix, only: e, w, s, n, t, b
      USE param1, only: zero


      use fldvar

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! diffusion through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: d_fe, d_fw
      DOUBLE PRECISION, INTENT(OUT) :: d_fn, d_fs
      DOUBLE PRECISION, INTENT(OUT) :: d_ft, d_fb
! ijk index
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: i, j, k, ip, jm, km
      INTEGER :: ijkc, ijke, ijkn, ijkne, ijks, ijkse
      INTEGER :: ijkt, ijkte, ijkb, ijkbe
! length terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
! avg voidage
      DOUBLE PRECISION :: EPGA
!---------------------------------------------------------------------//

      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IP = IP1(I)
      JM = JM1(J)
      KM = KM1(K)

      IJKE = EAST_OF(IJK)
      IF (WALL_AT(IJK)) THEN
         IJKC = IJKE
      ELSE
         IJKC = IJK
      ENDIF
      IJKN = NORTH_OF(IJK)
      IJKNE = EAST_OF(IJKN)
      IJKS = SOUTH_OF(IJK)
      IJKSE = EAST_OF(IJKS)

      IF(CUT_U_TREATMENT_AT(IJK)) THEN
         C_AE = ONEoDX_E_U(IJK)
         C_AW = ONEoDX_E_U(IMJK)
         C_AN = ONEoDY_N_U(IJK)
         C_AS = ONEoDY_N_U(IJMK)
         C_AT = ONEoDZ_T_U(IJK)
         C_AB = ONEoDZ_T_U(IJKM)
      ELSE
         C_AE = ODX(IP)
         C_AW = ODX(I)
         C_AN = ODY_N(J)
         C_AS = ODY_N(JM)
         C_AT = ODZ_T(K)
         C_AB = ODZ_T(KM)
      ENDIF

! East face (i+1, j, k)
      D_FE = MU_G(IJKE)*C_AE*AYZ_U(IJK)
! West face (i, j, k)
      D_FW = MU_G(IJKC)*C_AW*AYZ_U(IMJK)


! North face (i+1/2, j+1/2, k)
      D_FN = AVG_X_H(AVG_Y_H(MU_G(IJKC),MU_G(IJKN),J),&
                     AVG_Y_H(MU_G(IJKE),MU_G(IJKNE),J),I)*&
             C_AN*AXZ_U(IJK)
! South face (i+1/2, j-1/2, k)
      D_FS = AVG_X_H(AVG_Y_H(MU_G(IJKS),MU_G(IJKC),JM),&
                     AVG_Y_H(MU_G(IJKSE),MU_G(IJKE),JM),I)*&
             C_AS*AXZ_U(IJMK)

      D_FT = ZERO
      D_FB = ZERO
      IF (DO_K) THEN
         IJKT = TOP_OF(IJK)
         IJKTE = EAST_OF(IJKT)
         IJKB = BOTTOM_OF(IJK)
         IJKBE = EAST_OF(IJKB)

! Top face (i+1/2, j, k+1/2)
         D_FT = AVG_X_H(AVG_Z_H(MU_G(IJKC),MU_G(IJKT),K),&
                        AVG_Z_H(MU_G(IJKE),MU_G(IJKTE),K),I)*&
                OX_E(I)*C_AT*AXY_U(IJK)
! Bottom face (i+1/2, j, k-1/2)
         D_FB = AVG_X_H(AVG_Z_H(MU_G(IJKB),MU_G(IJKC),KM),&
                        AVG_Z_H(MU_G(IJKBE),MU_G(IJKE),KM),I)*&
                OX_E(I)*C_AB*AXY_U(IJKM)
      ENDIF


      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE GET_UCELL_GDIFF_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_U_g0                                            C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs. C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative. See      C
!  source_u_g.                                                         C
!  Implement FOUP discretization                                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-96  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_U_G0(A_U_G)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE functions, only: flow_at_e
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: do_k

      USE param, only: dimension_3, dimension_m
      USE param1, only: zero
      use matrix, only: e, w, n, s, t, b
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_U_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_g(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
      INTEGER :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Face mass flux
      DOUBLE PRECISION :: flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b
! Diffusion parameter
      DOUBLE PRECISION :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

!---------------------------------------------------------------------//

!$omp     parallel do default(none)                                &
!$omp     private(IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,         &
!$omp             D_fe, d_fw, d_fn, d_fs, d_ft, d_fb,              &
!$omp             flux_e, flux_w, flux_n, flux_s, flux_t, flux_b)  &
!$omp     shared(ijkstart3, ijkend3, do_k, a_u_g)

      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_E(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_UCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk)
            CALL GET_UCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1, j, k)
            IF (Flux_e >= ZERO) THEN
               A_U_G(IJK,E,0) = D_Fe
               A_U_G(IPJK,W,0) = D_Fe + Flux_e
            ELSE
               A_U_G(IJK,E,0) = D_Fe - Flux_e
               A_U_G(IPJK,W,0) = D_Fe
            ENDIF
! West face (i, j, k)
            IF (.NOT.FLOW_AT_E(IMJK)) THEN
               IF (Flux_w >= ZERO) THEN
                  A_U_G(IJK,W,0) = D_Fw + Flux_w
               ELSE
                  A_U_G(IJK,W,0) = D_Fw
               ENDIF
            ENDIF


! North face (i+1/2, j+1/2, k)
            IF (Flux_n >= ZERO) THEN
               A_U_G(IJK,N,0) = D_Fn
               A_U_G(IJPK,S,0) = D_Fn + Flux_n
            ELSE
               A_U_G(IJK,N,0) = D_Fn - Flux_n
               A_U_G(IJPK,S,0) = D_Fn
            ENDIF
! South face (i+1/2, j-1/2, k)
            IF (.NOT.FLOW_AT_E(IJMK)) THEN
               IF (Flux_s >= ZERO) THEN
                  A_U_G(IJK,S,0) = D_Fs + Flux_s
               ELSE
                  A_U_G(IJK,S,0) = D_Fs
               ENDIF
            ENDIF


            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)

! Top face (i+1/2, j, k+1/2)
               IF (Flux_t >= ZERO) THEN
                  A_U_G(IJK,T,0) = D_Ft
                  A_U_G(IJKP,B,0) = D_Ft + Flux_t
               ELSE
                  A_U_G(IJK,T,0) = D_Ft - Flux_t
                  A_U_G(IJKP,B,0) = D_Ft
               ENDIF
! Bottom face (i+1/2, j, k-1/2)
               IF (.NOT.FLOW_AT_E(IJKM)) THEN
                  IF (Flux_b >= ZERO) THEN
                     A_U_G(IJK,B,0) = D_Fb + Flux_b
                  ELSE
                     A_U_G(IJK,B,0) = D_Fb
                  ENDIF
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if (flow_at_e)
      ENDDO   !end do (ijk)
!$omp end parallel do

      RETURN
      END SUBROUTINE STORE_A_U_G0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_U_g1                                            C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative.          C
!  Implements higher order discretization.                             C
!  See source_u_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_U_G1(A_U_G)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE fldvar, only: u_g

      USE functions, only: flow_at_e
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: do_k

      USE param, only: dimension_3, dimension_m
      USE param1, only: one

      use matrix, only: e, w, n, s, t, b

      USE run, only: discretize

      USE tmp_array, only: U => Array1, V => Array2, WW => Array3
      USE tmp_array, only: lock_tmp_array, unlock_tmp_array

      USE xsi, only: calc_xsi
      USE xsi_array, only: xsi_e, xsi_n, xsi_t
      USE xsi_array, only: lock_xsi_array, unlock_xsi_array
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_U_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_g(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IPJK, IMJK, IJPK, IJMK, IJKP, IJKM
! indicator for shear
      INTEGER :: incr
! Diffusion parameter
      DOUBLE PRECISION :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      DOUBLE PRECISION :: Flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b

! temporary use of global arrays:
! array1 (locally u)  - the x directional velocity
!      DOUBLE PRECISION :: U(DIMENSION_3)
! array2 (locally v)  - the y directional velocity
!      DOUBLE PRECISION :: V(DIMENSION_3)
! array3 (locally ww) - the z directional velocity
!      DOUBLE PRECISION :: WW(DIMENSION_3)
!---------------------------------------------------------------------//

      call lock_tmp_array
      call lock_xsi_array

      CALL GET_UCELL_GVTERMS(U, V, WW)

! shear indicator:
      incr=1
      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, &
                     XSI_T, incr)


!!!$omp      parallel do                                                 &
!!!$omp&     private(IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,            &
!!!$omp&             d_fe, d_fw, d_fn, d_fs, d_ft, d_fb,                 &
!!!$omp&             flux_e, flux_w, flux_n, flux_s, flux_t, flux_b)
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_E(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_UCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk)
            CALL GET_UCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)

            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1, j, k)
            A_U_G(IJK,E,0) = D_Fe - XSI_E(IJK) * Flux_e
            A_U_G(IPJK,W,0) = D_Fe + (ONE - XSI_E(IJK)) * Flux_e
! West face (i, j, k)
            IF (.NOT.FLOW_AT_E(IMJK)) THEN
               A_U_G(IJK,W,0) = D_Fw + (ONE - XSI_E(IMJK)) * Flux_w
            ENDIF


! North face (i+1/2, j+1/2, k)
            A_U_G(IJK,N,0) = D_Fn - XSI_N(IJK) * Flux_n
            A_U_G(IJPK,S,0) = D_Fn + (ONE - XSI_N(IJK)) * Flux_n
! South face (i+1/2, j-1/2, k)
            IF (.NOT.FLOW_AT_E(IJMK)) THEN
               A_U_G(IJK,S,0) = D_Fs + (ONE - XSI_N(IJMK)) * Flux_s
            ENDIF


! Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)
               A_U_G(IJK,T,0) = D_Ft - XSI_T(IJK) * Flux_t
               A_U_G(IJKP,B,0) = D_Ft + (ONE - XSI_T(IJK)) * Flux_t
! Bottom face (i+1/2, j, k-1/2)
               IF (.NOT.FLOW_AT_E(IJKM)) THEN
                  A_U_G(IJK,B,0) = D_Fb + (ONE - XSI_T(IJKM)) * Flux_b
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if flow_at_e
      ENDDO   ! end do ijk

      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_U_G1

