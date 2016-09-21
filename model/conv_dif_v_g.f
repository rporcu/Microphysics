!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_V_g(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_DIF_V_G(A_M, B_M, IER)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3
      USE run, only: momentum_y_eq
      USE run, only: discretize
      use fldvar
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!---------------------------------------------------------------------//

      IF (.NOT.MOMENTUM_Y_EQ(0)) RETURN

! DO NOT USE DEFERRED CORRECTION TO SOLVE V_G
      IF (DISCRETIZE(4) == 0) THEN               ! 0 & 1 => FOUP
         CALL STORE_A_V_G0(A_M)
      ELSE
         CALL STORE_A_V_G1(A_M)
      ENDIF

      RETURN
      END SUBROUTINE CONV_DIF_V_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a v-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_VCELL_GVTERMS(U, V, WW)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE cutcell, only: cut_v_treatment_at
      USE cutcell, only: theta_vn, theta_vn_bar
      USE cutcell, only: theta_v_ne, theta_v_se
      USE cutcell, only: theta_v_nt, theta_v_st
      USE cutcell, only: alpha_ve_c, alpha_vn_c, alpha_vt_c

      USE fldvar, only: u_g, v_g, w_g

      USE fun_avg, only: avg_y_n, avg_y
      USE functions, only: jp_of
      USE geometry, only: do_k
      USE indices, only: j_of

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
      INTEGER :: IJK, J, IJPK
! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//

!!!$omp parallel do private(IJK,J,IJPK)
      DO IJK = ijkstart3, ijkend3
         J = J_OF(IJK)
         IJPK = JP_OF(IJK)

         IF(CUT_V_TREATMENT_AT(IJK)) THEN

! East face (i+1/2, j+1/2, k)
            U(IJK) = (Theta_V_se(IJK) * U_g(IJK) + &
                      Theta_V_ne(IJK) * U_g(IJPK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM', &
               ALPHA_Ve_c(IJK), AW, HW, VELW)
            U(IJK) = U(IJK) * AW

! North face (i, j+1, k)
            V(IJK) = (Theta_Vn_bar(IJK) * V_g(IJK) + &
                      Theta_Vn(IJK) * V_g(IJPK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
               alpha_Vn_c(IJK), AW, HW, VELW)
            V(IJK) = V(IJK) * AW

! Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN
               WW(IJK) = (Theta_V_nt(IJK) * W_g(IJK) + &
                          Theta_V_st(IJK) * W_g(IJPK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',&
                  ALPHA_Vt_c(IJK), AW, HW, VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF

         ELSE
            U(IJK) = AVG_Y(U_G(IJK),U_G(IJPK),J)
            V(IJK) = AVG_Y_N(V_G(IJK),V_G(IJPK),0)
            IF (DO_K) WW(IJK) = AVG_Y(W_G(IJK),W_G(IJPK),J)
         ENDIF
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE GET_VCELL_GVTERMS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the convective fluxes through the faces of a     C
!  v-momentum cell. Note the fluxes are calculated at all faces        C
!  regardless of flow_at_n of condition of the west, south, or         C
!  bottom face.                                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_VCELL_GCFLUX_TERMS(FLUX_E, FLUX_W, FLUX_N, &
         FLUX_S, FLUX_T, FLUX_B, IJK)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_v_treatment_at
      USE cutcell, only: theta_vn, theta_vn_bar
      USE cutcell, only: theta_v_ne, theta_v_se
      USE cutcell, only: theta_v_nt, theta_v_st
      USE cutcell, only: alpha_ve_c, alpha_vn_c, alpha_vt_c

      USE functions, only: jp_of, im_of, jm_of, km_of

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
! indices
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: ijpk, imjpk, ijpkm

! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//
! indices
      IJPK = JP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IJPKM = JP_OF(IJKM)
      IMJPK = JP_OF(IMJK)

      IF(CUT_V_TREATMENT_AT(IJK)) THEN
! East face (i+1/2, j+1/2, k)
         Flux_e = (Theta_V_se(IJK) * Flux_gE(IJK) + &
                 Theta_V_ne(IJK) * Flux_gE(IJPK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',&
            ALPHA_Ve_c(IJK), AW, HW, VELW)
! West face (i-1/2, j+1/2, k)
         Flux_e = Flux_e * AW
         Flux_w = (Theta_V_se(IMJK) * Flux_gE(IMJK) + &
                 Theta_V_ne(IMJK) * Flux_gE(IMJPK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',&
            ALPHA_Ve_c(IMJK), AW, HW, VELW)
         Flux_w = Flux_w * AW

! North face (i, j+1, k)
         Flux_n = (Theta_Vn_bar(IJK) * Flux_gN(IJK) + &
                 Theta_Vn(IJK) * Flux_gN(IJPK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',&
            alpha_Vn_c(IJK), AW, HW, VELW)
         Flux_n = Flux_n * AW
! South face (i, j, k)
         Flux_s = (Theta_Vn_bar(IJMK) * Flux_gN(IJMK) + &
                 Theta_Vn(IJMK) * Flux_gN(IJK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',&
            alpha_Vn_c(IJMK), AW, HW, VELW)
         Flux_s = Flux_s * AW

         IF (DO_K) THEN
! Top face (i, j+1/2, k+1/2)
            Flux_t = (Theta_V_nt(IJK) * Flux_gT(IJK) + &
                    Theta_V_st(IJK) * Flux_gT(IJPK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',&
               ALPHA_Vt_c(IJK), AW, HW, VELW)
            Flux_t = Flux_t * AW
! Bottom face (i, j+1/2, k-1/2)
            Flux_b = (Theta_V_nt(IJKM) * Flux_gT(IJKM) + &
                    Theta_V_st(IJKM) * Flux_gT(IJPKM))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',&
               ALPHA_Vt_c(IJKM), AW, HW, VELW)
            Flux_b = Flux_b * AW
         ENDIF
      ELSE
         Flux_e = HALF * (Flux_gE(IJK) + Flux_gE(IJPK))
         Flux_w = HALF * (Flux_gE(IMJK) + Flux_gE(IMJPK))

         Flux_n = HALF * (Flux_gN(IJK) + Flux_gN(IJPK))
         Flux_s = HALF * (Flux_gN(IJMK) + Flux_gN(IJK))

         IF (DO_K) THEN
            Flux_t = HALF * (Flux_gT(IJK) + Flux_gT(IJPK))
            Flux_b = HALF * (Flux_gT(IJKM) + Flux_gT(IJPKM))
         ENDIF
      ENDIF   ! end if/else cut_v_treatmeant_at

      RETURN
      END SUBROUTINE GET_VCELL_GCFLUX_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  the faces of a v-momentum cell. Note the fluxes are calculated at   C
!  all faces regardless of flow_at_n of condition of the west, south   C
!  or bottom face.                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_VCELL_GDIFF_TERMS(D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, IJK)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_v_treatment_at
      USE cutcell, only: oneodx_e_v, oneody_n_v, oneodz_t_v

      USE functions, only: wall_at
      USE functions, only: east_of, north_of, top_of
      USE functions, only: west_of, bottom_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: odx_e, ody, odz_t
      USE geometry, only: do_k
      USE geometry, only: ox
      USE geometry, only: ayz_v, axz_v, axy_v

      USE indices, only: i_of, j_of, k_of
      USE indices, only: jp1, im1, km1

      use matrix, only: e, w, n, s, t, b
      USE param1, only: zero
      use fldvar
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! diffusion through faces of given ijk v-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: d_fe, d_fw
      DOUBLE PRECISION, INTENT(OUT) :: d_fn, d_fs
      DOUBLE PRECISION, INTENT(OUT) :: d_ft, d_fb
! ijk index
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: i, j, k, jp, im, km
      INTEGER :: ijkc, ijkn, ijke, ijkne, ijkw, ijknw
      INTEGER :: ijkt, ijktn, ijkb, ijkbn
! length terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
!---------------------------------------------------------------------//

      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = IM1(I)
      JP = JP1(J)
      KM = KM1(K)

      IJKN = NORTH_OF(IJK)
      IF (WALL_AT(IJK)) THEN
         IJKC = IJKN
      ELSE
         IJKC = IJK
      ENDIF
      IJKE = EAST_OF(IJK)
      IJKNE = EAST_OF(IJKN)
      IJKW = WEST_OF(IJK)
      IJKNW = NORTH_OF(IJKW)

      IF(CUT_V_TREATMENT_AT(IJK)) THEN
         C_AE = ONEoDX_E_V(IJK)
         C_AW = ONEoDX_E_V(IMJK)
         C_AN = ONEoDY_N_V(IJK)
         C_AS = ONEoDY_N_V(IJMK)
         C_AT = ONEoDZ_T_V(IJK)
         C_AB = ONEoDZ_T_V(IJKM)
      ELSE
         C_AE = ODX_E(I)
         C_AW = ODX_E(IM)
         C_AN = ODY(JP)
         C_AS = ODY(J)
         C_AT = ODZ_T(K)
         C_AB = ODZ_T(KM)
      ENDIF

! East face (i+1/2, j+1/2, k)
      D_Fe = AVG_Y_H(AVG_X_H(MU_G(IJKC),MU_G(IJKE),I),&
                     AVG_X_H(MU_G(IJKN),MU_G(IJKNE),I),J)*&
             C_AE*AYZ_V(IJK)
! West face (i-1/2, j+1/2, k)
      D_Fw = AVG_Y_H(AVG_X_H(MU_G(IJKW),MU_G(IJKC),IM),&
                     AVG_X_H(MU_G(IJKNW),MU_G(IJKN),IM),J)*&
             C_AW*AYZ_V(IMJK)

! North face (i, j+1, k)
      D_Fn = MU_G(IJKN)*C_AN*AXZ_V(IJK)
! South face (i, j, k)
      D_Fs = MU_G(IJKC)*C_AS*AXZ_V(IJMK)

      D_FT = ZERO
      D_FB = ZERO
      IF (DO_K) THEN
         IJKT = TOP_OF(IJK)
         IJKTN = NORTH_OF(IJKT)
         IJKB = BOTTOM_OF(IJK)
         IJKBN = NORTH_OF(IJKB)

! Top face (i, j+1/2, k+1/2)
         D_Ft = AVG_Y_H(AVG_Z_H(MU_G(IJKC),MU_G(IJKT),K),&
                        AVG_Z_H(MU_G(IJKN),MU_G(IJKTN),K),J)*&
                OX(I)*C_AT*AXY_V(IJK)
! Bottom face (i, j+1/2, k-1/2)
         D_Fb = AVG_Y_H(AVG_Z_H(MU_G(IJKB),MU_G(IJKC),KM),&
                        AVG_Z_H(MU_G(IJKBN),MU_G(IJKN),KM),J)*&
                OX(I)*C_AB*AXY_V(IJKM)
      ENDIF   ! end if (do_k)


      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE GET_VCELL_GDIFF_TERMS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_V_g0                                            C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs. C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative. See      C
!  source_v_g.                                                         C
!  Implement FOUP discretization                                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-JUN-96   C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_V_G0(A_V_G)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE functions, only: flow_at_n
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: do_k

      USE param, only: dimension_3
      USE param1, only: zero
      use matrix, only: e, w, n, s, t, b

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_V_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_V_g(DIMENSION_3, -3:3)

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

!$omp     parallel do default(none)                              &
!$omp     private(IJK, IMJK, IJMK, IPJK, IJPK, IJKM, IJKP,       &
!$omp             D_fe, d_fw, d_fn, d_fs, d_ft, d_fb,            &
!$omp             flux_e, flux_w, flux_n, flux_s, flux_t,        &
!$omp             flux_b)                                        &
!$omp     shared(ijkstart3, ijkend3, do_k, a_v_g)
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_N(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_VCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk)
            CALL GET_VCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1/2, j+1/2, k)
            IF (Flux_e >= ZERO) THEN
               A_V_G(IJK,E) = D_Fe
               A_V_G(IPJK,W) = D_Fe + Flux_e
            ELSE
               A_V_G(IJK,E) = D_Fe - Flux_e
               A_V_G(IPJK,W) = D_Fe
            ENDIF
! West face (i-1/2, j+1/2, k)
            IF (.NOT.FLOW_AT_N(IMJK)) THEN
               IF (Flux_w >= ZERO) THEN
                  A_V_G(IJK,W) = D_Fw + Flux_w
               ELSE
                  A_V_G(IJK,W) = D_Fw
               ENDIF
            ENDIF


! North face (i, j+1, k)
            IF (Flux_n >= ZERO) THEN
               A_V_G(IJK,N) = D_Fn
               A_V_G(IJPK,S) = D_Fn + Flux_n
            ELSE
               A_V_G(IJK,N) = D_Fn - Flux_n
               A_V_G(IJPK,S) = D_Fn
            ENDIF
! South face (i, j, k)
            IF (.NOT.FLOW_AT_N(IJMK)) THEN
               IF (Flux_s >= ZERO) THEN
                  A_V_G(IJK,S) = D_Fs + Flux_s
               ELSE
                  A_V_G(IJK,S) = D_Fs
               ENDIF
            ENDIF


            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)

! Top face (i, j+1/2, k+1/2)
               IF (Flux_t >= ZERO) THEN
                  A_V_G(IJK,T) = D_Ft
                  A_V_G(IJKP,B) = D_Ft + Flux_t
               ELSE
                  A_V_G(IJK,T) = D_Ft - Flux_t
                  A_V_G(IJKP,B) = D_Ft
               ENDIF
! Bottom face (i, j+1/2, k-1/2)
               IF (.NOT.FLOW_AT_N(IJKM)) THEN
                  IF (Flux_b >= ZERO) THEN
                     A_V_G(IJK,B) = D_Fb + Flux_b
                  ELSE
                     A_V_G(IJK,B) = D_Fb
                  ENDIF
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if flow_at_n
      END DO   ! end do ijk
!$omp end parallel do

      RETURN
      END SUBROUTINE STORE_A_V_G0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_V_g1                                            C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative.          C
!  Implements higher order discretization.                             C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_V_G1(A_V_G)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE fldvar, only: v_g

      USE functions, only: flow_at_n
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: do_k

      USE param, only: dimension_3
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
! Septadiagonal matrix A_V_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_V_g(DIMENSION_3, -3:3)

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

      CALL GET_VCELL_GVTERMS(U, V, WW)

! shear indicator: y-momentum
      incr=2
      CALL CALC_XSI (DISCRETIZE(4), V_G, U, V, WW, XSI_E, XSI_N, &
         XSI_T, incr)


!!!$omp      parallel do                                             &
!!!$omp&     private(IJK, IMJK, IJMK, IPJK, IJPK, IJKM, IJKP,        &
!!!$omp&             d_fe, d_fw, d_fn, d_fs, d_ft, d_fb,             &
!!!$omp&             flux_e, flux_w, flux_n, flux_s, flux_t, flux_b)
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_N(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_VCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk)
            CALL GET_VCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)

            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

! East face (i+1/2, j+1/2, k)
            A_V_G(IJK,E) = D_Fe - XSI_E(IJK)*Flux_e
            A_V_G(IPJK,W) = D_Fe + (ONE - XSI_E(IJK))*Flux_e
! West face (i-1/2, j+1/2, k)
            IF (.NOT.FLOW_AT_N(IMJK)) THEN
               A_V_G(IJK,W) = D_Fw + (ONE - XSI_E(IMJK))*Flux_w
            ENDIF


! North face (i, j+1, k)
            A_V_G(IJK,N) = D_Fn - XSI_N(IJK)*Flux_n
            A_V_G(IJPK,S) = D_Fn + (ONE - XSI_N(IJK))*Flux_n
! South face (i, j, k)
            IF (.NOT.FLOW_AT_N(IJMK)) THEN
               A_V_G(IJK,S) = D_Fs + (ONE - XSI_N(IJMK))*Flux_s
            ENDIF

            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)
! Top face (i, j+1/2, k+1/2)
               A_V_G(IJK,T) = D_Ft - XSI_T(IJK)*Flux_t
               A_V_G(IJKP,B) = D_Ft + (ONE - XSI_T(IJK))*Flux_t
! Bottom face (i, j+1/2, k-1/2)
               IF (.NOT.FLOW_AT_N(IJKM)) THEN
                  A_V_G(IJK,B) = D_Fb + (ONE - XSI_T(IJKM))*Flux_b
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if (flow_at_n)
      ENDDO   ! end do (ijk)

      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_V_G1
