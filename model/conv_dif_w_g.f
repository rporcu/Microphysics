!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_W_g                                            C
!  Purpose: Determine convection diffusion terms for w_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_DIF_W_G(A_M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3
      USE run, only: momentum_z_eq
      USE run, only: discretize
      use fldvar
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
!---------------------------------------------------------------------//

      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN

      IF (DISCRETIZE(5) == 0) THEN               ! 0 & 1 => FOUP
         CALL STORE_A_W_G0 (A_M)
      ELSE
         CALL STORE_A_W_G1 (A_M)
      ENDIF

      RETURN
      END SUBROUTINE CONV_DIF_W_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a w-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_WCELL_GVTERMS(U, V, WW)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      USE fldvar, only: u_g, v_g, w_g

      USE fun_avg, only: avg_z_t, avg_z
      USE functions, only: funijk
      USE functions, only: kp_of
      USE indices, only: k_of

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
      INTEGER :: IJK, I, J, K, IJKP
! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//


      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = funijk(i,j,k)
         IJKP = KP_OF(IJK)

         U(IJK) = AVG_Z(U_G(IJK),U_G(IJKP),K)
         V(IJK) = AVG_Z(V_G(IJK),V_G(IJKP),K)
         WW(IJK) = AVG_Z_T(W_G(IJK),W_G(IJKP),0)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GET_WCELL_GVTERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the convective fluxes through the faces of a     C
!  w-momentum cell. Note the fluxes are calculated at all faces of     C
!  regardless of flow_at_t of condition of the west, south, or         C
!  bottom face.                                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_WCELL_GCFLUX_TERMS(FLUX_E, FLUX_W, FLUX_N, &
         FLUX_S, FLUX_T, FLUX_B, I, J, K)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: funijk, kp_of, im_of, jm_of, km_of

      USE fldvar, only: flux_ge, flux_gn, flux_gt

      USE param1, only: half
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! fluxes through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: flux_e, flux_w
      DOUBLE PRECISION, INTENT(OUT) :: flux_n, flux_s
      DOUBLE PRECISION, INTENT(OUT) :: flux_t, flux_b
! ijk index
      INTEGER, INTENT(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: ijk, imjk, ijmk, ijkm
      INTEGER :: ijkp, imjkp, ijmkp

! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//

      IJK = funijk(i,j,k)

      IJKP = KP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IMJKP = KP_OF(IMJK)
      IJMKP = KP_OF(IJMK)

      Flux_e = HALF * (Flux_gE(IJK) + Flux_gE(IJKP))
      Flux_w = HALF * (Flux_gE(IMJK) + Flux_gE(IMJKP))
      Flux_n = HALF * (Flux_gN(IJK) + Flux_gN(IJKP))
      Flux_s = HALF * (Flux_gN(IJMK) + Flux_gN(IJMKP))
      Flux_t = HALF * (Flux_gT(IJK) + Flux_gT(IJKP))
      Flux_b = HALF * (Flux_gT(IJKM) + Flux_gT(IJK))

      RETURN
      END SUBROUTINE GET_WCELL_GCFLUX_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a w-momentum cell. Note the fluxes are calculated at       C
!  all faces regardless of flow_at_t condition of the west, south      C
!  or bottom face.                                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_WCELL_GDIFF_TERMS(D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, I, J, K)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: funijk, wall_cell
      USE functions, only: ieast, iwest, jnorth, jsouth, ktop
      USE functions, only: iminus, jminus, kminus

      USE geometry, only: odx_e, ody_n, odz
      USE geometry, only: ox
      USE geometry, only: ayz_w, axz_w, axy_w

      USE indices, only: i_of, j_of, k_of
      USE indices, only: kp1, im1, jm1

      use matrix, only: e, w, n, s, t, b
      use fldvar
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! diffusion through faces of given ijk w-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: d_fe, d_fw
      DOUBLE PRECISION, INTENT(OUT) :: d_fn, d_fs
      DOUBLE PRECISION, INTENT(OUT) :: d_ft, d_fb

      INTEGER, INTENT(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: ijk, imjk, ijmk, ijkm
      INTEGER :: kp, im, jm
      INTEGER :: ijkc, ijkt, ijke, ijkte, ijkw, ijkwt
      INTEGER :: ijkn, ijktn, ijks, ijkst
      INTEGER :: itmp, jtmp, ktmp
! length terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
!---------------------------------------------------------------------//

      IJK = funijk(i,j,k)

      IMJK = FUNIJK(iminus(i,j,k),j,k)
      IJMK = FUNIJK(i,jminus(i,j,k),j)
      IJKM = FUNIJK(i,j,kminus(i,j,k))

      KP = KP1(K)
      IM = IM1(I)
      JM = JM1(J)

      ktmp = ktop(i,j,k)
      IJKT = FUNIJK(i,j,ktmp)

      IJKE  = FUNIJK(ieast(i,j,k),j,k)

      itmp = iwest(i,j,k)
      IJKW  = FUNIJK(itmp,j,k)
      IJKWT = FUNIJK(itmp,j,ktop(itmp,j,k))

      IJKTE = FUNIJK(ieast(i,j,ktmp),j,ktmp)
      IJKTN = FUNIJK(i,jnorth(i,j,ktmp),ktmp)

      IJKN = FUNIJK(i,jnorth(i,j,k),k)

      jtmp  = jsouth(i,j,k)
      IJKS  = FUNIJK(i,jtmp,k)
      IJKST = FUNIJK(i,jtmp,ktop(i,jtmp,k))

      IF (wall_cell(i,j,k)) THEN
         IJKC = IJKT
      ELSE
         IJKC = IJK
      ENDIF

      C_AE = ODX_E(I)
      C_AW = ODX_E(IM)
      C_AN = ODY_N(J)
      C_AS = ODY_N(JM)
      C_AT = ODZ(KP)
      C_AB = ODZ(K)

! East face (i+1/2, j, k+1/2)
      D_Fe = AVG_Z_H(AVG_X_H(MU_G(IJKC),MU_G(IJKE),I),&
                     AVG_X_H(MU_G(IJKT),MU_G(IJKTE),I),K)*&
               C_AE*AYZ_W(IJK)
! West face (i-1/2, j, k+1/2)
      D_Fw = AVG_Z_H(AVG_X_H(MU_G(IJKW),MU_G(IJKC),IM),&
                     AVG_X_H(MU_G(IJKWT),MU_G(IJKT),IM),K)*&
                C_AW*AYZ_W(IMJK)

! North face (i, j+1/2, k+1/2)
      D_Fn = AVG_Z_H(AVG_Y_H(MU_G(IJKC),MU_G(IJKN),J),&
                     AVG_Y_H(MU_G(IJKT),MU_G(IJKTN),J),K)*&
                C_AN*AXZ_W(IJK)
! South face (i, j-1/2, k+1/2)
      D_Fs = AVG_Z_H(AVG_Y_H(MU_G(IJKS),MU_G(IJKC),JM),&
                     AVG_Y_H(MU_G(IJKST),MU_G(IJKT),JM),K)*&
                C_AS*AXZ_W(IJMK)

! Top face (i, j, k+1)
      D_Ft = MU_G(IJKT)*OX(I)*C_AT*AXY_W(IJK)
! Bottom face (i, j, k)
      D_Fb = MU_G(IJK)*OX(I)*C_AB*AXY_W(IJKM)

      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE GET_WCELL_GDIFF_TERMS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_W_g0                                            C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs. C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative. See      C
!  source_w_g.                                                         C
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
      SUBROUTINE STORE_A_W_G0(A_W_G)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      USE functions, only: funijk
      USE functions, only: flow_at_t
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE param, only: dimension_3
      USE param1, only: zero
      use matrix, only: e, w, n, s, t, b

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_U_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_W_g(DIMENSION_3, -3:3)
! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Face mass flux
      DOUBLE PRECISION :: flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b
! Diffusion parameter
      DOUBLE PRECISION :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

!---------------------------------------------------------------------//

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = funijk(i,j,k)

         IF (FLOW_AT_T(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_WCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, i, j, k)

            CALL GET_WCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, i, j, k)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

! East face (i+1/2, j, k+1/2)
            IF (Flux_e >= ZERO) THEN
               A_W_G(IJK,E) = D_Fe
               A_W_G(IPJK,W) = D_Fe + Flux_e
            ELSE
               A_W_G(IJK,E) = D_Fe - Flux_e
               A_W_G(IPJK,W) = D_Fe
            ENDIF
! West face (i-1/2, j, k+1/2)
            IF (.NOT.FLOW_AT_T(IMJK)) THEN
               IF (Flux_w >= ZERO) THEN
                  A_W_G(IJK,W) = D_Fw + Flux_w
               ELSE
                  A_W_G(IJK,W) = D_Fw
               ENDIF
            ENDIF


! North face (i, j+1/2, k+1/2)
            IF (Flux_n >= ZERO) THEN
               A_W_G(IJK,N) = D_Fn
               A_W_G(IJPK,S) = D_Fn + Flux_n
            ELSE
               A_W_G(IJK,N) = D_Fn - Flux_n
               A_W_G(IJPK,S) = D_Fn
            ENDIF
! South face (i, j-1/2, k+1/2)
            IF (.NOT.FLOW_AT_T(IJMK)) THEN
              IF (Flux_s >= ZERO) THEN
                  A_W_G(IJK,S) = D_Fs + Flux_s
               ELSE
                  A_W_G(IJK,S) = D_Fs
               ENDIF
            ENDIF


! Top face (i, j, k+1)
            IF (Flux_T >= ZERO) THEN
               A_W_G(IJK,T) = D_Ft
               A_W_G(IJKP,B) = D_Ft + Flux_t
            ELSE
               A_W_G(IJK,T) = D_Ft - Flux_t
               A_W_G(IJKP,B) = D_Ft
            ENDIF
! Bottom face (i, j, k)
            IF (.NOT.FLOW_AT_T(IJKM)) THEN
               IF (Flux_b >= ZERO) THEN
                  A_W_G(IJK,B) = D_Fb + Flux_b
               ELSE
                  A_W_G(IJK,B) = D_Fb
               ENDIF
            ENDIF
         ENDIF   ! end if (flow_at_t)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE STORE_A_W_G0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_W_g1                                            C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative.          C
!  Implements higher order discretization.                             C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_W_G1(A_W_G)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3
      USE fldvar, only: w_g

      USE functions, only: funijk
      USE functions, only: flow_at_t
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE param, only: dimension_3
      USE param1, only: one

      use matrix, only: e, w, n, s, t, b

      USE run, only: discretize

      USE xsi, only: calc_xsi
      USE xsi_array, only: xsi_e, xsi_n, xsi_t
      USE xsi_array, only: lock_xsi_array, unlock_xsi_array

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_W_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_W_g(DIMENSION_3, -3:3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I,J,K, IJK, IPJK, IMJK, IJPK, IJMK, IJKP, IJKM
! indicator for shear
      INTEGER :: incr
! Diffusion parameter
      DOUBLE PRECISION :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      DOUBLE PRECISION :: Flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b

! temporary use of global arrays:
! x directional velocity
      DOUBLE PRECISION, allocatable :: U(:)
! y directional velocity
      DOUBLE PRECISION, allocatable :: V(:)
! z directional velocity
      DOUBLE PRECISION, allocatable :: WW(:)
!---------------------------------------------------------------------//

      call lock_xsi_array

      allocate(  U(DIMENSION_3) )
      allocate(  V(DIMENSION_3) )
      allocate( WW(DIMENSION_3) )

      CALL GET_WCELL_GVTERMS(U, V, WW)

! shear indicator:
      incr=0
      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N,&
         XSI_T, incr)

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = funijk(i,j,k)

         IF (FLOW_AT_T(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_WCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, i, j, k)
            CALL GET_WCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, i, j, k)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

! East face (i+1/2, j, k+1/2)
            A_W_G(IJK,E) = D_Fe - XSI_E(IJK)*Flux_e
            A_W_G(IPJK,W) = D_Fe + (ONE - XSI_E(IJK))*Flux_e
! West face (i-1/2, j, k+1/2)
            IF (.NOT.FLOW_AT_T(IMJK)) THEN
               A_W_G(IJK,W) = D_Fw + (ONE - XSI_E(IMJK))*Flux_w
            ENDIF


! North face (i, j+1/2, k+1/2)
            A_W_G(IJK,N) = D_Fn - XSI_N(IJK)*Flux_n
            A_W_G(IJPK,S) = D_Fn + (ONE - XSI_N(IJK))*Flux_n
! South face (i, j-1/2, k+1/2)
            IF (.NOT.FLOW_AT_T(IJMK)) THEN
               A_W_G(IJK,S) = D_Fs + (ONE - XSI_N(IJMK))*Flux_s
            ENDIF

! Top face (i, j, k+1)
            A_W_G(IJK,T) = D_Ft - XSI_T(IJK)*Flux_t
            A_W_G(IJKP,B) = D_Ft + (ONE - XSI_T(IJK))*Flux_t
! Bottom face (i, j, k)
            IF (.NOT.FLOW_AT_T(IJKM)) THEN
              A_W_G(IJK,B) = D_Fb + (ONE - XSI_T(IJKM))*Flux_b
            ENDIF

         ENDIF   ! end if flow_at_t
      ENDDO
      ENDDO
      ENDDO

      deallocate( U, V, WW )
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_W_G1
