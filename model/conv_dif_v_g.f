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
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      USE fldvar, only: u_g, v_g, w_g

      USE functions, only: funijk, iplus, iminus, jplus, jminus
      USE fun_avg, only: avg_y_n, avg_y
      USE geometry, only: do_k

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
      INTEGER :: IJK, I, J, K, IJPK
!---------------------------------------------------------------------//

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = funijk(i,j,k)
         IJPK = FUNIJK(i,jplus(i,j,k),k)

         U(IJK) = AVG_Y(U_G(IJK),U_G(IJPK),J)
         V(IJK) = AVG_Y_N(V_G(IJK),V_G(IJPK),0)
         IF (DO_K) WW(IJK) = AVG_Y(W_G(IJK),W_G(IJPK),J)
      ENDDO
      ENDDO
      ENDDO

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
         FLUX_S, FLUX_T, FLUX_B, I, J, K)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: funijk
      USE functions, only: iplus, iminus, jplus, jminus, kplus, kminus

      USE geometry, only: do_k

      USE fldvar, only: flux_ge, flux_gn, flux_gt

      USE param1, only: half
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! fluxes through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: flux_e, flux_w
      DOUBLE PRECISION, INTENT(OUT) :: flux_n, flux_s
      DOUBLE PRECISION, INTENT(OUT) :: flux_t, flux_b
! indices
      INTEGER, INTENT(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: ijk, imjk, ijmk, ijkm
      INTEGER :: ijpk, imjpk, ijpkm
      INTEGER :: itmp,ktmp

!---------------------------------------------------------------------//
! indices
      IJK  = funijk(i,j,k)

      IJPK = FUNIJK(i,jplus(i,j,k),k)
      IJMK = FUNIJK(i,jminus(i,j,k),k)

      itmp  = iminus(i,j,k)
      IMJK  = FUNIJK(itmp,j,k)
      IMJPK = FUNIJK(itmp,jplus(itmp,j,k),k)

      ktmp  = kminus(i,j,k)
      IJKM  = FUNIJK(i,j,ktmp)
      IJPKM = FUNIJK(i,jplus(i,j,ktmp),ktmp)


      Flux_e = HALF * (Flux_gE(IJK) + Flux_gE(IJPK))
      Flux_w = HALF * (Flux_gE(IMJK) + Flux_gE(IMJPK))

      Flux_n = HALF * (Flux_gN(IJK) + Flux_gN(IJPK))
      Flux_s = HALF * (Flux_gN(IJMK) + Flux_gN(IJK))

      IF (DO_K) THEN
         Flux_t = HALF * (Flux_gT(IJK) + Flux_gT(IJPK))
         Flux_b = HALF * (Flux_gT(IJKM) + Flux_gT(IJPKM))
      ENDIF

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
         D_FT, D_FB, I, J, K)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: funijk, wall_cell
      USE functions, only: iminus, jminus, kminus
      USE functions, only: ieast, iwest, jnorth, jsouth, kbot, ktop

      USE geometry, only: odx_e, ody, odz_t
      USE geometry, only: do_k
      USE geometry, only: ayz, axz, axy

      USE functions, only: im1, jp1, km1

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
      INTEGER, INTENT(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: ijk, imjk, ijmk, ijkm
      INTEGER :: jp, im, km
      INTEGER :: ijkc, ijkn, ijke, ijkne, ijkw, ijknw
      INTEGER :: ijkt, ijktn, ijkb, ijkbn
      INTEGER :: itmp, jtmp, ktmp
! length terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
!---------------------------------------------------------------------//

      IJK = funijk(i,j,k)

      IMJK = funijk(iminus(i,j,k),j,k)
      IJMK = funijk(i,jminus(i,j,k),k)
      IJKM = funijk(i,j,kminus(i,j,k))

      IM = IM1(I)
      JP = JP1(J)
      KM = KM1(K)

      jtmp  = jnorth(i,j,k)
      IJKN  = funijk(i,jtmp,k)
      IJKNE = funijk(ieast(i,jtmp,k),jtmp,k)

      IF (wall_cell(i,j,k)) THEN
         IJKC = IJKN
      ELSE
         IJKC = IJK
      ENDIF

      IJKE = funijk(ieast(i,j,k),j,k)

      itmp  = iwest(i,j,k)
      IJKW  = funijk(itmp,j,k)
      IJKNW = funijk(itmp,jnorth(itmp,j,k),k)

      C_AE = ODX_E(I)
      C_AW = ODX_E(IM)
      C_AN = ODY(JP)
      C_AS = ODY(J)
      C_AT = ODZ_T(K)
      C_AB = ODZ_T(KM)

! East face (i+1/2, j+1/2, k)
      D_Fe = AVG_Y_H(AVG_X_H(MU_G(IJKC),MU_G(IJKE),I),&
                     AVG_X_H(MU_G(IJKN),MU_G(IJKNE),I),J)*&
             C_AE*AYZ
! West face (i-1/2, j+1/2, k)
      D_Fw = AVG_Y_H(AVG_X_H(MU_G(IJKW),MU_G(IJKC),IM),&
                     AVG_X_H(MU_G(IJKNW),MU_G(IJKN),IM),J)*&
             C_AW*AYZ

! North face (i, j+1, k)
      D_Fn = MU_G(IJKN)*C_AN*AXZ
! South face (i, j, k)
      D_Fs = MU_G(IJKC)*C_AS*AXZ

      D_FT = ZERO
      D_FB = ZERO
      IF (DO_K) THEN

         ktmp  = ktop(i,j,k)
         IJKT  = funijk(i,j,ktmp)
         IJKTN = funijk(i,jnorth(i,j,ktmp),ktmp)

         ktmp  = kbot(i,j,k)
         IJKB  = funijk(i,j,ktmp)
         IJKBN = funijk(i,jnorth(i,j,ktmp),ktmp)

! Top face (i, j+1/2, k+1/2)
         D_Ft = AVG_Y_H(AVG_Z_H(MU_G(IJKC),MU_G(IJKT),K),&
                        AVG_Z_H(MU_G(IJKN),MU_G(IJKTN),K),J)*&
                C_AT*AXY
! Bottom face (i, j+1/2, k-1/2)
         D_Fb = AVG_Y_H(AVG_Z_H(MU_G(IJKB),MU_G(IJKC),KM),&
                        AVG_Z_H(MU_G(IJKBN),MU_G(IJKN),KM),J)*&
                C_AB*AXY
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
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      USE functions, only: funijk, iminus, iplus, jminus, jplus, kminus, kplus
      USE functions, only: flow_at_n

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

         IF (FLOW_AT_N(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_VCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, i, j, k)
            CALL GET_VCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, i, j, k)

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)

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
               IJKP = FUNIJK(i,j,kplus(i,j,k))
               IJKM = FUNIJK(i,j,kminus(i,j,k))

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
      ENDDO
      ENDDO
      ENDDO

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
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3
      USE fldvar, only: v_g

      USE functions, only: funijk, iminus, iplus, jminus, jplus, kminus, kplus
      USE functions, only: flow_at_n

      USE geometry, only: do_k

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
! Septadiagonal matrix A_V_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_V_g(DIMENSION_3, -3:3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IPJK, IMJK, IJPK, IJMK, IJKP, IJKM
      INTEGER :: I, J, K
! indicator for shear
      INTEGER :: incr
! Diffusion parameter
      DOUBLE PRECISION :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      DOUBLE PRECISION :: Flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b

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
      CALL GET_VCELL_GVTERMS(U, V, WW)

! shear indicator: y-momentum
      incr=2
      CALL CALC_XSI (DISCRETIZE(4), V_G, U, V, WW, XSI_E, XSI_N, &
         XSI_T, incr)


      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = funijk(i,j,k)

         IF (FLOW_AT_N(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_VCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, i, j, k)
            CALL GET_VCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, i, j, k)

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)

            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)

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
               IJKM = FUNIJK(i,j,kminus(i,j,k))
               IJKP = FUNIJK(i,j,kplus(i,j,k))
! Top face (i, j+1/2, k+1/2)
               A_V_G(IJK,T) = D_Ft - XSI_T(IJK)*Flux_t
               A_V_G(IJKP,B) = D_Ft + (ONE - XSI_T(IJK))*Flux_t
! Bottom face (i, j+1/2, k-1/2)
               IF (.NOT.FLOW_AT_N(IJKM)) THEN
                  A_V_G(IJK,B) = D_Fb + (ONE - XSI_T(IJKM))*Flux_b
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if (flow_at_n)
      ENDDO
      ENDDO
      ENDDO

      deallocate( U, V, WW )
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_V_G1
