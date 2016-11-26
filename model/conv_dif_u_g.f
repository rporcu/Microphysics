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
      USE param, only: dimension_3
      USE run, only: momentum_x_eq
      USE run, only: discretize
      use fldvar
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
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
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      USE fldvar, only: u_g, v_g, w_g

      USE fun_avg, only: avg
      USE functions, only: funijk, iplus
      USE geometry, only: do_k
      USE functions, only:  ip1

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
      INTEGER :: IJK, I,j,k, IP, IPJK
!---------------------------------------------------------------------//

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3

            IJK = funijk(i,j,k)
            IP = IP1(I)
            IPJK = FUNIJK(iplus(i,j,k),j,k)

            U(IJK) = AVG(U_G(IJK),U_G(IPJK))
            V(IJK) = AVG(V_G(IJK),V_G(IPJK))
            IF (DO_K) WW(IJK) = AVG(W_G(IJK),W_G(IPJK))
          ENDDO
        ENDDO
      ENDDO

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
         FLUX_S, FLUX_T, FLUX_B, I, J, K)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: funijk
      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus
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

      INTEGER, INTENT(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: ijk
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: ipjk, ipjmk, ipjkm

!---------------------------------------------------------------------//

! indices
      IJK  = funijk(I,J,K)

      IMJK = FUNIJK(iminus(i,j,k),j,k)
      IPJK = FUNIJK(iplus(i,j,k),j,k)
      IJMK = FUNIJK(i,jminus(i,j,k),k)
      IJKM = FUNIJK(i,j,kminus(i,j,k))
      IPJMK = FUNIJK(iplus(i,jminus(i,j,k),k),jminus(i,j,k),k)
      IPJKM = FUNIJK(iplus(i,j,kminus(i,j,k)),j,kminus(i,j,k))

! First calculate the fluxes at the faces
      Flux_e = HALF * (Flux_gE(IJK) + Flux_gE(IPJK))
      Flux_w = HALF * (Flux_gE(IMJK) + Flux_gE(IJK))
      Flux_n = HALF * (Flux_gN(IJK) + Flux_gN(IPJK))
      Flux_s = HALF * (Flux_gN(IJMK) + Flux_gN(IPJMK))

      IF (DO_K) THEN
         Flux_t = HALF * (Flux_gT(IJK) + Flux_gT(IPJK))
         Flux_b = HALF * (Flux_gT(IJKM) + Flux_gT(IPJKM))
      ENDIF
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
         D_FT, D_FB, I, J, K)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: funijk, wall_at
      USE functions, only: ieast, jnorth, ktop
      USE functions, only: ieast, jnorth, jsouth, ktop, kbot
      USE functions, only: iminus, jminus, kminus

      USE geometry, only: odx, ody, odz
      USE geometry, only: do_k
      USE geometry, only: ayz, axz, axy

      USE functions, only: ip1, jm1, km1

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
      INTEGER, INTENT(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: ijk
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: ip, jm, km
      INTEGER :: ijkc, ijke, ijkn, ijkne, ijks, ijkse
      INTEGER :: ijkt, ijkte, ijkb, ijkbe
      INTEGER :: jtmp,ktmp
! length terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
!---------------------------------------------------------------------//

      IJK = funijk(i,j,k)

      IMJK = FUNIJK(iminus(i,j,k),j,k)
      IJMK = FUNIJK(i,jminus(i,j,k),k)
      IJKM = FUNIJK(i,j,kminus(i,j,k))

      IP = IP1(I)
      JM = JM1(J)
      KM = KM1(K)

      IJKE = FUNIJK(ieast(i,j,k),j,k)

      IF (wall_at(i,j,k))  THEN
         IJKC = IJKE
      ELSE
         IJKC = IJK
      ENDIF

      jtmp = jnorth(i,j,k)
      IJKN = FUNIJK(i,jtmp,k)
      IJKNE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

      jtmp = jsouth(i,j,k)
      IJKS = FUNIJK(i,jtmp,k)
      IJKSE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

      C_AE = ODX
      C_AW = ODX
      C_AN = ODY
      C_AS = ODY
      C_AT = ODZ
      C_AB = ODZ

! East face (i+1, j, k)
      D_FE = MU_G(IJKE)*C_AE*AYZ
! West face (i, j, k)
      D_FW = MU_G(IJKC)*C_AW*AYZ


! North face (i+1/2, j+1/2, k)
      D_FN = AVG_H(AVG_H(MU_G(IJKC),MU_G(IJKN)),&
                   AVG_H(MU_G(IJKE),MU_G(IJKNE)))*C_AN*AXZ
! South face (i+1/2, j-1/2, k)
      D_FS = AVG_H(AVG_H(MU_G(IJKS),MU_G(IJKC)),&
                   AVG_H(MU_G(IJKSE),MU_G(IJKE)))*C_AS*AXZ

      D_FT = ZERO
      D_FB = ZERO
      IF (DO_K) THEN

         ktmp = ktop(i,j,k)
         IJKT  = FUNIJK(i,j,ktmp)
         IJKTE = FUNIJK(ieast(i,j,ktmp),j,ktmp)

         ktmp  = kbot(i,j,k)
         IJKB  = FUNIJK(i,j,ktmp)
         IJKBE = FUNIJK(ieast(i,j,ktmp),j,ktmp)

! Top face (i+1/2, j, k+1/2)
         D_FT = AVG_H(AVG_H(MU_G(IJKC),MU_G(IJKT)),&
                      AVG_H(MU_G(IJKE),MU_G(IJKTE)))*C_AT*AXY
! Bottom face (i+1/2, j, k-1/2)
         D_FB = AVG_H(AVG_H(MU_G(IJKB),MU_G(IJKC)),&
                      AVG_H(MU_G(IJKBE),MU_G(IJKE)))*C_AB*AXY
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
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      USE functions, only: funijk
      USE functions, only: flow_at_e
      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus

      USE geometry, only: do_k

      USE param, only: dimension_3
      USE param1, only: zero
      use matrix, only: e, w, n, s, t, b
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_U_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)

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

         IF (FLOW_AT_E(i,j,k)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_UCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, i, j, k)
            CALL GET_UCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, i, j, k)

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)

! East face (i+1, j, k)
            IF (Flux_e >= ZERO) THEN
               A_U_G(I,J,K,E) = D_Fe
               A_U_G(iplus(i,j,k),j,k,W) = D_Fe + Flux_e
            ELSE
               A_U_G(I,J,K,E) = D_Fe - Flux_e
               A_U_G(iplus(i,j,k),j,k,W) = D_Fe
            ENDIF
! West face (i, j, k)
            IF (.NOT.FLOW_AT_E(iminus(i,j,k),j,k)) THEN
               IF (Flux_w >= ZERO) THEN
                  A_U_G(I,J,K,W) = D_Fw + Flux_w
               ELSE
                  A_U_G(I,J,K,W) = D_Fw
               ENDIF
            ENDIF


! North face (i+1/2, j+1/2, k)
            IF (Flux_n >= ZERO) THEN
               A_U_G(I,J,K,N) = D_Fn
               A_U_G(i,jplus(i,j,k),k,S) = D_Fn + Flux_n
            ELSE
               A_U_G(I,J,K,N) = D_Fn - Flux_n
               A_U_G(i,jplus(i,j,k),k,S) = D_Fn
            ENDIF
! South face (i+1/2, j-1/2, k)
            IF (.NOT.FLOW_AT_E(i,jminus(i,j,k),k)) THEN
               IF (Flux_s >= ZERO) THEN
                  A_U_G(I,J,K,S) = D_Fs + Flux_s
               ELSE
                  A_U_G(I,J,K,S) = D_Fs
               ENDIF
            ENDIF


            IF (DO_K) THEN

               IJKM = FUNIJK(i,j,kminus(i,j,k))
               IJKP = FUNIJK(i,j,kplus(i,j,k))

! Top face (i+1/2, j, k+1/2)
               IF (Flux_t >= ZERO) THEN
                  A_U_G(I,J,K,T) = D_Ft
                  A_U_G(i,j,kplus(i,j,k),B) = D_Ft + Flux_t
               ELSE
                  A_U_G(I,J,K,T) = D_Ft - Flux_t
                  A_U_G(i,j,kplus(i,j,k),B) = D_Ft
               ENDIF
! Bottom face (i+1/2, j, k-1/2)
               IF (.NOT.FLOW_AT_E(i,j,kminus(i,j,k))) THEN
                  IF (Flux_b >= ZERO) THEN
                     A_U_G(I,J,K,B) = D_Fb + Flux_b
                  ELSE
                     A_U_G(I,J,K,B) = D_Fb
                  ENDIF
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if (flow_at_e)
      ENDDO
      ENDDO
      ENDDO

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
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3
      USE fldvar, only: u_g

      USE functions, only: funijk
      USE functions, only: flow_at_e
      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus

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
! Septadiagonal matrix A_U_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I,J,K,IJK, IPJK, IMJK, IJPK, IJMK, IJKP, IJKM
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

      CALL GET_UCELL_GVTERMS(U, V, WW)

! shear indicator:
      incr=1
      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, &
                     XSI_T, incr)

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = funijk(i,j,k)

         IF (FLOW_AT_E(i,j,k)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_UCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, i, j, k)
            CALL GET_UCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, i, j, k)

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)

! East face (i+1, j, k)
            A_U_G(I,J,K,E) = D_Fe - XSI_E(IJK) * Flux_e
            A_U_G(iplus(i,j,k),j,k,W) = D_Fe + (ONE - XSI_E(IJK)) * Flux_e
! West face (i, j, k)
            IF (.NOT.FLOW_AT_E(iminus(i,j,k),j,k)) THEN
               A_U_G(I,J,K,W) = D_Fw + (ONE - XSI_E(IMJK)) * Flux_w
            ENDIF


! North face (i+1/2, j+1/2, k)
            A_U_G(I,J,K,N) = D_Fn - XSI_N(IJK) * Flux_n
            A_U_G(i,jplus(i,j,k),k,S) = D_Fn + (ONE - XSI_N(IJK)) * Flux_n
! South face (i+1/2, j-1/2, k)
            IF (.NOT.FLOW_AT_E(i,jminus(i,j,k),k)) THEN
               A_U_G(I,J,K,S) = D_Fs + (ONE - XSI_N(IJMK)) * Flux_s
            ENDIF

! Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN

               IJKM = FUNIJK(i,j,kminus(i,j,k))
               IJKP = FUNIJK(i,j,kplus(i,j,k))

               A_U_G(I,J,K,T) = D_Ft - XSI_T(IJK) * Flux_t
               A_U_G(i,j,kplus(i,j,k),B) = D_Ft + (ONE - XSI_T(IJK)) * Flux_t
! Bottom face (i+1/2, j, k-1/2)
               IF (.NOT.FLOW_AT_E(i,j,kminus(i,j,k))) THEN
                  A_U_G(I,J,K,B) = D_Fb + (ONE - XSI_T(IJKM)) * Flux_b
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if flow_at_e
      ENDDO
      ENDDO
      ENDDO

      deallocate( U, V, WW )

      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_U_G1
