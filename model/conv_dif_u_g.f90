module u_g_conv_dif

   use bl_fort_module, only : c_real

   private
   public :: conv_dif_u_g

   contains

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
      SUBROUTINE CONV_DIF_U_G(A_M, mu_g, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use run, only: discretize
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3

      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)

      real(c_real), intent(in   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      real(c_real), intent(in   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      real(c_real), intent(in   ) :: flux_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: flux_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: flux_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer    , intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
      real(c_real), intent(in   ) :: dt, dx, dy, dz
!---------------------------------------------------------------------//

      IF (DISCRETIZE(3) == 0) THEN
         CALL STORE_A_U_G0(A_M,MU_G,flux_ge,flux_gn,flux_gt, &
                           flag, dx, dy, dz)
      ELSE
         CALL STORE_A_U_G1(A_M,MU_G,u_g,v_g,w_g,flux_ge,flux_gn,flux_gt, &
                           flag, dt, dx, dy, dz)
      ENDIF

      END SUBROUTINE CONV_DIF_U_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a u-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GVTERMS(U, V, WW, u_g, v_g, w_g)

! Modules
!---------------------------------------------------------------------//
      use compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      use functions, only: avg
      use functions, only: iplus
      use functions, only:  ip1

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      real(c_real), intent(OUT) :: U&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(OUT) :: V&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(OUT) :: WW&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I,J,K
!---------------------------------------------------------------------//

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
            U(I,J,K) = AVG(U_G(I,J,K), U_G(iplus(i,j,k),j,k))
            V(I,J,K) = AVG(V_G(I,J,K), V_G(iplus(i,j,k),j,k))
            WW(I,J,K) = AVG(W_G(I,J,K), W_G(iplus(i,j,k),j,k))
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GET_UCELL_GVTERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the convective fluxes through the faces of a     C
!  u-momentum cell. Note the fluxes are calculated at all faces.       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GCFLUX_TERMS(&
         FLUX_E, FLUX_W, FLUX_N, &
         FLUX_S, FLUX_T, FLUX_B, &
         flux_ge, flux_gn, flux_gt, i, j, k)

! Modules
!---------------------------------------------------------------------//
      use compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use functions, only: iminus, iplus, jminus, jplus, kminus, kplus
      use param1   , only: half

      implicit none

!---------------------------------------------------------------------//
! Fluxes through faces of given u-momentum cell

      real(c_real), intent(OUT) :: flux_e, flux_w
      real(c_real), intent(OUT) :: flux_n, flux_s
      real(c_real), intent(OUT) :: flux_t, flux_b

      real(c_real), intent(in   ) :: flux_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: flux_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: flux_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      INTEGER         , intent( IN) :: i, j, k

!---------------------------------------------------------------------//

! First calculate the fluxes at the faces
      Flux_e = HALF * (Flux_gE(i,j,k) + Flux_gE(iplus(i,j,k),j,k))
      Flux_w = HALF * (Flux_gE(iminus(i,j,k),j,k) + Flux_gE(i,j,k))
      Flux_n = HALF * (Flux_gN(i,j,k) + Flux_gN(iplus(i,j,k),j,k))
      Flux_s = HALF * (Flux_gN(i,jminus(i,j,k),k) + &
         Flux_gN(iplus(i,jminus(i,j,k),k),jminus(i,j,k),k))

      Flux_t = HALF * (Flux_gT(i,j,k) + Flux_gT(iplus(i,j,k),j,k))
      Flux_b = HALF * (Flux_gT(i,j,kminus(i,j,k)) + &
         Flux_gT(iplus(i,j,kminus(i,j,k)),j,kminus(i,j,k)))

      END SUBROUTINE GET_UCELL_GCFLUX_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a u-momentum cell. Note the fluxes are calculated at       C
!  all faces.                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GDIFF_TERMS(&
         D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, MU_G, I, J, K, flag, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! diffusion through faces of given ijk u-momentum cell
      real(c_real), intent(OUT) :: d_fe, d_fw
      real(c_real), intent(OUT) :: d_fn, d_fs
      real(c_real), intent(OUT) :: d_ft, d_fb

      real(c_real), intent( IN) :: MU_G&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, intent( IN) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

      real(c_real), intent(in) :: dx, dy, dz

! ijk index
      INTEGER, intent(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: ip, jm, km, ic
! length terms
      real(c_real) :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
      real(c_real) :: odx, ody, odz
      real(c_real) :: axy, axz, ayz
!---------------------------------------------------------------------//

      odx = 1.0 / dx
      ody = 1.0 / dy
      odz = 1.0 / dz

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      IP = IP1(I)
      JM = JM1(J)
      KM = KM1(K)

      IF (flag(i,j,k,1)>=100)  THEN
         IC = ieast(i,j,k)
      ELSE
         IC = i
      ENDIF

      C_AE = ODX
      C_AW = ODX
      C_AN = ODY
      C_AS = ODY
      C_AT = ODZ
      C_AB = ODZ

! East face (i+1, j, k)
      D_FE = MU_G(ieast(i,j,k),j,k)*C_AE*AYZ
! West face (i, j, k)
      D_FW = MU_G(ic,j,k)*C_AW*AYZ

! North face (i+1/2, j+1/2, k)
      D_FN = AVG_H(AVG_H(MU_G(IC,J,K),MU_G(i,jnorth(i,j,k),k)),&
                   AVG_H(MU_G(ieast(i,j,k),j,k),MU_G(ieast(i,jnorth(i,j,k),k),jnorth(i,j,k),k)))*C_AN*AXZ
! South face (i+1/2, j-1/2, k)
      D_FS = AVG_H(AVG_H(MU_G(i,jsouth(i,j,k),k),MU_G(IC,J,K)),&
                   AVG_H(MU_G(ieast(i,jsouth(i,j,k),k),jsouth(i,j,k),k),MU_G(ieast(i,j,k),j,k)))*C_AS*AXZ

! Top face (i+1/2, j, k+1/2)
      D_FT = AVG_H(AVG_H(MU_G(IC,J,K),MU_G(i,j,ktop(i,j,k))),&
         AVG_H(MU_G(ieast(i,j,k),j,k),MU_G(ieast(i,j,ktop(i,j,k)),j,ktop(i,j,k))))*C_AT*AXY
! Bottom face (i+1/2, j, k-1/2)
         D_FB = AVG_H(AVG_H(MU_G(i,j,kbot(i,j,k)),MU_G(IC,J,K)),&
                      AVG_H(MU_G(ieast(i,j,kbot(i,j,k)),j,kbot(i,j,k)),MU_G(ieast(i,j,k),j,k)))*C_AB*AXY

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

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
      SUBROUTINE STORE_A_U_G0(A_U_G,MU_G,flux_ge,flux_gn,flux_gt, flag, &
                              dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      use functions, only: iminus, iplus, jminus, jplus, kminus, kplus

      use param1, only: zero
      use matrix, only: e, w, n, s, t, b
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_U_g
      real(c_real), intent(inout) :: A_U_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)

      real(c_real), intent(in   ) :: MU_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
      real(c_real), intent(in   ) :: flux_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: flux_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: flux_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K
! Face mass flux
      real(c_real) :: flux_e, flux_w, flux_n, flux_s
      real(c_real) :: flux_t, flux_b
! Diffusion parameter
      real(c_real) :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

!---------------------------------------------------------------------//

      DO K = kstart3, kend3
         DO J = jstart3, jend3
            DO I = istart3, iend3

               IF(flag(i,j,k,2) >= 2000 .and. &
                  flag(i,j,k,2) <= 2011) THEN

! Calculate convection-diffusion fluxes through each of the faces
                  CALL GET_UCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
                     flux_s, flux_t, flux_b, &
                     flux_ge, flux_gn, flux_gt, i, j, k)
                  CALL GET_UCELL_GDIFF_TERMS(&
                     d_fe, d_fw, d_fn, d_fs, &
                     d_ft, d_fb, mu_g, i, j, k, flag, &
                     dx, dy, dz)

! East face (i+1, j, k)
                  IF (Flux_e >= ZERO) THEN
                     A_U_G(I,J,K,E) = D_Fe
                     A_U_G(iplus(i,j,k),j,k,W) = D_Fe + Flux_e
                  ELSE
                     A_U_G(I,J,K,E) = D_Fe - Flux_e
                     A_U_G(iplus(i,j,k),j,k,W) = D_Fe
                  ENDIF
! West face (i, j, k)
                  IF(flag(iminus(i,j,k),j,k,2) < 2000 .or. &
                     flag(iminus(i,j,k),j,k,2) > 2011) THEN
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
                  IF(flag(i,jminus(i,j,k),k,2) < 2000 .or. &
                     flag(i,jminus(i,j,k),k,2) > 2011) THEN
                     IF (Flux_s >= ZERO) THEN
                        A_U_G(I,J,K,S) = D_Fs + Flux_s
                     ELSE
                        A_U_G(I,J,K,S) = D_Fs
                     ENDIF
                  ENDIF

! Top face (i+1/2, j, k+1/2)
                  IF (Flux_t >= ZERO) THEN
                     A_U_G(I,J,K,T) = D_Ft
                     A_U_G(i,j,kplus(i,j,k),B) = D_Ft + Flux_t
                  ELSE
                     A_U_G(I,J,K,T) = D_Ft - Flux_t
                     A_U_G(i,j,kplus(i,j,k),B) = D_Ft
                  ENDIF
! Bottom face (i+1/2, j, k-1/2)
                  IF(flag(i,j,kminus(i,j,k),2) < 2000 .or. &
                     flag(i,j,kminus(i,j,k),2) > 2011) THEN
                     IF (Flux_b >= ZERO) THEN
                        A_U_G(I,J,K,B) = D_Fb + Flux_b
                     ELSE
                        A_U_G(I,J,K,B) = D_Fb
                     ENDIF
                  ENDIF
               ENDIF
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
      SUBROUTINE STORE_A_U_G1(A_U_G, MU_G, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      use functions, only: iminus, iplus, jminus, jplus, kminus, kplus
      use param1   , only: one
      use matrix   , only: e, w, n, s, t, b
      use run      , only: discretize

      use xsi, only: calc_xsi

      IMPLICIT NONE

      real(c_real), intent(in   ) :: MU_G&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      real(c_real), intent(in   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      real(c_real), intent(in   ) :: flux_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: flux_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: flux_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer, intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

      real(c_real), intent(in   ) :: dx, dy, dz

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_U_g
      real(c_real), intent(inout) :: A_U_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)

      real(c_real), intent(in   ) :: dt
! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I,J,K
! indicator for shear
      INTEGER :: incr
! Diffusion parameter
      real(c_real) :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      real(c_real) :: Flux_e, flux_w, flux_n, flux_s
      real(c_real) :: flux_t, flux_b

! x, y, z directional velocity
      real(c_real), allocatable :: U(:,:,:), V(:,:,:), WW(:,:,:)
      real(c_real), allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)
!---------------------------------------------------------------------//


      allocate(  U(istart3:iend3, jstart3:jend3, kstart3:kend3) )
      allocate(  V(istart3:iend3, jstart3:jend3, kstart3:kend3) )
      allocate( WW(istart3:iend3, jstart3:jend3, kstart3:kend3) )

      allocate(xsi_e(istart3:iend3, jstart3:jend3, kstart3:kend3) )
      allocate(xsi_n(istart3:iend3, jstart3:jend3, kstart3:kend3) )
      allocate(xsi_t(istart3:iend3, jstart3:jend3, kstart3:kend3) )

      CALL GET_UCELL_GVTERMS(U, V, WW, u_g, v_g, w_g)

! shear indicator:
      incr=1
      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, XSI_T, incr, &
                     dt, dx, dy, dz)

      DO K = kstart3, kend3
         DO J = jstart3, jend3
            DO I = istart3, iend3

               IF(flag(i,j,k,2) >= 2000 .and. &
                  flag(i,j,k,2) <= 2011) THEN

! Calculate convection-diffusion fluxes through each of the faces
                  CALL GET_UCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
                     flux_s, flux_t, flux_b, &
                     flux_ge, flux_gn, flux_gt, i, j, k)
                  CALL GET_UCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
                     d_ft, d_fb, mu_g, i, j, k, flag, &
                     dx, dy, dz)

! East face (i+1, j, k)
                  A_U_G(I,J,K,E) = D_Fe - XSI_E(i,j,k) * Flux_e
                  A_U_G(iplus(i,j,k),j,k,W) = D_Fe + flux_e*&
                     (ONE - XSI_E(i,j,k))
! West face (i, j, k)
                  IF(flag(iminus(i,j,k),j,k,2) < 2000 .or. &
                     flag(iminus(i,j,k),j,k,2) > 2011) THEN
                     A_U_G(I,J,K,W) = D_Fw + flux_w*&
                        (ONE - XSI_E(iminus(i,j,k),j,k))
                  ENDIF


! North face (i+1/2, j+1/2, k)
                  A_U_G(I,J,K,N) = D_Fn - XSI_N(i,j,k) * Flux_n
                  A_U_G(i,jplus(i,j,k),k,S) = D_Fn + flux_n*&
                     (ONE - XSI_N(i,j,k))
! South face (i+1/2, j-1/2, k)
                  IF(flag(i,jminus(i,j,k),k,2) < 2000 .or. &
                     flag(i,jminus(i,j,k),k,2) > 2011) THEN
                     A_U_G(I,J,K,S) = D_Fs + flux_s*&
                        (ONE - XSI_N(i,jminus(i,j,k),k))
                  ENDIF

! Top face (i+1/2, j, k+1/2)
                  A_U_G(I,J,K,T) = D_Ft - XSI_T(i,j,k) * Flux_t
                  A_U_G(i,j,kplus(i,j,k),B) = D_Ft + flux_t*&
                     (ONE - XSI_T(i,j,k))
! Bottom face (i+1/2, j, k-1/2)
                  IF(flag(i,j,kminus(i,j,k),2) < 2000 .or. &
                     flag(i,j,kminus(i,j,k),2) > 2011) THEN
                     A_U_G(I,J,K,B) = D_Fb + flux_b*&
                        (ONE - XSI_T(i,j,kminus(i,j,k)))
                  ENDIF

               ENDIF
            ENDDO
         ENDDO
      ENDDO

      deallocate( U, V, WW )
      deallocate( xsi_e, xsi_n, xsi_t)


      RETURN
      END SUBROUTINE STORE_A_U_G1
end module u_g_conv_dif
