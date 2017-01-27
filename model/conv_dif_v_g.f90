module v_g_conv_dif

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   private
   public :: conv_dif_v_g

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_V_g(A_m)                                 C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_DIF_V_G(slo, shi, lo, hi, A_M, mu_g, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE run, only: discretize

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), INTENT(IN   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), INTENT(IN   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), INTENT(IN   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, INTENT(IN) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dt, dx, dy, dz
!---------------------------------------------------------------------//

! DO NOT USE DEFERRED CORRECTION TO SOLVE V_G
      IF (DISCRETIZE(4) == 0) THEN               ! 0 & 1 => FOUP
         CALL STORE_A_V_G0(slo,shi,lo,hi,A_M,MU_G,flux_ge,flux_gn,flux_gt,flag,dx,dy,dz)
      ELSE
         CALL STORE_A_V_G1(slo,shi,lo,hi,A_M,MU_G,u_g,v_g,w_g,flux_ge,flux_gn,flux_gt,flag,dt,dx,dy,dz)
      ENDIF

      RETURN
      END SUBROUTINE CONV_DIF_V_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a v-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_VCELL_GVTERMS(slo, shi, lo, hi, U, V, WW, u_g, v_g, w_g)

! Modules
!---------------------------------------------------------------------//

      USE functions, only: jplus
      USE functions, only: avg

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      real(c_real), INTENT(OUT) :: U&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(OUT) :: V&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(OUT) :: WW&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), INTENT( in) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT( in) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT( in) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K
!---------------------------------------------------------------------//

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)
            U(I,J,K) = AVG(U_G(I,J,K),U_G(i,jplus(i,j,k),k))
            V(I,J,K) = AVG(V_G(I,J,K),V_G(i,jplus(i,j,k),k))
            WW(I,J,K) = AVG(W_G(I,J,K),W_G(i,jplus(i,j,k),k))
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GET_VCELL_GVTERMS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the convective fluxes through the faces of a     C
!  v-momentum cell. Note the fluxes are calculated at all faces.       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_VCELL_GCFLUX_TERMS(&
         slo, shi, lo, hi, &
         FLUX_E, FLUX_W, FLUX_N, &
         FLUX_S, FLUX_T, FLUX_B, &
         flux_ge, flux_gn, flux_gt, i, j, k)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: iplus, iminus, jplus, jminus, kplus, kminus

      USE param1, only: half
      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! fluxes through faces of given u-momentum cell
      real(c_real), INTENT(OUT) :: flux_e, flux_w
      real(c_real), INTENT(OUT) :: flux_n, flux_s
      real(c_real), INTENT(OUT) :: flux_t, flux_b

      real(c_real), INTENT(IN   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
! indices
      INTEGER, INTENT(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//

      Flux_e = HALF * (Flux_gE(i,j,k) + Flux_gE(i,jplus(i,j,k),k))
      Flux_w = HALF * (Flux_gE(iminus(i,j,k),j,k) + &
                       Flux_gE(iminus(i,j,k),jplus(iminus(i,j,k),j,k),k))

      Flux_n = HALF * (Flux_gN(i,j,k) + Flux_gN(i,jplus(i,j,k),k))
      Flux_s = HALF * (Flux_gN(i,jminus(i,j,k),k) + Flux_gN(i,j,k))

      Flux_t = HALF * (Flux_gT(i,j,k) + Flux_gT(i,jplus(i,j,k),k))
      Flux_b = HALF * (Flux_gT(i,j,kminus(i,j,k)) + &
                       Flux_gT(i,jplus(i,j,kminus(i,j,k)),kminus(i,j,k)))

      END SUBROUTINE GET_VCELL_GCFLUX_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  the faces of a v-momentum cell. Note the fluxes are calculated at   C
!  all faces.                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_VCELL_GDIFF_TERMS(&
         slo, shi, lo, hi, &
         D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, MU_G, I, J, K,flag, &
         dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero
      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! diffusion through faces of given v-momentum cell
      real(c_real), INTENT(  OUT) :: d_fe, d_fw
      real(c_real), INTENT(  OUT) :: d_fn, d_fs
      real(c_real), INTENT(  OUT) :: d_ft, d_fb

      real(c_real), INTENT(IN   ) :: MU_G&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), INTENT(IN   ) :: dx, dy, dz

      INTEGER, INTENT(IN) :: i, j, k

! Local variables
      INTEGER :: jc, jn
      real(c_real) :: odx, ody, odz
      real(c_real) :: axy, axz, ayz
!---------------------------------------------------------------------//
      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      jn = jnorth(i,j,k)

      IF (flag(i,j,k,1)>=100) THEN
         jc = jn
      ELSE
         jc = j
      ENDIF

! East face (i+1/2, j+1/2, k)
      D_Fe = AVG_H(AVG_H(MU_G(i,jc,k),MU_G(ieast(i,j,k),j,k)),&
                   AVG_H(MU_G(i,jn,k), &
                         MU_G(ieast(i,jn,k),jn,k)))*ODX*AYZ
! West face (i-1/2, j+1/2, k)
      D_Fw = AVG_H(AVG_H(MU_G(i-1,j,k),MU_G(i,jc,k)),&
                   AVG_H(MU_G(i-1,jnorth(i-1,j,k),k),&
                         MU_G(i           ,jn,k)))*ODX*AYZ

      ! North face (i, j+1, k)
      D_Fn = MU_G(i,jn,k)*ODY*AXZ

      ! South face (i, j, k)
      D_Fs = MU_G(i,jc,k)*ODY*AXZ

      D_FT = ZERO
      D_FB = ZERO

! Top face (i, j+1/2, k+1/2)
      D_Ft = AVG_H(AVG_H(MU_G(i,jc,k),MU_G(i,j,ktop(i,j,k))),&
             AVG_H(MU_G(i,jn,k),MU_G(i,jnorth(i,j,ktop(i,j,k)),ktop(i,j,k))))*ODZ*AXY
! Bottom face (i, j+1/2, k-1/2)
      D_Fb = AVG_H(AVG_H(MU_G(i,j,kbot(i,j,k)),MU_G(i,jc,k)),&
             AVG_H(MU_G(i,jnorth(i,j,kbot(i,j,k)),kbot(i,j,k)),MU_G(i,jn,k)))*ODZ*AXY

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

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
      SUBROUTINE STORE_A_V_G0(&
         slo, shi, lo, hi, &
         A_V_G,mu_g,flux_ge,flux_gn,flux_gt,flag,dx,dy,dz)

! Modules
!---------------------------------------------------------------------//
      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus
      USE param1, only: zero
      use matrix, only: e, w, n, s, t, b

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_V_g
      real(c_real), INTENT(INOUT) :: A_V_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), INTENT(IN   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), INTENT(IN   ) :: dx,dy,dz

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

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)

               IF(flag(i,j,k,3) >= 2000 .and. &
                  flag(i,j,k,3) <= 2011) THEN

! Calculate convection-diffusion fluxes through each of the faces
                  CALL GET_VCELL_GCFLUX_TERMS(&
                     slo, shi, lo, hi, &
                     flux_e, flux_w, flux_n, &
                     flux_s, flux_t, flux_b, &
                     flux_ge, flux_gn, flux_gt, i, j, k)
                  CALL GET_VCELL_GDIFF_TERMS(&
                     slo, shi, lo, hi, &
                     d_fe, d_fw, d_fn, d_fs, &
                     d_ft, d_fb, mu_g, i, j, k, flag, dx, dy, dz)

! East face (i+1/2, j+1/2, k)
                  IF (Flux_e >= ZERO) THEN
                     A_V_G(I,J,K,E) = D_Fe
                     A_V_G(iplus(i,j,k),j,k,W) = D_Fe + Flux_e
                  ELSE
                     A_V_G(I,J,K,E) = D_Fe - Flux_e
                     A_V_G(iplus(i,j,k),j,k,W) = D_Fe
                  ENDIF
! West face (i-1/2, j+1/2, k)
                  IF(flag(iminus(i,j,k),j,k,3) < 2000 .or. &
                     flag(iminus(i,j,k),j,k,3) > 2011) THEN
                     IF (Flux_w >= ZERO) THEN
                        A_V_G(I,J,K,W) = D_Fw + Flux_w
                     ELSE
                        A_V_G(I,J,K,W) = D_Fw
                     ENDIF
                  ENDIF


! North face (i, j+1, k)
                  IF (Flux_n >= ZERO) THEN
                     A_V_G(I,J,K,N) = D_Fn
                     A_V_G(i,jplus(i,j,k),k,S) = D_Fn + Flux_n
                  ELSE
                     A_V_G(I,J,K,N) = D_Fn - Flux_n
                     A_V_G(i,jplus(i,j,k),k,S) = D_Fn
                  ENDIF
! South face (i, j, k)
                  IF(flag(i,jminus(i,j,k),k,3) < 2000 .or. &
                     flag(i,jminus(i,j,k),k,3) > 2011) THEN
                     IF (Flux_s >= ZERO) THEN
                        A_V_G(I,J,K,S) = D_Fs + Flux_s
                     ELSE
                        A_V_G(I,J,K,S) = D_Fs
                     ENDIF
                  ENDIF


! Top face (i, j+1/2, k+1/2)
                  IF (Flux_t >= ZERO) THEN
                     A_V_G(I,J,K,T) = D_Ft
                     A_V_G(i,j,kplus(i,j,k),B) = D_Ft + Flux_t
                  ELSE
                     A_V_G(I,J,K,T) = D_Ft - Flux_t
                     A_V_G(i,j,kplus(i,j,k),B) = D_Ft
                  ENDIF
! Bottom face (i, j+1/2, k-1/2)
                  IF(flag(i,j,kminus(i,j,k),3) < 2000 .or. &
                     flag(i,j,kminus(i,j,k),3) > 2011) THEN
                     IF (Flux_b >= ZERO) THEN
                        A_V_G(I,J,K,B) = D_Fb + Flux_b
                     ELSE
                        A_V_G(I,J,K,B) = D_Fb
                     ENDIF
                  ENDIF
               ENDIF
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
      SUBROUTINE STORE_A_V_G1(&
         slo, shi, lo, hi, A_V_G, mu_g, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt,flag, dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//

      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus
      USE param1, only: one
      use matrix, only: e, w, n, s, t, b
      USE run, only: discretize
      USE xsi, only: calc_xsi

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_V_g
      real(c_real), INTENT(INOUT) :: A_V_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), INTENT(IN   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT( in) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT( in) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT( in) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K
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

      allocate(  U(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate(  V(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate( WW(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )

      allocate(xsi_e(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate(xsi_n(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate(xsi_t(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )

      CALL GET_VCELL_GVTERMS(slo, shi, lo, hi, U, V, WW, u_g, v_g, w_g)

! shear indicator: y-momentum
      incr=2
      CALL CALC_XSI (DISCRETIZE(4), slo, shi, lo, hi, &
                     V_G, U, V, WW, XSI_E, XSI_N, XSI_T, incr, &
                     dt, dx, dy, dz)

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)

               IF(flag(i,j,k,3) >= 2000 .and. &
                  flag(i,j,k,3) <= 2011) THEN

! Calculate convection-diffusion fluxes through each of the faces
                  CALL GET_VCELL_GCFLUX_TERMS(&
                     slo, shi, lo, hi, &
                     flux_e, flux_w, flux_n, &
                     flux_s, flux_t, flux_b, &
                     flux_ge, flux_gn, flux_gt, i, j, k)
                  CALL GET_VCELL_GDIFF_TERMS(&
                     slo, shi, lo, hi, &
                     d_fe, d_fw, d_fn, d_fs, &
                     d_ft, d_fb, mu_g, i, j, k, flag, dx, dy, dz)

! East face (i+1/2, j+1/2, k)
                  A_V_G(I,J,K,E) = D_Fe - XSI_E(i,j,k)*Flux_e
                  A_V_G(iplus(i,j,k),j,k,W) = D_Fe + flux_e*&
                     (ONE - XSI_E(i,j,k))
! West face (i-1/2, j+1/2, k)
                  IF(flag(iminus(i,j,k),j,k,3) < 2000 .or. &
                     flag(iminus(i,j,k),j,k,3) > 2011) THEN
                     A_V_G(I,J,K,W) = D_Fw + flux_w*&
                        (ONE - XSI_E(iminus(i,j,k),j,k))
                  ENDIF

! North face (i, j+1, k)
                  A_V_G(I,J,K,N) = D_Fn - XSI_N(i,j,k)*Flux_n
                  A_V_G(i,jplus(i,j,k),k,S) = D_Fn + flux_n*&
                     (ONE - XSI_N(i,j,k))
! South face (i, j, k)
                  IF(flag(i,jminus(i,j,k),k,3) < 2000 .or. &
                     flag(i,jminus(i,j,k),k,3) > 2011) THEN
                     A_V_G(I,J,K,S) = D_Fs + flux_s*&
                        (ONE - XSI_N(i,jminus(i,j,k),k))
                  ENDIF

! Top face (i, j+1/2, k+1/2)
                  A_V_G(I,J,K,T) = D_Ft - XSI_T(i,j,k)*Flux_t
                  A_V_G(i,j,kplus(i,j,k),B) = D_Ft + flux_t*&
                     (ONE - XSI_T(i,j,k))
! Bottom face (i, j+1/2, k-1/2)
                  IF(flag(i,j,kminus(i,j,k),3) < 2000 .or. &
                     flag(i,j,kminus(i,j,k),3) > 2011) THEN
                     A_V_G(I,J,K,B) = D_Fb + flux_b*&
                        (ONE - XSI_T(i,j,kminus(i,j,k)))
                  ENDIF

               ENDIF
            ENDDO
         ENDDO
      ENDDO

      deallocate( U, V, WW )
      deallocate( xsi_e, xsi_n, xsi_t)

      RETURN
      END SUBROUTINE STORE_A_V_G1
end module v_g_conv_dif
