module u_g_conv_dif

   use bl_fort_module, only: c_real
   use geometry      , only: domlo, domhi
   use param1        , only: half, one, zero

   implicit none

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
   subroutine conv_dif_u_g(&
      slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
      A_M, mu_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt,&
      flag, dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use run, only: discretize

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      real(c_real), intent(inout) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer    , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!---------------------------------------------------------------------//

      if (discretize(3) == 0) then
         call store_a_u_g0(&
            slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
            A_m, mu_g, flux_ge, flux_gn, flux_gt, flag, dx, dy, dz)
      else
         call store_a_u_g1(&
            slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
            A_m, mu_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt, &
            flag, dt, dx, dy, dz)
      ENDIF

      END SUBROUTINE CONV_DIF_U_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a u-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GVTERMS(slo, shi, lo, hi, U, V, WW, u_g, v_g, w_g)

      use functions, only: avg
      use functions, only: iplus

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      real(c_real), intent(OUT) :: U&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(OUT) :: V&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(OUT) :: WW&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I,J,K
!---------------------------------------------------------------------//

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)
            U(I,J,K) = AVG(U_G(I,J,K), U_G(iplus(i,j,k),j,k))
            V(I,J,K) = AVG(V_G(I,J,K), V_G(iplus(i,j,k),j,k))
            WW(I,J,K) = AVG(W_G(I,J,K), W_G(iplus(i,j,k),j,k))
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE GET_UCELL_GVTERMS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a u-momentum cell. Note the fluxes are calculated at       C
!  all faces.                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GDIFF_TERMS(&
         slo, shi, lo, hi, &
         D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, MU_G, I, J, K, flag, dx, dy, dz)

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! diffusion through faces of given ijk u-momentum cell
      real(c_real), intent(OUT) :: d_fe, d_fw
      real(c_real), intent(OUT) :: d_fn, d_fs
      real(c_real), intent(OUT) :: d_ft, d_fb

      real(c_real), intent( IN) :: MU_G&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, intent( IN) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

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

      IP = min(domhi(1)+1, i+1)
      JM = max(domlo(2)-1, j-1)
      KM = max(domlo(3)-1, k-1)

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
   subroutine store_a_u_g0(&
      slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
      A_U_g, mu_g, flux_ge, flux_gn, flux_gt, flag, dx, dy, dz)

      use functions, only: iminus, iplus, jminus, jplus, kminus, kplus

      use matrix, only: e, w, n, s, t, b

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      real(c_real), intent(in   ) :: dx, dy, dz

      real(c_real), intent(inout) :: A_U_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(in   ) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer, intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

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

      DO K = lo(3),hi(3)
        DO J = lo(2),hi(2)
          DO I = slo(1),hi(1)

               IF(flag(i,j,k,2) >= 2000 .and. &
                  flag(i,j,k,2) <= 2011) THEN

! Calculate convection-diffusion fluxes through each of the faces
                  flux_e = HALF * (flux_gE(i  ,j,k) + flux_gE(i+1,j,k))
                  if (i.eq.domlo(1)-1) then
                     flux_w = flux_gE(i  ,j,k)
                  else
                     flux_w = HALF * (flux_gE(i-1,j,k) + flux_gE(i  ,j,k))
                  end if
                  flux_n = HALF * (flux_gN(i,j  ,k) + flux_gN(i+1,j  ,k))
                  flux_s = HALF * (flux_gN(i,j-1,k) + flux_gN(i+1,j-1,k))
                  flux_t = HALF * (flux_gT(i,j,k  ) + flux_gT(i+1,j,k  ))
                  flux_b = HALF * (flux_gT(i,j,k-1) + flux_gT(i+1,j,k-1))
                  CALL GET_UCELL_GDIFF_TERMS(&
                     slo, shi, lo, hi, &
                     d_fe, d_fw, d_fn, d_fs, &
                     d_ft, d_fb, mu_g, i, j, k, flag, &
                     dx, dy, dz)

! East face (i+1, j, k)
                  IF (flux_e >= ZERO) THEN
                     A_U_G(I,J,K,E) = D_Fe
                     A_U_G(iplus(i,j,k),j,k,W) = D_Fe + flux_e
                  ELSE
                     A_U_G(I,J,K,E) = D_Fe - flux_e
                     A_U_G(iplus(i,j,k),j,k,W) = D_Fe
                  ENDIF
! West face (i, j, k)
                  IF(flag(iminus(i,j,k),j,k,2) < 2000 .or. &
                     flag(iminus(i,j,k),j,k,2) > 2011) THEN
                     IF (flux_w >= ZERO) THEN
                        A_U_G(I,J,K,W) = D_Fw + flux_w
                     ELSE
                        A_U_G(I,J,K,W) = D_Fw
                     ENDIF
                  ENDIF

! North face (i+1/2, j+1/2, k)
                  IF (flux_n >= ZERO) THEN
                     A_U_G(I,J,K,N) = D_Fn
                     A_U_G(i,jplus(i,j,k),k,S) = D_Fn + flux_n
                  ELSE
                     A_U_G(I,J,K,N) = D_Fn - flux_n
                     A_U_G(i,jplus(i,j,k),k,S) = D_Fn
                  ENDIF

! South face (i+1/2, j-1/2, k)
                  IF(flag(i,jminus(i,j,k),k,2) < 2000 .or. &
                     flag(i,jminus(i,j,k),k,2) > 2011) THEN
                     IF (flux_s >= ZERO) THEN
                        A_U_G(I,J,K,S) = D_Fs + flux_s
                     ELSE
                        A_U_G(I,J,K,S) = D_Fs
                     ENDIF
                  ENDIF

! Top face (i+1/2, j, k+1/2)
                  IF (flux_t >= ZERO) THEN
                     A_U_G(I,J,K,T) = D_Ft
                     A_U_G(i,j,kplus(i,j,k),B) = D_Ft + flux_t
                  ELSE
                     A_U_G(I,J,K,T) = D_Ft - flux_t
                     A_U_G(i,j,kplus(i,j,k),B) = D_Ft
                  ENDIF
! Bottom face (i+1/2, j, k-1/2)
                  IF(flag(i,j,kminus(i,j,k),2) < 2000 .or. &
                     flag(i,j,kminus(i,j,k),2) > 2011) THEN
                     IF (flux_b >= ZERO) THEN
                        A_U_G(I,J,K,B) = D_Fb + flux_b
                     ELSE
                        A_U_G(I,J,K,B) = D_Fb
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      return
   end subroutine store_a_u_g0

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
      SUBROUTINE STORE_A_U_G1(&
         slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
         A_U_g, mu_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt, &
         flag, dt, dx, dy, dz)

      use functions, only: iminus, iplus, jminus, jplus, kminus, kplus
      use matrix   , only: e, w, n, s, t, b
      use run      , only: discretize

      use xsi, only: calc_xsi

! Dummy arguments
!---------------------------------------------------------------------//
      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      real(c_real), intent(in   ) :: dx, dy, dz, dt


      real(c_real), intent(inout) :: A_U_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer, intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)


! Local variables
!---------------------------------------------------------------------//
! Indices
      integer :: i,j,k
! indicator for shear
      integer :: incr
! Diffusion parameter
      real(c_real) :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      real(c_real) :: flux_e, flux_w, flux_n, flux_s
      real(c_real) :: flux_t, flux_b

! x, y, z directional velocity
      real(c_real), allocatable :: U(:,:,:), V(:,:,:), WW(:,:,:)
      real(c_real), allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)

      allocate(  U(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)) )
      allocate(  V(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)) )
      allocate( WW(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)) )

      allocate(xsi_e(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)) )
      allocate(xsi_n(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)) )
      allocate(xsi_t(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)) )

      CALL GET_UCELL_GVTERMS(slo, shi, lo, hi, U, V, WW, u_g, v_g, w_g)

! shear indicator:
      incr=1
      call calc_xsi (discretize(3), slo, shi, lo, hi, &
         u_g, u, v, ww, xsi_e, xsi_n, xsi_t, &
         dt, dx, dy, dz)

      do K = lo(3),hi(3)
         do J = lo(2),hi(2)
            do I = slo(1),hi(1)

               IF(flag(i,j,k,2) >= 2000 .and. &
                  flag(i,j,k,2) <= 2011) THEN

! Calculate convection-diffusion fluxes through each of the faces
                  flux_e = HALF * (flux_gE(i  ,j,k) + flux_gE(i+1,j,k))
                  flux_n = HALF * (flux_gN(i,j  ,k) + flux_gN(i+1,j  ,k))
                  flux_s = HALF * (flux_gN(i,j-1,k) + flux_gN(i+1,j-1,k))
                  flux_t = HALF * (flux_gT(i,j,k  ) + flux_gT(i+1,j,k  ))
                  flux_b = HALF * (flux_gT(i,j,k-1) + flux_gT(i+1,j,k-1))
                  if (i.eq.domlo(1)-1) then
                     flux_w = flux_gE(i  ,j,k)
                  else
                     flux_w = HALF * (flux_gE(i-1,j,k) + flux_gE(i  ,j,k))
                  end if
                  CALL GET_UCELL_GDIFF_TERMS(&
                     slo, shi, lo, hi, &
                     d_fe, d_fw, d_fn, d_fs, &
                     d_ft, d_fb, mu_g, i, j, k, flag, &
                     dx, dy, dz)

! East face (i+1, j, k)
                  A_U_G(I,J,K,E) = D_Fe - XSI_E(i,j,k) * flux_e
                  A_U_G(iplus(i,j,k),j,k,W) = D_Fe + flux_e*&
                     (ONE - XSI_E(i,j,k))
! West face (i, j, k)
                  IF(flag(iminus(i,j,k),j,k,2) < 2000 .or. &
                     flag(iminus(i,j,k),j,k,2) > 2011) THEN
                     A_U_G(I,J,K,W) = D_Fw + flux_w*&
                        (ONE - XSI_E(iminus(i,j,k),j,k))
                  ENDIF


! North face (i+1/2, j+1/2, k)
                  A_U_G(I,J,K,N) = D_Fn - XSI_N(i,j,k) * flux_n
                  A_U_G(i,jplus(i,j,k),k,S) = D_Fn + flux_n*&
                     (ONE - XSI_N(i,j,k))
! South face (i+1/2, j-1/2, k)
                  IF(flag(i,jminus(i,j,k),k,2) < 2000 .or. &
                     flag(i,jminus(i,j,k),k,2) > 2011) THEN
                     A_U_G(I,J,K,S) = D_Fs + flux_s*&
                        (ONE - XSI_N(i,jminus(i,j,k),k))
                  ENDIF

! Top face (i+1/2, j, k+1/2)
                  A_U_G(I,J,K,T) = D_Ft - XSI_T(i,j,k) * flux_t
                  A_U_G(i,j,kplus(i,j,k),B) = D_Ft + flux_t*&
                     (ONE - XSI_T(i,j,k))
! Bottom face (i+1/2, j, k-1/2)
                  IF(flag(i,j,kminus(i,j,k),2) < 2000 .or. &
                     flag(i,j,kminus(i,j,k),2) > 2011) THEN
                     A_U_G(I,J,K,B) = D_Fb + flux_b*&
                        (ONE - XSI_T(i,j,kminus(i,j,k)))
                  endif

               endif
            enddo
         enddo
      enddo

      deallocate( u, v, ww )
      deallocate( xsi_e, xsi_n, xsi_t)


      RETURN
      END SUBROUTINE STORE_A_U_G1
end module u_g_conv_dif
