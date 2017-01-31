module w_g_conv_dif

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use geometry      , only: domlo, domhi
   use param1        , only: half, one, zero

   implicit none

   private
   public :: conv_dif_w_g

   contains
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
      SUBROUTINE CONV_DIF_W_G(slo, shi, lo, hi, A_m, MU_G, u_g, v_g, w_g, &
                              flux_ge, flux_gn, flux_gt, flag, dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE run, only: discretize

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(IN   ) :: MU_G&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer     , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dt, dx, dy, dz
!---------------------------------------------------------------------//

      IF (DISCRETIZE(5) == 0) THEN               ! 0 & 1 => FOUP
         CALL STORE_A_W_G0 (slo, shi, lo, hi, A_m, MU_G, flux_ge, flux_gn, flux_gt, flag, &
                            dx, dy, dz)
      ELSE
         CALL STORE_A_W_G1 (slo, shi, lo, hi, A_m, MU_G, u_g, v_g, w_g, &
                            flux_ge, flux_gn, flux_gt,flag, dt, dx, dy, dz)
      ENDIF

      END SUBROUTINE CONV_DIF_W_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a w-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_WCELL_GVTERMS(slo, shi, U, V, WW, u_g, v_g, w_g)

      USE functions, only: avg
      USE functions, only: kplus

      integer     , intent(in   ) :: slo(3),shi(3)

      real(c_real), intent(OUT) :: U&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(OUT) :: V&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(OUT) :: WW&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent( in) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent( in) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent( in) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K
!---------------------------------------------------------------------//


      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)
            U(I,J,K)  = avg(u_g(I,J,K),U_G(i,j,kplus(i,j,k)))
            V(I,J,K)  = avg(v_g(I,J,K),V_G(i,j,kplus(i,j,k)))
            WW(I,J,K) = avg(w_g(I,J,K),W_G(i,j,kplus(i,j,k)))
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GET_WCELL_GVTERMS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a w-momentum cell. Note the fluxes are calculated at       C
!  all faces.                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_WCELL_GDIFF_TERMS(&
         slo, shi, lo, hi, &
         D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, MU_G, I, J, K, flag, &
         dx, dy, dz)

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) ::  lo(3), hi(3)

      ! diffusion through faces of given w-momentum cell
      real(c_real), intent(OUT) :: d_fe, d_fw
      real(c_real), intent(OUT) :: d_fn, d_fs
      real(c_real), intent(OUT) :: d_ft, d_fb

      real(c_real), intent(IN   ) :: MU_G&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, intent(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      INTEGER, intent(IN) :: i, j, k
      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: kp, im, jm, kc
      INTEGER :: itmp, jtmp, ktmp
! length terms
      real(c_real) :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB

      real(c_real) :: odx, ody, odz
      real(c_real) :: axy, axz, ayz
!---------------------------------------------------------------------//

      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      im = min(domlo(1)-1,i-1)
      jm = min(domlo(2)-1,j-1)
      kp = max(domhi(3)+1,k+1)

      ktmp = ktop(i,j,k)
      itmp = i-1
      jtmp  = jsouth(i,j,k)

      IF (flag(i,j,k,1)>=100) THEN
         kc = ktop(i,j,k)
      ELSE
         kc = k
      ENDIF

      C_AE = ODX
      C_AW = ODX
      C_AN = ODY
      C_AS = ODY
      C_AT = ODZ
      C_AB = ODZ

! East face (i+1/2, j, k+1/2)
      D_Fe = AVG_H(AVG_H(MU_G(i,j,kc),MU_G(ieast(i,j,k),j,k)),&
                   AVG_H(MU_G(i,j,ktmp),MU_G(ieast(i,j,ktmp),j,ktmp)))*C_AE*AYZ
! West face (i-1/2, j, k+1/2)
      D_Fw = AVG_H(AVG_H(MU_G(itmp,j,k),MU_G(i,j,kc)),&
                   AVG_H(MU_G(itmp,j,ktop(itmp,j,k)),MU_G(i,j,ktmp)))*C_AW*AYZ

! North face (i, j+1/2, k+1/2)
      D_Fn = AVG_H(AVG_H(MU_G(i,j,kc),MU_G(i,jnorth(i,j,k),k)),&
                   AVG_H(MU_G(i,j,ktmp),MU_G(i,jnorth(i,j,ktmp),ktmp)))*C_AN*AXZ
! South face (i, j-1/2, k+1/2)
      D_Fs = AVG_H(AVG_H(MU_G(i,jtmp,k),MU_G(i,j,kc)),&
                   AVG_H(MU_G(i,jtmp,ktop(i,jtmp,k)),MU_G(i,j,ktmp)))*C_AS*AXZ

! Top face (i, j, k+1)
      D_Ft = MU_G(i,j,ktmp)*C_AT*AXY
! Bottom face (i, j, k)
      D_Fb = MU_G(i,j,k)*C_AB*AXY

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

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
      SUBROUTINE STORE_A_W_G0(slo, shi, lo, hi, &
                              A_W_G, MU_G, flux_ge, flux_gn, flux_gt, flag, &
                              dx, dy, dz)

! Modules
!---------------------------------------------------------------------//

      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus

      use matrix, only: e, w, n, s, t, b

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_U_g
      real(c_real), intent(inout) :: A_W_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(in   ) :: MU_G&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer     , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dx, dy, dz
! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K
! Face mass flux
      real(c_real) :: flux_e, flux_n, flux_t
! Diffusion parameter
      real(c_real) :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

!---------------------------------------------------------------------//

      do k = slo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

! Calculate convection-diffusion fluxes through each of the faces
               flux_e = HALF * (flux_gE(i,j,k) + flux_gE(i,j,k+1))
               flux_n = HALF * (flux_gN(i,j,k) + flux_gN(i,j,k+1))
               flux_t = HALF * (flux_gT(i,j,k) + flux_gT(i,j,k+1))

               CALL GET_WCELL_GDIFF_TERMS(&
                  slo, shi, lo, hi, &
                  d_fe, d_fw, d_fn, d_fs, &
                  d_ft, d_fb, mu_g, i, j, k, flag, dx, dy, dz)

! East face (i+1/2, j, k+1/2)
               if (flux_e >= zero) then
                  a_w_g(i,  j,k,e) = d_fe
                  a_w_g(i+1,j,k,w) = d_fe + flux_e
               else
                  a_w_g(i,  j,k,e) = d_fe - flux_e
                  a_w_g(i+1,j,k,w) = d_fe
               endif

! North face (i, j+1/2, k+1/2)
               if (flux_n >= zero) then
                  a_w_g(i,j,  k,n) = d_fn
                  a_w_g(i,j+1,k,s) = d_fn + flux_n
               else
                  a_w_g(i,j,  k,n) = d_fn - flux_n
                  a_w_g(i,j+1,k,s) = d_fn
               endif

! Top face (i, j, k+1)
               if (flux_t >= zero) then
                  a_w_g(i,j,k,  t) = d_ft
                  a_w_g(i,j,k+1,b) = d_ft + flux_t
               else
                  a_w_g(i,j,k,  t) = d_ft - flux_t
                  a_w_g(i,j,k+1,b) = d_ft
               endif

! West face (i-1/2, j, k+1/2)
               if(i==lo(1)) a_w_g(i,j,k,w) = d_fw

! South face (i, j-1/2, k+1/2)
               if(j==lo(2)) a_w_g(i,j,k,s) = d_fs

            enddo
         enddo
      enddo

      return
   end subroutine store_a_w_g0


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
      SUBROUTINE STORE_A_W_G1(slo, shi, lo, hi, A_W_G, MU_G, u_g, v_g, w_g, &
                              flux_ge, flux_gn, flux_gt, flag,  &
                              dt, dx, dy, dz)

      USE functions, only: iplus, iminus, jplus, jminus, kplus, kminus

      use matrix, only: e, w, n, s, t, b

      USE run, only: discretize

      USE xsi, only: calc_xsi

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_W_g
      real(c_real), intent(INOUT) :: A_W_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(IN   ) :: MU_G&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(IN   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(IN   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(IN   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, intent(IN) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I,J,K
! Diffusion parameter
      real(c_real) :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      real(c_real) :: flux_e, flux_w, flux_n, flux_s
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

      CALL GET_WCELL_GVTERMS(slo, shi, U, V, WW, u_g, v_g, w_g )


      call calc_xsi (discretize(3), slo, shi, hi, &
         w_g, u, v, ww, xsi_e, xsi_n, xsi_t, &
         dt, dx, dy, dz)

      do k = slo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               flux_e = half * (flux_ge(i,j,k) + flux_ge(i,j,k+1))
               flux_n = half * (flux_gn(i,j,k) + flux_gn(i,j,k+1))
               flux_t = half * (flux_gt(i,j,k) + flux_gt(i,j,k+1))

               CALL GET_WCELL_GDIFF_TERMS(&
                  slo, shi, lo, hi, &
                  d_fe, d_fw, d_fn, d_fs, &
                  d_ft, d_fb, mu_g, i, j, k, flag, dx, dy, dz)

! East face (i+1/2, j, k+1/2)
               a_w_g(i,  j,k,e) = d_fe - flux_e*(xsi_e(i,j,k))
               a_w_g(i+1,j,k,w) = d_fe + flux_e*(one - xsi_e(i,j,k))

! North face (i, j+1/2, k+1/2)
               a_w_g(i,j,  k,n) = d_fn - flux_n*(xsi_n(i,j,k))
               a_w_g(i,j+1,k,s) = d_fn + flux_n*(one - xsi_n(i,j,k))

! Top face (i, j, k+1)
               a_w_g(i,j,k,  t) = d_ft - flux_t*(xsi_t(i,j,k))
               a_w_g(i,j,k+1,b) = d_ft + flux_t*(one - xsi_t(i,j,k))

! West face (i-1/2, j, k+1/2)
               if(i==lo(1)) a_w_g(i,j,k,w) = d_fw
! South face (i, j-1/2, k+1/2)
               if(j==lo(2)) a_w_g(i,j,k,s) = d_fs

            enddo
         enddo
      enddo

      deallocate( U, V, WW )
      deallocate( xsi_e, xsi_n, xsi_t)

      return
   end subroutine store_a_w_g1
end module w_g_conv_dif
