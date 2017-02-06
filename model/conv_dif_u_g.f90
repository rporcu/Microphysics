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
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a u-momentum cell. Note the fluxes are calculated at       C
!  all faces.                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_GDIFF_TERMS(&
         slo, shi, &
         D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, mu_g, I, J, K, flag, dx, dy, dz)

      use functions, only: avg, avg_h

      integer     , intent(in   ) :: slo(3),shi(3)

      ! diffusion through faces of given ijk u-momentum cell
      real(c_real), intent(OUT) :: d_fe, d_fw
      real(c_real), intent(OUT) :: d_fn, d_fs
      real(c_real), intent(OUT) :: d_ft, d_fb

      real(c_real), intent( IN) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, intent( IN) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in) :: dx, dy, dz

! ijk index
      INTEGER, intent(IN) :: i, j, k

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: ic
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

!     IF (flag(i,j,k,1)>=100)  THEN
!        IC = ieast(i,j,k)
!     ELSE
!        IC = i
!     ENDIF

      ic = i

      C_AE = ODX
      C_AW = ODX
      C_AN = ODY
      C_AS = ODY
      C_AT = ODZ
      C_AB = ODZ

      ! East face (i+1, j, k)
      D_FE = mu_g(i+1,j,k)*C_AE*AYZ

      ! West face (i, j, k)
      D_FW = mu_g(ic,j,k)*C_AW*AYZ

      ! North face (i+1/2, j+1/2, k)
      D_FN = avg_h(avg_h(mu_g(IC,J,K),mu_g(i,j+1,k)),&
                   avg_h(mu_g(i+1,j,k),mu_g(i+1,j+1,k)))*C_AN*AXZ

      ! South face (i+1/2, j-1/2, k)
      D_FS = avg_h(avg_h(mu_g(i  ,j-1,k),mu_g(IC,J,K)),&
                   avg_h(mu_g(i+1,j-1,k),mu_g(i+1,j,k)))*C_AS*AXZ

      ! Top face (i+1/2, j, k+1/2)
      D_FT = avg_h(avg_h(mu_g(IC,J,K),mu_g(i,j,k+1)),&
                   avg_h(mu_g(i+1,j,k),mu_g(i+1,j,k+1)))*C_AT*AXY

      ! Bottom face (i+1/2, j, k-1/2)
      D_FB = avg_h(avg_h(mu_g(i  ,j,k-1),mu_g(IC ,j,k)),&
                   avg_h(mu_g(i+1,j,k-1),mu_g(i+1,j,k)))*C_AB*AXY

    CONTAINS

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
      real(c_real) :: flux_e, flux_n, flux_t
! Diffusion parameter
      real(c_real) :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

!---------------------------------------------------------------------//

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = slo(1),hi(1)

! Calculate convection-diffusion fluxes through each of the faces
               flux_e = HALF * (flux_gE(i,j,k) + flux_gE(i+1,j,k))
               flux_n = HALF * (flux_gN(i,j,k) + flux_gN(i+1,j,k))
               flux_t = HALF * (flux_gT(i,j,k) + flux_gT(i+1,j,k))

               CALL GET_UCELL_GDIFF_TERMS(&
                  slo, shi, &
                  d_fe, d_fw, d_fn, d_fs, &
                  d_ft, d_fb, mu_g, i, j, k, flag, &
                  dx, dy, dz)

! East face (i+1, j, k)
               if (flux_e >= zero) then
                  a_u_g(i,  j,k,e) = d_fe
                  a_u_g(i+1,j,k,w) = d_fe + flux_e
               else
                  a_u_g(i,  j,k,e) = d_fe - flux_e
                  a_u_g(i+1,j,k,w) = d_fe
               endif

! North face (i+1/2, j+1/2, k)
               if (flux_n >= zero) then
                  a_u_g(i,j,  k,n) = d_fn
                  a_u_g(i,j+1,k,s) = d_fn + flux_n
               else
                  a_u_g(i,j,  k,n) = d_fn - flux_n
                  a_u_g(i,j+1,k,s) = d_fn
               endif

! Top face (i+1/2, j, k+1/2)
               if (flux_t >= zero) then
                  a_u_g(i,j,k,  t) = d_ft
                  a_u_g(i,j,k+1,b) = d_ft + flux_t
               else
                  a_u_g(i,j,k,  t) = d_ft - flux_t
                  a_u_g(i,j,k+1,b) = d_ft
               endif

! South face (i+1/2, j-1/2, k)
               if(j==lo(2)) a_u_g(i,j,k,s) = d_fs

! Bottom face (i+1/2, j, k-1/2)
               if(k==lo(3)) a_u_g(i,j,k,b) = d_fb

            enddo
         enddo
      enddo

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

      use functions, only: avg
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
      integer :: i,j,k

      ! indicator for shear
      integer :: incr

      ! Diffusion parameter
      real(c_real) :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb

      ! Face mass flux
      real(c_real) :: flux_e, flux_n, flux_t

      real(c_real), allocatable :: u(:,:,:), v(:,:,:), ww(:,:,:)
      real(c_real), allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)

      allocate(  u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)) )
      allocate(  v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)) )
      allocate( ww(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)) )

      allocate(xsi_e(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)) )
      allocate(xsi_n(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)) )
      allocate(xsi_t(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)) )

      !  Calculate the components of velocity on the east, north, 
      !  and top face of a u-momentum cell.
      DO k = slo(3),shi(3)
        DO j = slo(2),shi(2)
          DO i = slo(1),hi(1)
             u(I,J,K) = avg(u_g(i,j,k), U_G(i+1,j,k))
             v(I,J,K) = avg(v_g(i,j,k), V_G(i+1,j,k))
            ww(I,J,K) = avg(w_g(i,j,k), W_G(i+1,j,k))
          ENDDO
        ENDDO
      ENDDO

      ! shear indicator:
      incr=1
      call calc_xsi (discretize(3), slo, shi, hi, &
         u_g, u, v, ww, xsi_e, xsi_n, xsi_t, &
         dt, dx, dy, dz)


      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = slo(1),hi(1)

! Calculate convection-diffusion fluxes through each of the faces
               flux_e = HALF * (flux_gE(i  ,j,k) + flux_gE(i+1,j,k))
               flux_n = HALF * (flux_gN(i,j  ,k) + flux_gN(i+1,j  ,k))
               flux_t = HALF * (flux_gT(i,j,k  ) + flux_gT(i+1,j,k  ))

               CALL GET_UCELL_GDIFF_TERMS(&
                  slo, shi, &
                  d_fe, d_fw, d_fn, d_fs, &
                  d_ft, d_fb, mu_g, i, j, k, flag, &
                  dx, dy, dz)

               ! East face (i+1, j, k)
               a_u_g(i,  j,k,e) = d_fe - flux_e*(      xsi_e(i,j,k))
               a_u_g(i+1,j,k,w) = d_fe + flux_e*(one - xsi_e(i,j,k))

               ! North face (i+1/2, j+1/2, k)
               a_u_g(i,j,  k,n) = d_fn - flux_n*(      xsi_n(i,j,k))
               a_u_g(i,j+1,k,s) = d_fn + flux_n*(one - xsi_n(i,j,k))

               ! Top face (i+1/2, j, k+1/2)
               a_u_g(i,j,k,  t) = d_ft - flux_t*(      xsi_t(i,j,k))
               a_u_g(i,j,k+1,b) = d_ft + flux_t*(one - xsi_t(i,j,k))

               ! South face (i+1/2, j-1/2, k)
               if(j==lo(2)) a_u_g(i,j,k,s) = d_fs

               ! Bottom face (i+1/2, j, k-1/2)
               if(k==lo(3)) a_u_g(i,j,k,b) = d_fb

            enddo
         enddo
      enddo

      deallocate( u, v, ww )
      deallocate( xsi_e, xsi_n, xsi_t)

      return
   end subroutine store_a_u_g1
end module u_g_conv_dif
