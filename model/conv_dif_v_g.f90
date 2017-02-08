module v_g_conv_dif

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   use param1        , only: half, one, zero

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
   subroutine conv_dif_v_g(&
      slo, shi, lo, hi, ulo, uhi, vlo, vhi, wlo, whi, &
      A_m, mu_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt,&
      dt, dx, dy, dz)


! Modules
!---------------------------------------------------------------------//
      use run, only: discretize

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: dt, dx, dy, dz
!---------------------------------------------------------------------//

! DO NOT use DEFERRED CORRECTION TO SOLVE V_G
      IF (DISCRETIZE(4) == 0) THEN               ! 0 & 1 => FOUP
         CALL STORE_A_V_G0(&
              slo,shi,lo,hi,A_M,mu_g,flux_ge,flux_gn,flux_gt,dx,dy,dz)
      ELSE
         CALL STORE_A_V_G1(&
              slo,shi,lo,hi,A_M,mu_g,u_g,v_g,w_g,flux_ge,flux_gn,flux_gt,dt,dx,dy,dz)
      endIF

      end subroutine conv_dif_v_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  the faces of a v-momentum cell. Note the fluxes are calculated at   C
!  all faces.                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine get_vcell_gdiff_terms(&
         slo, shi, &
         D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, mu_g, I, J, K, &
         dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use functions, only: avg, avg_h
      use param1, only: zero
      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)

      ! diffusion through faces of given v-momentum cell
      real(c_real), intent(  out) :: d_fe, d_fw
      real(c_real), intent(  out) :: d_fn, d_fs
      real(c_real), intent(  out) :: d_ft, d_fb

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: dx, dy, dz

      integer, intent(in) :: i, j, k

! Local variables
      real(c_real) :: odx, ody, odz
      real(c_real) :: axy, axz, ayz
!---------------------------------------------------------------------//
      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      ! East face (i+1/2, j+1/2, k)
      D_Fe = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i+1,j  ,k)),&
                   avg_h(mu_g(i,j+1,k),mu_g(i+1,j+1,k)))*ODX*AYZ

      ! West face (i-1/2, j+1/2, k)
      D_Fw = avg_h(avg_h(mu_g(i-1,j  ,k),mu_g(i,j  ,k)),&
                   avg_h(mu_g(i-1,j+1,k),mu_g(i,j+1,k)))*ODX*AYZ

      ! North face (i, j+1, k)
      D_Fn = mu_g(i,j+1,k)*ODY*AXZ

      ! South face (i, j, k)
      D_Fs = mu_g(i,j  ,k)*ODY*AXZ

      D_FT = ZERO
      D_FB = ZERO

      ! Top face (i, j+1/2, k+1/2)
      D_Ft = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i,j  ,k+1)),&
                   avg_h(mu_g(i,j+1,k),mu_g(i,j+1,k+1)))*ODZ*AXY

      ! Bottom face (i, j+1/2, k-1/2)
      D_Fb = avg_h(avg_h(mu_g(i,j  ,k+1),mu_g(i,j  ,k)),&
                   avg_h(mu_g(i,j+1,k-1),mu_g(i,j+1,k)))*ODZ*AXY

    end subroutine get_vcell_gdiff_terms

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
      subroutine store_A_V_G0(&
         slo, shi, lo, hi, &
         A_V_G,mu_g,flux_ge,flux_gn,flux_gt,dx,dy,dz)

! Modules
!---------------------------------------------------------------------//
      use param1, only: zero
      use matrix, only: e, w, n, s, t, b

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_V_g
      real(c_real), intent(inout) :: A_V_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: dx,dy,dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      integer :: I, J, K
! Face mass flux
      real(c_real) :: flux_e, flux_n
      real(c_real) :: flux_t
! Diffusion parameter
      real(c_real) :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

!---------------------------------------------------------------------//

      do k = lo(3),hi(3)
         do j = slo(2),hi(2)
            do i = lo(1),hi(1)

! Calculate convection-diffusion fluxes through each of the faces
               flux_e = half * (flux_ge(i,j,k) + flux_ge(i,j+1,k))
               flux_n = half * (flux_gn(i,j,k) + flux_gn(i,j+1,k))
               flux_t = half * (flux_gt(i,j,k) + flux_gt(i,j+1,k))

               CALL GET_VCELL_GDIFF_TERMS(&
                  slo, shi, &
                  d_fe, d_fw, d_fn, d_fs, &
                  d_ft, d_fb, mu_g, i, j, k, dx, dy, dz)

! East face (i+1/2, j+1/2, k)
               if (flux_e >= zero) then
                  a_v_g(i,  j,k,e) = d_fe
                  a_v_g(i+1,j,k,w) = d_fe + flux_e
               else
                  a_v_g(i,  j,k,e) = d_fe - flux_e
                  a_v_g(i+1,j,k,w) = d_fe
               endif

! North face (i, j+1, k)
               if (flux_n >= zero) then
                  a_v_g(i,j,  k,n) = d_fn
                  a_v_g(i,j+1,k,s) = d_fn + flux_n
               else
                  a_v_g(i,j,  k,n) = d_fn - flux_n
                  a_v_g(i,j+1,k,s) = d_fn
               endif

! Top face (i, j+1/2, k+1/2)
               if (flux_t >= zero) then
                  a_v_g(i,j,k,  t) = d_ft
                  a_v_g(i,j,k+1,b) = d_ft + flux_t
               else
                  a_v_g(i,j,k,  t) = d_ft - flux_t
                  a_v_g(i,j,k+1,b) = d_ft
               endif

! West face (i-1/2, j+1/2, k)
               if(i==lo(1)) a_v_g(i,j,k,w) = d_fw

! Bottom face (i+1/2, j, k-1/2)
               if(k==lo(3)) a_v_g(i,j,k,b) = d_fb

            enddo
         enddo
      enddo

      end subroutine store_A_V_G0


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
      subroutine store_A_V_G1(&
         slo, shi, lo, hi, A_V_G, mu_g, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt,dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//

      use functions, only: avg, avg_h
      use matrix, only: e, w, n, s, t, b
      use run, only: discretize
      use xsi, only: calc_xsi

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_V_g
      real(c_real), intent(inout) :: A_V_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent( in) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent( in) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent( in) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: dt, dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      integer :: I, J, K
! indicator for shear
      integer :: incr
! Diffusion parameter
      real(c_real) :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      real(c_real) :: Flux_e, flux_n
      real(c_real) :: flux_t

! x, y, z directional velocity
      real(c_real), allocatable :: u(:,:,:), v(:,:,:), ww(:,:,:)
      real(c_real), allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)
!---------------------------------------------------------------------//

      allocate(  U(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate(  V(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate( WW(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )

      allocate(xsi_e(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate(xsi_n(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      allocate(xsi_t(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )

      !  Calculate the components of velocity on the east, north,
      !  and top face of a v-momentum cell
      do k = slo(3),shi(3)
        do j = slo(2),hi(2)
          do i = slo(1),shi(1)
             u(I,J,K) = avg(u_g(i,j,k),u_g(i,j+1,k))
             v(I,J,K) = avg(v_g(i,j,k),v_g(i,j+1,k))
            ww(I,J,K) = avg(w_g(i,j,k),w_g(i,j+1,k))
          end do
        end do
      end do

      ! shear indicator: y-momentum
      incr=2
      call calc_xsi (discretize(4), slo, shi, hi, &
                     v_g, u, v, ww, XSI_E, XSI_N, XSI_T, &
                     dt, dx, dy, dz)

      do k = lo(3),hi(3)
         do j = slo(2),hi(2)
            do i = lo(1),hi(1)

! Calculate convection-diffusion fluxes through each of the faces
               flux_e = half * (flux_ge(i,j,k) + flux_ge(i,j+1,k))
               flux_n = half * (flux_gn(i,j,k) + flux_gn(i,j+1,k))
               flux_t = half * (flux_gt(i,j,k) + flux_gt(i,j+1,k))

               CALL GET_VCELL_GDIFF_TERMS(&
                  slo, shi, &
                  d_fe, d_fw, d_fn, d_fs, &
                  d_ft, d_fb, mu_g, i, j, k, dx, dy, dz)

! East face (i+1/2, j+1/2, k)
               a_v_g(i,  j,k,e) = d_fe - flux_e*(xsi_e(i,j,k))
               a_v_g(i+1,j,k,w) = d_fe + flux_e*(one - xsi_e(i,j,k))

! North face (i, j+1, k)
               a_v_g(i,j  ,k,n) = d_fn - flux_n*(xsi_n(i,j,k))
               a_v_g(i,j+1,k,s) = d_fn + flux_n*(one - xsi_n(i,j,k))

! Top face (i, j+1/2, k+1/2)
               a_v_g(i,j,k,  t) = d_ft - flux_t*(xsi_t(i,j,k))
               a_v_g(i,j,k+1,b) = d_ft + flux_t*(one - xsi_t(i,j,k))


! West face (i-1/2, j+1/2, k)
               if(i==lo(1)) a_v_g(i,j,k,w) = d_fw

! Bottom face (i+1/2, j, k-1/2)
               if(k==lo(3)) a_v_g(i,j,k,b) = d_fb

            enddo
         enddo
      enddo

      deallocate( U, V, WW )
      deallocate( xsi_e, xsi_n, xsi_t)

      end subroutine store_A_V_G1
end module v_g_conv_dif
