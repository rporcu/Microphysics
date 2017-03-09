module v_g_conv_dif

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int
   use param1        , only: half, one, zero

   implicit none

   private
   public :: conv_dif_v_g

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: conv_div_v_g(A_m)                                      !
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  !
!  The off-diagonal coefficients calculated here must be positive. The !
!  center coefficient and the source vector are negative;              !
!  See source_v_g                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine conv_dif_v_g(&
      slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
      A_m, mu_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt,&
      dt, dx, dy, dz)


! Modules
!---------------------------------------------------------------------//
      use run, only: discretize

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

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

!---------------------------------------------------------------------//

      ! DO NOT use DEFERRED CORRECTION TO SOLVE V_G
      IF (discretize(4) == 0) THEN               ! 0 & 1 => FOUR
         call store_a_v_g0(&
              slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
              A_m, mu_g, flux_ge, flux_gn, flux_gt, dx, dy, dz)
      else
         call store_a_v_g1(&
              slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
              A_m, mu_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt, &
              dt, dx, dy, dz)
      endIF

      end subroutine conv_dif_v_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: store_a_v_g0                                            !
!  Purpose: Determine convection diffusion terms for V_g momentum eqs. !
!  The off-diagonal coefficients calculated here must be positive.     !
!  The center coefficient and the source vector are negative. See      !
!  source_v_g.                                                         !
!                                                                      !
!  Implement FOUP discretization                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine store_a_v_g0(&
         slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
         A_m, mu_g, flux_ge, flux_gn, flux_gt, dx, dy, dz)

      use functions, only: avg_h
      use matrix, only: e, w, n, s, t, b

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

!     Local variables
!---------------------------------------------------------------------//

      integer :: i, j, k

      ! Face mass flux
      real(c_real) :: flux_e, flux_n, flux_t

      ! Diffusion parameter
      real(c_real) :: d_fe, d_fn, d_ft

      real(c_real) :: c_ae,c_aw,c_an,c_as,c_at,c_ab

      c_ae = dy*dz / dx
      c_aw = dy*dz / dx
      c_an = dx*dz / dy
      c_as = dx*dz / dy
      c_at = dx*dy / dz
      c_ab = dx*dy / dz

!---------------------------------------------------------------------//


      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1)-1,ahi(1)

               ! Calculate convection-diffusion fluxes through each of the faces
               flux_e = HALF * (flux_gE(i,j,k) + flux_gE(i,j+1,k))

               d_fe = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i+1,j  ,k)),&
                            avg_h(mu_g(i,j+1,k),mu_g(i+1,j+1,k))) * c_ae

               ! East face (i+1, j, k)
               if (flux_e >= zero) then
                  if (i.ge.alo(1)) A_m(i,  j,k,e) = d_fe
                  if (i.lt.ahi(1)) A_m(i+1,j,k,w) = d_fe + flux_e
               else
                  if (i.ge.alo(1)) A_m(i,  j,k,e) = d_fe - flux_e
                  if (i.lt.ahi(1)) A_m(i+1,j,k,w) = d_fe
               endif

            enddo
         enddo
      enddo

      do k = alo(3),ahi(3)
         do j = alo(2)-1,ahi(2)
            do i = alo(1),ahi(1)

               flux_n = HALF * (flux_gN(i,j,k) + flux_gN(i,j+1,k))
               d_fn = mu_g(i,j+1,k) * c_an

               ! North face (i+1/2, j+1/2, k)
               if (flux_n >= zero) then
                  if (j.ge.alo(2)) A_m(i,j,  k,n) = d_fn
                  if (j.lt.ahi(2)) A_m(i,j+1,k,s) = d_fn + flux_n
               else
                  if (j.ge.alo(2)) A_m(i,j,  k,n) = d_fn - flux_n
                  if (j.lt.ahi(2)) A_m(i,j+1,k,s) = d_fn
               endif

            enddo
         enddo
      enddo

      do k = alo(3)-1,ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               flux_t = HALF * (flux_gT(i,j,k) + flux_gT(i,j+1,k))

               d_ft = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i,j  ,k+1)),&
                            avg_h(mu_g(i,j+1,k),mu_g(i,j+1,k+1))) * c_at

               if (flux_t >= zero) then
                  if (k.ge.alo(3)) A_m(i,j,k,  t) = d_ft
                  if (k.lt.ahi(3)) A_m(i,j,k+1,b) = d_ft + flux_t
               else
                  if (k.ge.alo(3)) A_m(i,j,k,  t) = d_ft - flux_t
                  if (k.lt.ahi(3)) A_m(i,j,k+1,b) = d_ft
               endif

            enddo
         enddo
      enddo

      end subroutine store_a_v_g0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: store_a_v_g1                                            !
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  !
!  The off-diagonal coefficients calculated here must be positive.     !
!  The center coefficient and the source vector are negative.          !
!  Implements higher order discretization.                             !
!  See source_v_g                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine store_a_v_g1(&
         slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
         A_m, mu_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt, &
         dt, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//

      use functions, only: avg, avg_h
      use matrix, only: e, w, n, s, t, b
      use run, only: discretize

      use xsi, only: calc_xsi_e, calc_xsi_n, calc_xsi_t

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx, dy, dz, dt

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent( in) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent( in) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent( in) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Local variables
!---------------------------------------------------------------------//
      integer :: xlo(3), xhi(3)
      integer :: vello(3),velhi(3)
      integer :: i, j, k

      ! Diffusion parameter
      real(c_real) :: d_fe, d_fn, d_ft

      ! Face mass flux
      real(c_real) :: flux_e, flux_n, flux_t

      real(c_real), allocatable :: u(:,:,:), v(:,:,:), ww(:,:,:)
      real(c_real), allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)

      real(c_real) :: c_ae,c_aw,c_an,c_as,c_at,c_ab

      allocate(  u(slo(1)-2:shi(1)+2,slo(2)-2:shi(2)+2,slo(3)-2:shi(3)+2) )
      allocate(  v(slo(1)-2:shi(1)+2,slo(2)-2:shi(2)+2,slo(3)-2:shi(3)+2) )
      allocate( ww(slo(1)-2:shi(1)+2,slo(2)-2:shi(2)+2,slo(3)-2:shi(3)+2) )

      c_ae = dy*dz / dx
      c_aw = dy*dz / dx
      c_an = dx*dz / dy
      c_as = dx*dz / dy
      c_at = dx*dy / dz
      c_ab = dx*dy / dz

      vello(1) = slo(1)-2
      vello(2) = slo(2)-2
      vello(3) = slo(3)-2
      velhi(1) = shi(1)+2
      velhi(2) = shi(2)+2
      velhi(3) = shi(3)+2

      xhi(1) = ahi(1)
      xhi(2) = ahi(2)
      xhi(3) = ahi(3)

!---------------------------------------------------------------------//

       u(:,:,:) = 0.d0
       v(:,:,:) = 0.d0
      ww(:,:,:) = 0.d0

      do k = vlo(3),vhi(3)
        do j = vlo(2)+1,vhi(2)-1
          do i = vlo(1),vhi(1)
             u(i,j,k) = avg(u_g(i,j,k), u_g(i,j+1,k))
          end do
        end do
      end do

      xlo(1) = alo(1)-1
      xlo(2) = alo(2)
      xlo(3) = alo(3)

      allocate( xsi_e(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_e (discretize(4), v_g, vlo, vhi, u, vello, velhi, xsi_e, xlo, xhi, &
                       dt, dx, dy, dz)

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1)-1,ahi(1)

               flux_e = half * (flux_ge(i,j,k) + flux_ge(i,j+1,k))

               d_fe = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i+1,j  ,k)),&
                            avg_h(mu_g(i,j+1,k),mu_g(i+1,j+1,k))) * c_ae

               if (i.ge.alo(1)) A_m(i,  j,k,e) = d_fe - flux_e*(xsi_e(i,j,k))
               if (i.lt.ahi(1)) A_m(i+1,j,k,w) = d_fe + flux_e*(one - xsi_e(i,j,k))

            enddo
         enddo
      enddo

!---------------------------------------------------------------------//

       u(:,:,:) = 0.d0
       v(:,:,:) = 0.d0
      ww(:,:,:) = 0.d0

      do k = vlo(3),vhi(3)
        do j = vlo(2)+1,vhi(2)-1
          do i = vlo(1),vhi(1)
             v(i,j,k) = avg(v_g(i,j,k), v_g(i,j+1,k))
          end do
        end do
      end do

      xlo(1) = alo(1)
      xlo(2) = alo(2)-1
      xlo(3) = alo(3)
      allocate( xsi_n(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_n (discretize(4), v_g, vlo, vhi, v, vello, velhi, xsi_n, xlo, xhi, &
                       dt, dx, dy, dz)

     do k = alo(3),ahi(3)
         do j = alo(2)-1,ahi(2)
            do i = alo(1),ahi(1)

               flux_n = half * (flux_gn(i,j,k) + flux_gn(i,j+1,k))

               d_fn = mu_g(i,j+1,k) * c_an

               if (j.ge.alo(2)) A_m(i,j  ,k,n) = d_fn - flux_n*(      xsi_n(i,j,k))
               if (j.lt.ahi(2)) A_m(i,j+1,k,s) = d_fn + flux_n*(one - xsi_n(i,j,k))

            enddo
         enddo
      enddo

!---------------------------------------------------------------------//

       u(:,:,:) = 0.d0
       v(:,:,:) = 0.d0
      ww(:,:,:) = 0.d0

      do k = vlo(3),vhi(3)
        do j = vlo(2)+1,vhi(2)-1
          do i = vlo(1),vhi(1)
             ww(i,j,k) = avg(w_g(i,j,k), w_g(i,j+1,k))
          end do
        end do
      end do

      xlo(1) = alo(1)
      xlo(2) = alo(2)
      xlo(3) = alo(3)-1

      allocate( xsi_t(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_t (discretize(4), v_g, vlo, vhi, ww, vello, velhi, xsi_t, xlo, xhi, &
                       dt, dx, dy, dz)

     do k = alo(3)-1,ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               flux_t = half * (flux_gt(i,j,k) + flux_gt(i,j+1,k))

               d_ft = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i,j  ,k+1)),&
                            avg_h(mu_g(i,j+1,k),mu_g(i,j+1,k+1))) * c_at

               if (k.ge.alo(3)) A_m(i,j,k,  t) = d_ft - flux_t*(xsi_t(i,j,k))
               if (k.lt.ahi(3)) A_m(i,j,k+1,b) = d_ft + flux_t*(one - xsi_t(i,j,k))

            end do
         end do
      end do

      deallocate( u, v, ww )
      deallocate( xsi_e, xsi_n, xsi_t)

      end subroutine store_a_v_g1

end module v_g_conv_dif
