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
      A_m, mu_g, u_g, v_g, w_g, fluxX, fluxY, fluxZ,&
      dt, dx, dy, dz, domlo, domhi)


! Modules
!---------------------------------------------------------------------//
      use run, only: discretize

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: fluxZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

!---------------------------------------------------------------------//

      ! DO NOT use DEFERRED CORRECTION TO SOLVE V_G
      IF (discretize(4) == 0) THEN               ! 0 & 1 => FOUR
         call store_a_v_g0(&
              slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
              A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)
      else
         call store_a_v_g1(&
              slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
              A_m, mu_g, u_g, v_g, w_g, fluxX, fluxY, fluxZ, &
              dt, dx, dy, dz, domlo, domhi)
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
         A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)

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

      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

!     Local variables
!---------------------------------------------------------------------//

      integer :: i, j, k

      ! Face mass flux
      real(c_real) :: lflux

      ! Diffusion parameter
      real(c_real) :: d_f

      real(c_real) :: ayz_x, axz_y, axy_z

      ayz_x = dy*dz / dx
      axz_y = dx*dz / dy
      axy_z = dx*dy / dz

!---------------------------------------------------------------------//


      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1)-1,ahi(1)

               ! Calculate convection-diffusion fluxes through each of the faces
               lflux = HALF * (fluxX(i,j,k) + fluxX(i,j+1,k))

               d_f = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i+1,j  ,k)),&
                            avg_h(mu_g(i,j+1,k),mu_g(i+1,j+1,k))) * ayz_x

               ! East face (i+1, j, k)
               if (lflux >= zero) then
                  if (i.ge.alo(1)) A_m(i,  j,k,e) = d_f
                  if (i.lt.ahi(1)) A_m(i+1,j,k,w) = d_f + lflux
               else
                  if (i.ge.alo(1)) A_m(i,  j,k,e) = d_f - lflux
                  if (i.lt.ahi(1)) A_m(i+1,j,k,w) = d_f
               endif

            enddo
         enddo
      enddo

      do k = alo(3),ahi(3)
         do j = alo(2)-1,ahi(2)
            do i = alo(1),ahi(1)

               lflux = HALF * (fluxY(i,j,k) + fluxY(i,j+1,k))
               d_f = mu_g(i,j+1,k) * axz_y

               ! North face (i+1/2, j+1/2, k)
               if (lflux >= zero) then
                  if (j.ge.alo(2)) A_m(i,j,  k,n) = d_f
                  if (j.lt.ahi(2)) A_m(i,j+1,k,s) = d_f + lflux
               else
                  if (j.ge.alo(2)) A_m(i,j,  k,n) = d_f - lflux
                  if (j.lt.ahi(2)) A_m(i,j+1,k,s) = d_f
               endif

            enddo
         enddo
      enddo

      do k = alo(3)-1,ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               lflux = HALF * (fluxZ(i,j,k) + fluxZ(i,j+1,k))

               d_f = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i,j  ,k+1)),&
                            avg_h(mu_g(i,j+1,k),mu_g(i,j+1,k+1))) * axy_z

               if (lflux >= zero) then
                  if (k.ge.alo(3)) A_m(i,j,k,  t) = d_f
                  if (k.lt.ahi(3)) A_m(i,j,k+1,b) = d_f + lflux
               else
                  if (k.ge.alo(3)) A_m(i,j,k,  t) = d_f - lflux
                  if (k.lt.ahi(3)) A_m(i,j,k+1,b) = d_f
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
         A_m, mu_g, u_g, v_g, w_g, fluxX, fluxY, fluxZ, &
         dt, dx, dy, dz, domlo, domhi)

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
      integer     , intent(in   ) :: domlo(3),domhi(3)
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

      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Local variables
!---------------------------------------------------------------------//
      integer :: xlo(3), xhi(3)
      integer :: vello(3),velhi(3)
      integer :: i, j, k

      ! Diffusion parameter
      real(c_real) :: d_f

      ! Face mass flux
      real(c_real) :: lflux

      real(c_real), allocatable :: vel(:,:,:)
      real(c_real), allocatable :: xsi_(:,:,:)

      real(c_real) :: ayz_x, axz_y, axy_z

      ayz_x = dy*dz / dx
      axz_y = dx*dz / dy
      axy_z = dx*dy / dz

      vello(1) = slo(1)-2
      vello(2) = slo(2)-2
      vello(3) = slo(3)-2
      velhi(1) = shi(1)+2
      velhi(2) = shi(2)+2
      velhi(3) = shi(3)+2

      xhi(1) = ahi(1)
      xhi(2) = ahi(2)
      xhi(3) = ahi(3)

      allocate( vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3)))
!---------------------------------------------------------------------//


      xlo(1) = alo(1)-1
      xlo(2) = alo(2)
      xlo(3) = alo(3)

      vel(:,:,:) = 0.d0
      do k = vlo(3),vhi(3)
        do j = vlo(2)+1,vhi(2)-1
          do i = vlo(1),vhi(1)
             vel(i,j,k) = avg(u_g(i,j,k), u_g(i,j+1,k))
          end do
        end do
      end do

      allocate( xsi_(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_e (discretize(4), v_g, vlo, vhi, vel, vello, velhi, &
         xsi_, xlo, xhi, dt, dx, dy, dz, domlo, domhi)

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1)-1,ahi(1)

               lflux = half * (fluxX(i,j,k) + fluxX(i,j+1,k))

               d_f = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i+1,j  ,k)),&
                            avg_h(mu_g(i,j+1,k),mu_g(i+1,j+1,k))) * ayz_x

               if (i.ge.alo(1)) A_m(i,  j,k,e) = d_f - lflux*(xsi_(i,j,k))
               if (i.lt.ahi(1)) A_m(i+1,j,k,w) = d_f + lflux*(one - xsi_(i,j,k))

            enddo
         enddo
      enddo
      deallocate(xsi_)

!---------------------------------------------------------------------//

      xlo(1) = alo(1)
      xlo(2) = alo(2)-1
      xlo(3) = alo(3)

      vel(:,:,:) = 0.d0
      do k = vlo(3),vhi(3)
        do j = vlo(2)+1,vhi(2)-1
          do i = vlo(1),vhi(1)
             vel(i,j,k) = avg(v_g(i,j,k), v_g(i,j+1,k))
          end do
        end do
      end do

      allocate( xsi_(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_n (discretize(4), v_g, vlo, vhi, vel, vello, velhi, &
         xsi_, xlo, xhi, dt, dx, dy, dz, domlo, domhi)

     do k = alo(3),ahi(3)
         do j = alo(2)-1,ahi(2)
            do i = alo(1),ahi(1)

               lflux = half * (fluxY(i,j,k) + fluxY(i,j+1,k))

               d_f = mu_g(i,j+1,k) * axz_y

               if (j.ge.alo(2)) A_m(i,j  ,k,n) = d_f - lflux*(      xsi_(i,j,k))
               if (j.lt.ahi(2)) A_m(i,j+1,k,s) = d_f + lflux*(one - xsi_(i,j,k))

            enddo
         enddo
      enddo
      deallocate(xsi_)

!---------------------------------------------------------------------//

      xlo(1) = alo(1)
      xlo(2) = alo(2)
      xlo(3) = alo(3)-1

      vel(:,:,:) = 0.d0
      do k = vlo(3),vhi(3)
        do j = vlo(2)+1,vhi(2)-1
          do i = vlo(1),vhi(1)
             vel(i,j,k) = avg(w_g(i,j,k), w_g(i,j+1,k))
          end do
        end do
      end do

      allocate( xsi_(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_t (discretize(4), v_g, vlo, vhi, vel, vello, velhi, &
         xsi_, xlo, xhi, dt, dx, dy, dz, domlo, domhi)

     do k = alo(3)-1,ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               lflux = half * (fluxZ(i,j,k) + fluxZ(i,j+1,k))

               d_f = avg_h(avg_h(mu_g(i,j  ,k),mu_g(i,j  ,k+1)),&
                            avg_h(mu_g(i,j+1,k),mu_g(i,j+1,k+1))) * axy_z

               if (k.ge.alo(3)) A_m(i,j,k,  t) = d_f - lflux*(xsi_(i,j,k))
               if (k.lt.ahi(3)) A_m(i,j,k+1,b) = d_f + lflux*(one - xsi_(i,j,k))

            end do
         end do
      end do
      deallocate(xsi_)

!---------------------------------------------------------------------//

      deallocate(vel)

      end subroutine store_a_v_g1

end module v_g_conv_dif
