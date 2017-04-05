module u_g_conv_dif

   use amrex_fort_module, only : c_real => amrex_real
   use param1        , only: half, one, zero

   implicit none

   private
   public :: conv_dif_u_g

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: conv_dif_u_g                                            !
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  !
!  The off-diagonal coefficients calculated here must be positive. The !
!  center coefficient and the source vector are negative;              !
!  See source_u_g                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine conv_dif_u_g(&
      slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
      A_m, mu_g, u_g, v_g, w_g, fluxX, fluxY, fluxZ, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use run, only: discretize

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)

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

      real(c_real), intent(in   ) :: dx, dy, dz
!---------------------------------------------------------------------//

      if (discretize(3) == 0) then

         call store_a_u_g0(&
            slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
            A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)
      else
         call store_a_u_g1(&
            slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
            A_m, mu_g, u_g, v_g, w_g, fluxX, fluxY, fluxZ, &
            dx, dy, dz)
      end if

      end subroutine conv_dif_u_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: store_a_u_g0                                            !
!  Purpose: Determine convection diffusion terms for U_g momentum eqs. !
!  The off-diagonal coefficients calculated here must be positive.     !
!  The center coefficient and the source vector are negative. See      !
!  source_u_g.                                                         !
!                                                                      !
!  Implement FOUP discretization                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine store_a_u_g0(&
      slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
      A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)

      use functions, only: avg_h
      use matrix   , only: e, w, n, s, t, b

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx, dy, dz

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

! Local variables
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

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1)-1,ahi(1)

               ! Calculate convection-diffusion fluxes through each of the faces
               lflux = HALF * (fluxX(i  ,j,k) + fluxX(i+1,j,k))

               d_f = mu_g(i,j,k) * ayz_x

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

               lflux = HALF * (fluxY(i-1,j+1,k) + fluxY(i  ,j+1,k))

               d_f = avg_h(avg_h(mu_g(i-1,j,k),mu_g(i-1,j+1,k)),&
                           avg_h(mu_g(i  ,j,k),mu_g(i  ,j+1,k))) * axz_y

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

               lflux = HALF * (fluxZ(i-1,j,k+1) + fluxZ(i  ,j,k+1))

               d_f = avg_h(avg_h(mu_g(i-1,j,k),mu_g(i-1,j,k+1)),&
                           avg_h(mu_g(i  ,j,k),mu_g(i  ,j,k+1))) * axy_z

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

      return
   end subroutine store_a_u_g0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: store_a_u_g1                                            !
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  !
!  The off-diagonal coefficients calculated here must be positive.     !
!  The center coefficient and the source vector are negative.          !
!  Implements higher order discretization.                             !
!  See source_u_g                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine store_a_u_g1(&
         slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
         A_m, mu_g, u_g, v_g, w_g, fluxX, fluxY, fluxZ, &
         dx, dy, dz)

      use functions, only: avg, avg_h
      use matrix   , only: e, w, n, s, t, b

      use xsi, only: calc_xsi_x, calc_xsi_y, calc_xsi_z

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx, dy, dz


      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Local variables
!---------------------------------------------------------------------//
      integer :: xlo(3),xhi(3)
      integer :: vello(3),velhi(3)
      integer :: i,j,k

      ! Diffusion parameter
      real(c_real) :: d_f,d_f_s,d_f_n,d_f_b,d_f_t

      ! Face mass flux
      real(c_real) :: lflux,lflux_s,lflux_n,lflux_b,lflux_t

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

      allocate( vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3)))

!---------------------------------------------------------------------//

      xlo(1) = alo(1)-1
      xlo(2) = alo(2)
      xlo(3) = alo(3)

      xhi(1) = ahi(1)
      xhi(2) = ahi(2)
      xhi(3) = ahi(3)

      vel(:,:,:) = 0.d0
      do k = ulo(3),uhi(3)
        do j = ulo(2),uhi(2)
          do i = ulo(1),uhi(1)-1
             vel(i,j,k) = avg(u_g(i,j,k), u_g(i+1,j,k))
          end do
        end do
      end do

      ! NOTES:   u_g  lives on x-faces    :   (lo(1): hi(1)+1, lo(2):hi(2)  , lo(3):hi(3))
      !          A_m  lives on x-faces    :   (lo(1): hi(1)+1, lo(2):hi(2)  , lo(3):hi(3))
      !          vel  lives on cell ctrs  :   (lo(1): hi(1)  , lo(2):hi(2)  , lo(3):hi(3))
      !          xsi  lives on cell ctrs  :   (lo(1): hi(1)  , lo(2):hi(2)  , lo(3):hi(3))
      !               (so xsi(i-1,j,k) and xsi(i,j,k) contribute to A_m(i,j,k))

      allocate(xsi_(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_x (u_g, ulo, uhi, vel, vello, velhi, xsi_, xlo, xhi, .false.)

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1)-1,ahi(1)

               lflux = half * (fluxX(i  ,j,k) + fluxX(i+1,j,k))
               d_f = ayz_x*mu_g(i,j,k)

               if (i.ge.alo(1)) A_m(i,  j,k,e) = d_f - lflux*(      xsi_(i,j,k))
               if (i.lt.ahi(1)) A_m(i+1,j,k,w) = d_f + lflux*(one - xsi_(i,j,k))

            enddo
         enddo
      enddo
      deallocate(xsi_)

!---------------------------------------------------------------------//

      xhi(1) = ahi(1)
      xhi(2) = ahi(2)+1
      xhi(3) = ahi(3)

      vel(:,:,:) = 0.d0
      do k = ulo(3),uhi(3)
        do j = ulo(2),uhi(2)
          do i = ulo(1)+1,uhi(1)-1
             vel(i,j,k) = avg(v_g(i-1,j,k), v_g(i,j,k))
          end do
        end do
      end do

      allocate(xsi_(alo(1):xhi(1),alo(2):xhi(2),alo(3):xhi(3)) )
      call calc_xsi_y (u_g, ulo, uhi, vel, vello, velhi, xsi_, alo, xhi, .true.)
  
      ! NOTES:   u_g  lives on x-faces   :   (lo(1): hi(1)+1, lo(2):hi(2)  , lo(3):hi(3))
      !          A_m  lives on x-faces   :   (lo(1): hi(1)+1, lo(2):hi(2)  , lo(3):hi(3))
      !          vel  lives on x-y edges :   (lo(1): hi(1)+1, lo(2):hi(2)+1, lo(3):hi(3))
      !          xsi  lives on x-y edges :   (lo(1): hi(1)+1, lo(2):hi(2)+1, lo(3):hi(3))

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               lflux_n = HALF * (fluxY(i-1,j+1,k) + fluxY(i  ,j+1,k))
               lflux_s = HALF * (fluxY(i-1,j  ,k) + fluxY(i  ,j  ,k))

               d_f_n = avg_h(avg_h(mu_g(i-1,j,k),mu_g(i-1,j+1,k)),&
                             avg_h(mu_g(i  ,j,k),mu_g(i  ,j+1,k))) * axz_y

               d_f_s = avg_h(avg_h(mu_g(i-1,j-1,k),mu_g(i-1,j,k)),&
                             avg_h(mu_g(i  ,j-1,k),mu_g(i  ,j,k))) * axz_y

               A_m(i,j,k,n) = d_f_n - lflux_n*(      xsi_(i,j+1,k))
               A_m(i,j,k,s) = d_f_s + lflux_s*(one - xsi_(i,j  ,k))

            enddo
         enddo
      enddo
      deallocate(xsi_)

!---------------------------------------------------------------------//

      xhi(1) = ahi(1)
      xhi(2) = ahi(2)
      xhi(3) = ahi(3)+1

      vel(:,:,:) = 0.d0
      do k = ulo(3),uhi(3)
        do j = ulo(2),uhi(2)
          do i = ulo(1)+1,uhi(1)-1
            vel(i,j,k) = avg(w_g(i-1,j,k), w_g(i,j,k))
          end do
        end do
      end do
  
      ! NOTES:   u_g  lives on x-faces   :   (lo(1): hi(1)+1, lo(2):hi(2), lo(3):hi(3)  )
      !          A_m  lives on x-faces   :   (lo(1): hi(1)+1, lo(2):hi(2), lo(3):hi(3)  )
      !          vel  lives on x-z edges :   (lo(1): hi(1)+1, lo(2):hi(2), lo(3):hi(3)+1)
      !          xsi  lives on x-z edges :   (lo(1): hi(1)+1, lo(2):hi(2), lo(3):hi(3)+1)

      allocate(xsi_(alo(1):xhi(1),alo(2):xhi(2),alo(3):xhi(3)) )
      call calc_xsi_z (u_g, ulo, uhi, vel, vello, velhi, xsi_, alo, xhi, .true.)

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               lflux_t = HALF * (fluxZ(i-1,j,k+1) + fluxZ(i  ,j,k+1))
               lflux_b = HALF * (fluxZ(i-1,j,k  ) + fluxZ(i  ,j,k  ))

               d_f_t = avg_h(avg_h(mu_g(i-1,j,k),mu_g(i-1,j,k+1)),&
                             avg_h(mu_g(i  ,j,k),mu_g(i  ,j,k+1))) * axy_z

               d_f_b = avg_h(avg_h(mu_g(i-1,j,k-1),mu_g(i-1,j,k)),&
                             avg_h(mu_g(i  ,j,k-1),mu_g(i  ,j,k))) * axy_z

               A_m(i,j,k,t) = d_f_t - lflux_t*(      xsi_(i,j,k+1))
               A_m(i,j,k,b) = d_f_b + lflux_b*(one - xsi_(i,j,k  ))

            enddo
         enddo
      enddo
      deallocate(xsi_)

!---------------------------------------------------------------------//

      deallocate(vel )

   end subroutine store_a_u_g1
end module u_g_conv_dif
