module w_g_conv_dif

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int
   use param1        , only: half, one, zero

   implicit none

   private
   public :: conv_dif_w_g

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: conv_dif_w_g                                            !
!  Purpose: Determine convection diffusion terms for w_g momentum eqs  !
!  The off-diagonal coefficients calculated here must be positive. The !
!  center coefficient and the source vector are negative;              !
!  See source_w_g                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine conv_dif_w_g(&
      slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
      A_m, mu_g, u_g, v_g, w_g, fluxX, fluxY, fluxZ,&
      dx, dy, dz)


! Modules
!---------------------------------------------------------------------//
      use run, only: discretize

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real) :: A_m&
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

      if (discretize(5) == 0) then               ! 0 & 1 => FOUR
         call store_a_w_g0 (slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
                            A_m, mu_g, fluxX, fluxY, fluxZ, dx, dy, dz)
      else
         call store_a_w_g1 (slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
                            A_m, mu_g, u_g, v_g, w_g, fluxX, fluxY, fluxZ, &
                            dx, dy, dz)
      end if

      end subroutine conv_dif_w_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: store_a_w_g0                                            !
!  Purpose: Determine convection diffusion terms for W_g momentum eqs. !
!  The off-diagonal coefficients calculated here must be positive.     !
!  The center coefficient and the source vector are negative. See      !
!  source_w_g.                                                         !
!  Implement FOUP discretization                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine store_a_w_g0(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
                              A_m, mu_g, fluxX, fluxY, fluxZ, &
                              dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      use functions, only: avg_h
      use matrix   , only: e, w, n, s, t, b

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

      real(c_real) :: ayz_x,axz_y,axy_z

      ayz_x = dy*dz / dx
      axz_y = dx*dz / dy
      axy_z = dx*dy / dz

!---------------------------------------------------------------------//

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1)-1,ahi(1)

               ! Calculate convection-diffusion fluxes through each of the faces
               lflux = HALF * (fluxX(i+1,j,k-1) + fluxX(i+1,j,k  ))

               d_f = avg_h(avg_h(mu_g(i,j,k-1),mu_g(i+1,j,k-1)),&
                           avg_h(mu_g(i,j,k  ),mu_g(i+1,j,k  ))) * ayz_x

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

               lflux = HALF * (fluxY(i,j+1,k-1) + fluxY(i,j+1,k  ))

               d_f = avg_h(avg_h(mu_g(i,j-1,k),mu_g(i,j-1,k+1)),&
                           avg_h(mu_g(i,j  ,k),mu_g(i,j  ,k+1))) * axz_y

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

               lflux = HALF * (fluxZ(i,j,k) + fluxZ(i,j,k+1))

               d_f = mu_g(i,j,k) * axy_z

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

   end subroutine store_a_w_g0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: store_a_w_g1                                            !
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  !
!  The off-diagonal coefficients calculated here must be positive.     !
!  The center coefficient and the source vector are negative.          !
!  Implements higher order discretization.                             !
!  See source_w_g                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine store_a_w_g1(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
                              A_m, mu_g, u_g, v_g, w_g, &
                              fluxX, fluxY, fluxZ, &
                              dx, dy, dz)

      use functions, only: avg, avg_h
      use matrix, only: e, w, n, s, t, b

      use xsi, only: calc_xsi_x, calc_xsi_y, calc_xsi_z

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx, dy, dz

      ! Septadiagonal matrix A_W_g
      real(c_real), intent(INOUT) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(IN   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(IN   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(IN   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(IN   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: fluxZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Local variables
!---------------------------------------------------------------------//
      integer :: xlo(3),xhi(3)
      integer :: vello(3),velhi(3)
      integer :: i,j,k

      ! Diffusion parameter
      real(c_real) :: d_f, d_f_e, d_f_w, d_f_n, d_f_s

      ! Face mass flux
      real(c_real) :: lflux, lflux_e, lflux_w, lflux_n, lflux_s

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

      xlo(1) = alo(1)
      xlo(2) = alo(2)
      xlo(3) = alo(3)

      xhi(1) = ahi(1)+1
      xhi(2) = ahi(2)
      xhi(3) = ahi(3)

      vel(:,:,:) = 0.d0
      do k = wlo(3)+1,whi(3)-1
        do j = wlo(2),whi(2)
          do i = wlo(1),whi(1)
             vel(i,j,k) = avg(u_g(i+1,j,k-1), u_g(i+1,j,k  ))
          end do
        end do
      end do

      ! NOTES:   w_g  lives on z-faces    :   (lo(1): hi(1)  , lo(2):hi(2), lo(3):hi(3)+1)
      !          A_m  lives on z-faces    :   (lo(1): hi(1)  , lo(2):hi(2), lo(3):hi(3)+1)
      !          vel  lives on x-z edges  :   (lo(1): hi(1)+1, lo(2):hi(2), lo(3):hi(3)+1)
      !          xsi  lives on x-z edges  :   (lo(1): hi(1)+1, lo(2):hi(2), lo(3):hi(3)+1)

      allocate(xsi_(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_x (w_g, wlo, whi, vel, vello, velhi, xsi_, xlo, xhi, .true.)

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               lflux_e = half * (fluxX(i+1,j,k-1) + fluxX(i+1,j,k))
               lflux_w = half * (fluxX(i  ,j,k-1) + fluxX(i  ,j,k))

               d_f_e = avg_h(avg_h(mu_g(i,j,k-1),mu_g(i+1,j,k-1)),&
                             avg_h(mu_g(i,j,k  ),mu_g(i+1,j,k  ))) * ayz_x

               d_f_w = avg_h(avg_h(mu_g(i-1,j,k-1),mu_g(i,j,k-1)),&
                             avg_h(mu_g(i-1,j,k  ),mu_g(i,j,k  ))) * ayz_x

               A_m(i,j,k,e) = d_f_e - lflux_e*(      xsi_(i+1,j,k))
               A_m(i,j,k,w) = d_f_w + lflux_w*(one - xsi_(i  ,j,k))

            enddo
         enddo
      enddo
      deallocate(xsi_)

!---------------------------------------------------------------------//

      xlo(1) = alo(1)
      xlo(2) = alo(2)
      xlo(3) = alo(3)

      xhi(1) = ahi(1)
      xhi(2) = ahi(2)+1
      xhi(3) = ahi(3)

      vel(:,:,:) = 0.d0
      do k = wlo(3)+1,whi(3)-1
        do j = wlo(2),whi(2)
          do i = wlo(1),whi(1)
             vel(i,j,k) = avg(v_g(i,j,k-1), v_g(i,j,k  ))
          end do
        end do
      end do

      ! NOTES:   w_g  lives on z-faces   :   (lo(1): hi(1)  , lo(2):hi(2), lo(3):hi(3)+1)
      !          A_m  lives on z-faces   :   (lo(1): hi(1)  , lo(2):hi(2), lo(3):hi(3)+1)
      !          vel  lives on y-z edges :   (lo(1): hi(1), lo(2):hi(2)+1, lo(3):hi(3)+1)
      !          xsi  lives on y-z edges :   (lo(1): hi(1), lo(2):hi(2)+1, lo(3):hi(3)+1)

      allocate(xsi_(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_y (w_g, wlo, whi, vel, vello, velhi, xsi_, xlo, xhi, .true.)

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               lflux_n = half * (fluxY(i,j+1,k-1) + fluxY(i,j+1,k  ))
               lflux_s = half * (fluxY(i,j  ,k-1) + fluxY(i,j  ,k  ))

               d_f_n = avg_h(avg_h(mu_g(i,j,k-1),mu_g(i,j+1,k-1)),&
                             avg_h(mu_g(i,j,k  ),mu_g(i,j+1,k  ))) * axz_y

               d_f_s = avg_h(avg_h(mu_g(i,j-1,k-1),mu_g(i,j,k-1)),&
                             avg_h(mu_g(i,j-1,k  ),mu_g(i,j,k  ))) * axz_y

               A_m(i,j,k,n) = d_f_n - lflux_n*(      xsi_(i,j+1,k))
               A_m(i,j,k,s) = d_f_s + lflux_s*(one - xsi_(i,j  ,k))

            end do
         end do
      end do
      deallocate(xsi_)

!---------------------------------------------------------------------//

      xlo(1) = alo(1)
      xlo(2) = alo(2)
      xlo(3) = alo(3)-1

      xhi(1) = ahi(1)
      xhi(2) = ahi(2)
      xhi(3) = ahi(3)

      vel(:,:,:) = 0.d0
      do k = wlo(3),whi(3)-1
        do j = wlo(2),whi(2)
          do i = wlo(1),whi(1)
             vel(i,j,k) = avg(w_g(i,j,k), w_g(i,j,k+1))
          end do
        end do
      end do

      ! NOTES:   w_g  lives on z-faces   :   (lo(1): hi(1)  , lo(2):hi(2), lo(3):hi(3)+1)
      !          A_m  lives on z-faces   :   (lo(1): hi(1)  , lo(2):hi(2), lo(3):hi(3)+1)
      !          vel  lives on cell ctrs :   (lo(1): hi(1)  , lo(2):hi(2)  , lo(3):hi(3))
      !          xsi  lives on cell ctrs :   (lo(1): hi(1)  , lo(2):hi(2)  , lo(3):hi(3))
      !               (so xsi(i,j,k-1) and xsi(i,j,k) contribute to A_m(i,j,k))

      allocate( xsi_(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) )
      call calc_xsi_z (w_g, wlo, whi, vel, vello, velhi, xsi_, xlo, xhi, .false.)

      do k = alo(3)-1,ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               lflux = half * (fluxZ(i,j,k) + fluxZ(i,j,k+1))
               d_f = mu_g(i,j,k) * axy_z

               if (k.ge.alo(3)) A_m(i,j,k,  t) = d_f - lflux*(      xsi_(i,j,k))
               if (k.lt.ahi(3)) A_m(i,j,k+1,b) = d_f + lflux*(one - xsi_(i,j,k))

            enddo
         enddo
       enddo
      deallocate(xsi_)

!---------------------------------------------------------------------//

      deallocate(vel)

   end subroutine store_a_w_g1

end module w_g_conv_dif
