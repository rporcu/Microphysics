! 
!              
!  Velocity field initialization routines for testing of 
!  projection method.
!  
!  Author: Michele Rosso
! 
!  Date: October 19, 2017
!
! 
module projection_tests_mod
   
   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   
   implicit none
   private

contains


   !
   ! plane = 1  --> init xy
   ! plane = 2  --> init yz
   ! plane = 3  --> init xz
   ! 

   subroutine init_periodic_vorteces ( utlo, uthi, u, ulo, uhi,  &
        & vtlo, vthi, v, vlo, vhi, wtlo, wthi, w, wlo, whi, dx,  &
        & domlo,  plane )  bind(C)

      ! Array bounds
      integer(c_int),   intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),   intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),   intent(in   ) :: wlo(3), whi(3)

      ! Tile bounds
      integer(c_int),   intent(in   ) :: utlo(3), uthi(3)
      integer(c_int),   intent(in   ) :: vtlo(3), vthi(3)
      integer(c_int),   intent(in   ) :: wtlo(3), wthi(3)
      
      ! Grid and domain lower bound
      integer(c_int),   intent(in   ) :: domlo(3)
      real(ar),         intent(in   ) :: dx(3)

      ! Plane
      integer(c_int),   intent(in   ) :: plane
      
      ! Arrays
      real(ar),         intent(inout) ::                   &
           & u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), & 
           & v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           & w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      ! Local variables
      integer(c_int)                  :: i, j, k
      real(ar)                        :: x, y, z
      real(ar)                        :: twopi = 8.0_ar * atan(one) 


      select case ( plane )

      case (1)  ! x-y plane
         
         ! x-direction
         do k = utlo(3), uthi(3)
            do j = utlo(2), uthi(2)
               y =  ( real(j,ar) + half ) * dx(2)
               do i = utlo(1), uthi(1)
                  u(i,j,k) = tanh ( 30.0_ar * (0.25_ar - abs ( y - 0.5_ar ) ) )
               end do
            end do
         end do
         
         ! y-direction
         do k = vtlo(3), vthi(3) 
            do j = vtlo(2), vthi(2)
               do i = vtlo(1), vthi(1)
                  x =  ( real(i,ar) + half ) * dx(1)
                  v(i,j,k) = 0.05_ar * sin ( twopi * x )
               end do
            end do
         end do
         
         ! z-direction
         do k = wtlo(3), wthi(3) 
            do j = wtlo(2), wthi(2)
               do i = wtlo(1), wthi(1)
                  w(i,j,k) = zero
               end do
            end do
         end do
         
      case (2)  ! x-z plane
         
         ! x-direction
         do k = utlo(3), uthi(3)
            z =  ( real(k,ar) + half ) * dx(3)
            do j = utlo(2), uthi(2)
               do i = utlo(1), uthi(1)
                  u(i,j,k) = tanh ( 30.0_ar * (0.25_ar - abs ( z - 0.5_ar ) ) )
               end do
            end do
         end do
         
         ! y-direction
         do k = vtlo(3), vthi(3) 
            do j = vtlo(2), vthi(2)
               do i = vtlo(1), vthi(1)
                  v(i,j,k) = zero
               end do
            end do
         end do
         
         ! z-direction
         do k = wtlo(3), wthi(3) 
            do j = wtlo(2), wthi(2)
               do i = wtlo(1), wthi(1)
                  x =  ( real(i,ar) + half ) * dx(1)
                  w(i,j,k) = 0.05_ar * sin ( twopi * x )
               end do
            end do
         end do

         
      case (3)  ! y-z plane
         
         ! x-direction
         do k = utlo(3), uthi(3)
            do j = utlo(2), uthi(2)
               do i = utlo(1), uthi(1)
                  u(i,j,k) = zero 
               end do
            end do
         end do
         
         ! y-direction
         do k = vtlo(3), vthi(3) 
            z =  ( real(k,ar) + half ) * dx(3)
            do j = vtlo(2), vthi(2)
               do i = vtlo(1), vthi(1)
                  v(i,j,k) = tanh ( 30.0_ar * (0.25_ar - abs ( z - 0.5_ar ) ) )
               end do
            end do
         end do
         
         ! z-direction
         do k = wtlo(3), wthi(3) 
            do j = wtlo(2), wthi(2)
               y =  ( real(j,ar) + half ) * dx(2)
               do i = wtlo(1), wthi(1)
                  w(i,j,k) = 0.05_ar * sin ( twopi * y )
               end do
            end do
         end do

      end select 
         
   end subroutine init_periodic_vorteces

end module projection_tests_mod
