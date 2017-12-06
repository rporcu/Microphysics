! 
!              
!  This module contains the subroutines to perform some of the steps of the
!  projection method.
!  
!  Author: Michele Rosso
! 
!  Date: October 13, 2017
!
! 
module projection_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one

   implicit none
   private

   ! Define here the unit vectors
   ! This is used to shift index  based on how the variable is staggered
   ! Check e_x, e_y and e_z in mfix_level.H
   integer(c_int), parameter :: e_i(3,3) = reshape ( [1,0,0,0,1,0,0,0,1], [3,3] )  

contains

   !
   ! Compute new dt by using the formula derived in 
   ! "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
   ! by Kang et al. (JCP).
   !
   !  dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )
   !
   ! where
   ! 
   ! C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection
   ! 
   ! V = 2 * max(mu/ro) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion
   !
   ! Fx, Fy, Fz = net acceleration due to external forces
   ! 
   ! WARNING: We use a slightly modified version of C in the implementation below
   ! 
   subroutine compute_new_dt ( umax, vmax, wmax, romin, mumax, dx, cfl, dt ) &
        & bind(C)
      
      ! subroutine compute_new_dt ( umax, vmax, wmax, fgdsumax, fgdsvmax fgdswmax, &
      !      dragumax, dragvmax, dragwmax, mumax, romin, dx, dt ) &
      

      use constant, only: gravity 

      real(ar),       intent(in   ) :: umax, vmax, wmax
      ! real(ar),       intent(in   ) :: fgdsumax, fgdsvmax fgdswmax
      ! real(ar),       intent(in   ) :: dragumax, dragvmax, dragwmax
      real(ar),       intent(in   ) :: mumax, romin
      real(ar),       intent(in   ) :: dx(3), cfl
      real(ar),       intent(inout) :: dt
      real(ar)                      :: old_dt
      real(ar)                      :: c_cfl, v_cfl, f_cfl
      real(ar)                      :: odx, ody, odz
      real(ar)                      :: fp_x, fp_y, fp_z, tmp
      real(ar),       parameter     :: two = 2.0_ar, four = two*two
      real(ar),       parameter     :: eps = epsilon (zero)
      
      odx    = one / dx(1)
      ody    = one / dx(2)
      odz    = one / dx(3)
      c_cfl  = zero
      v_cfl  = zero
      f_cfl  = zero
      old_dt = dt
      
      ! Convection
      ! c_cfl = umax*odx + vmax*ody + wmax*odz  <-- Too restricting
      c_cfl = max ( umax*odx, vmax*ody, wmax*odz)
      
      ! Viscous 
      v_cfl = two * ( mumax / romin ) * ( odx**2 + ody**2 + odz**2 )
              
      ! Gravity
      f_cfl = f_cfl + gravity(1) * odx + gravity(2) * ody &
           &        + gravity(3) * odz
      
      ! Put all together
      tmp = (c_cfl + v_cfl)  + sqrt ( (c_cfl + v_cfl)**2 + four * f_cfl )
      dt  = cfl * two / tmp

      ! Protect against tmp very small
      ! This may happen, for example, when the initial velocity field
      ! is zero for an inviscid flow with no external forcing  
      if ( tmp <= eps ) then
         dt = old_dt
      else 
         dt = min ( dt, old_dt )
      end if


      
   end subroutine compute_new_dt




   !
   ! Computes the term RHS = - div (uu) + div (tau)/rop
   ! along direction "dir"
   !
   ! dir = 1, 2, 3 ( 1=x, 2=y, 3=z ) 
   !
   subroutine compute_fluid_acceleration ( lo, hi, rhs, rlo, rhi, sl, &
        & u, ulo, uhi, v, vlo, vhi, w, wlo, whi, mu, slo, shi, rop,   &
        & dx, dir, order ) bind(C)

      use convection_mod
      use diffusion_mod

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: rlo(3), rhi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: wlo(3), whi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! Scheme order
      integer(c_int), intent(in   ) :: order

      ! Grid 
      real(ar),       intent(in   ) :: dx(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) ::                       &
           &   u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           &   v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           &   w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           &  mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           & rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           &  sl(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

      real(ar),       intent(inout) ::                       &
           & rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

      ! Local working arrays
      real(ar)                      ::                        &
           & conv(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3)), &
           & diff(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

      ! Local variables
      integer                       :: i, j, k

      ! Compute convection term
      select case ( dir )
      case (1)
         call compute_ugradu_x ( lo, hi, u, ulo, uhi, v, vlo, vhi, &
              & w, wlo, whi, sl, conv, dx, order )  

         call compute_divtau_x ( lo, hi, u, ulo, uhi, v, vlo, vhi, &
              & w, wlo, whi, mu, slo, shi, diff, dx )
      case(2)
         call compute_ugradu_y ( lo, hi, u, ulo, uhi, v, vlo, vhi,  &
              & w, wlo, whi, sl, conv, dx, order )

         call compute_divtau_y ( lo, hi, u, ulo, uhi, v, vlo, vhi, &
              & w, wlo, whi, mu, slo, shi, diff, dx )
      case(3)         
         call compute_ugradu_z ( lo, hi, u, ulo, uhi, v, vlo, vhi,  &
              & w, wlo, whi, sl, conv, dx, order )

         call compute_divtau_z ( lo, hi, u, ulo, uhi, v, vlo, vhi, &
              & w, wlo, whi, mu, slo, shi, diff, dx )
      end select

      !
      !
      !  REMEMBER TO DIVIDE div(tau) by rho !!!!!
      !
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhs(i,j,k) = -conv(i,j,k) + diff(i,j,k)
            end do
         end do
      end do

   end subroutine compute_fluid_acceleration


   !
   ! Computes  u_i = u_i + C * (beta) * (dphi/dx_i)
   !
   ! u_i  = i-th component of a staggered vector field u.
   !
   ! beta = scalar field defined at u_i location.
   !
   ! phi  = scalar field defined at cell centers
   ! 
   ! C    = real constant
   !
   ! dir  = 1,2,3 indicates x, y, z direction respectively.
   !  
   subroutine add_gradient ( lo, hi, u_i, ulo, uhi, beta,    &
        & phi, slo, shi, dx, c, dir ) bind (C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)

      ! Grid and time spacing
      real(ar),        intent(in   ) :: c, dx(3)

      ! Direction
      integer(c_int),  intent(in   ) :: dir

      ! Arrays
      real(ar),        intent(in   ) ::                        &
           & beta(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),  &
           &  phi(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),        intent(inout) ::                           &
           & u_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k, i0, j0, k0
      real(ar)                       :: codx

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      codx = c / dx(dir) 

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               u_i(i,j,k) = u_i(i,j,k) + codx * beta(i,j,k) * &
                    &      ( phi(i,j,k) - phi(i-i0,j-j0,k-k0) )


            end do
         end do
      end do

   end subroutine add_gradient



   !
   ! Compute the coefficients of the PPE, i.e. 1 / ro_g = eps_g/rho_g,
   ! at the faces of the pressure cells along the "dir"-axis.
   ! 
   subroutine compute_oro_g ( lo, hi, oro_g, alo, ahi, &
        rop_g, slo, shi, ep_g, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) :: &
           rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
           oro_g(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      integer      :: i, j, k, i0, j0, k0

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)


      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               oro_g(i,j,k) = half * ( ep_g(i,j,k) / rop_g(i,j,k) + &
                    & ep_g(i-i0,j-j0,k-k0) / rop_g(i-i0,j-j0,k-k0) )
            end do
         end do
      end do

   end subroutine compute_oro_g





end module projection_mod
