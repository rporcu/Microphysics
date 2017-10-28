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
   ! Compute new dt
   !
   subroutine compute_new_dt ( umax, vmax, wmax, romin, mumax, dx, dt ) &

   ! subroutine compute_new_dt ( umax, vmax, wmax, fgdsumax, fgdsvmax fgdswmax, &
   !      dragumax, dragvmax, dragwmax, mumax, romin, dx, dt ) &
        bind(C, name = "compute_new_dt")

      use constant, only: gravity 
      
      real(ar),       intent(in   ) :: umax, vmax, wmax
      ! real(ar),       intent(in   ) :: fgdsumax, fgdsvmax fgdswmax
      ! real(ar),       intent(in   ) :: dragumax, dragvmax, dragwmax
      real(ar),       intent(in   ) :: mumax, romin
      real(ar),       intent(in   ) :: dx(3)
      real(ar),       intent(inout) :: dt

      real(ar)                      :: uodx, vody, wodz
      real(ar)                      :: odx, ody, odz
      real(ar)                      :: fp_x, fp_y, fp_z
      real(ar)                      :: dt_c, dt_v
      real(ar)                      :: dt_gx, dt_gy, dt_gz, dt_g
      real(ar)                      :: dt_px, dt_py, dt_pz, dt_p
      real(ar),       parameter     :: cfl = 0.9_ar
      real(ar),       parameter     :: two = 2.0_ar, four = two*two
      real(ar),       parameter     :: eps = epsilon (zero)
      real(ar),       parameter     :: small = 1.0D-8
      
      odx  = one / dx(1)
      ody  = one / dx(2)
      odz  = one / dx(3)
      uodx = umax * odx
      vody = vmax * ody
      wodz = wmax * odz

      ! Convection
      ! Smaller components of velocity are not accounted for in the
      ! computation of the new time step (see IAMR )
      dt_c = dt
      if ( umax > small ) dt_c = min ( dt_c, cfl / uodx )
      if ( vmax > small ) dt_c = min ( dt_c, cfl / vody )
      if ( wmax > small ) dt_c = min ( dt_c, cfl / wodz )

      ! Viscous
      dt_v  = cfl * romin / ( two * ( mumax + eps ) * &
             & max ( odx*odx, ody*ody, odz*odz ) )

      !  Gravity
      dt_gx = cfl * two / helper ( uodx, gravity(1), odx ) 

      dt_gy = cfl * two / helper ( vody, gravity(2), ody ) 

      dt_gz = cfl * two / helper ( wodz, gravity(3), odz )

      dt_g  = min ( dt_gx, dt_gy, dt_gz )

      ! ! Particle-fluid momentum exchange
      ! fp_x  = fgdsumax * umax - dragumax


      ! dt_px = cfl * two / helper ( uodx, fp_x, odx ) 

      ! dt_py = cfl * two / helper ( vody, fp_y, ody ) 

      ! dt_pz = cfl * two / helper ( wodz, fp_z, odz )

      ! dt_p  = min ( dt_px, dt_py, dt_pz )

      ! print*, "umax, vmax, wmax = ", umax, vmax, wmax
      ! print*, "uodx, uody, uodz = ", uodx, vody, wodz
      ! print*, "dx               = ", dx
      ! print*, "dt, dt_c         = ", dt, dt_c

      dt = min ( dt, dt_c ) !, dt_v, dt_g )

   contains


      ! Compute root = b + sqrt ( b^2 + 4*f*odx )
      
      function helper ( b, f, odx )  result (res)

         real(ar), intent(in   ) :: b, f, odx
         real(ar)                :: res

         res = b + sqrt ( b*b + four * f * odx ) 

      end function helper

   end subroutine compute_new_dt



   ! 
   ! Given a pressure increment phi = scale * ( p^(n+1) - p^(n) )
   ! and the pressure at the previous time step, p^(n), this routine
   ! performs the update
   !
   !                   p^(n+1) = phi/scale +  p^(n)
   !
   !  where scale is a chosen scaling factor
   !
   subroutine update_pressure ( lo, hi, phi, slo, shi, pg, scale ) bind(C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: slo(3), shi(3)

      ! Scaling factor
      real(ar),        intent(in   ) :: scale
      
      ! Arrays
      real(ar),        intent(inout) ::                           &
           & pg(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
   
      real(ar),        intent(in   ) ::                           &
           & phi(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: oscale

      oscale = one / scale
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               pg(i,j,k) = oscale * phi(i,j,k) + pg(i,j,k)              
            end do
         end do
      end do
 
      
   end subroutine update_pressure



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

               u_i(i,j,k) = u_i(i,j,k) - codx * beta(i,j,k) * &
                    &      ( phi(i,j,k) - phi(i-i0,j-j0,k-k0) )
         
              
            end do
         end do
      end do

   end subroutine add_gradient

   

   
   !
   ! Computes  u_i* = u_i + dt * RHS
   ! where RHS is the Right-Hand Side of the momentum equation
   !
   ! Here u_i is the i-th component of velocity.
   ! Previous pressure gradient is not included.
   !
   ! dir = 1,2,3 indicates x, y, z direction respectively.
   !  
   subroutine compute_intermediate_velocity ( lo, hi, u_i_star, ulo, uhi,     &
        & u_i, ugradu_i,  divtau_i, drag_i, f_gds_i, rop, slo, shi, dt, dir ) &
        & bind(C, name="compute_intermediate_velocity")

      use constant, only: gravity
      
      integer(c_int),  intent(in)    :: lo(3),  hi(3)   ! Tile indeces for mf associated to ug
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      integer(c_int),  intent(in   ) :: dir
      real(ar),        intent(in   ) :: dt
      
      
      real(ar),        intent(inout) ::                           &
           & u_i_star(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & u_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),      &
           & divtau_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & ugradu_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), & 
           & drag_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),   &
           & f_gds_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),  &
           & rop(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
           
      integer(c_int)                 :: i, j, k, i0, j0, k0
      real(ar)                       :: rhs, irho, f 

      
      i0 = e_i(dir,1)

      j0 = e_i(dir,2)

      k0 = e_i(dir,3)
      
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               irho = half / rop(i,j,k) + half / rop(i-i0,j-j0,k-k0) 

               
               ! Remeber the minus sign between the two parts
               f    = f_gds_i(i,j,k) * u_i(i,j,k) - drag_i(i,j,k)
               
               rhs  = - ugradu_i(i,j,k) ! + gravity(dir)  + &
                     ! & irho * ( divtau_i(i,j,k) ) !+ f )

               u_i_star(i,j,k) = u_i(i,j,k) + dt * rhs      
              
               
            end do
         end do
      end do

   end subroutine compute_intermediate_velocity



   !
   ! Compute the RHS of the pressure poisson equation: rhs = - div(u*)/dt
   ! 
   subroutine compute_ppe_rhs ( lo, hi, rhs, slo, shi, u_g, ulo, uhi, v_g, vlo, vhi, &
        & w_g, wlo, whi, dx, dt, sum_rhsdv )  bind(C, name="compute_ppe_rhs")

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      real(ar),       intent(in   ) :: dx(3), dt

      real(ar),       intent(  out) :: sum_rhsdv
      
      real(ar),       intent(  out) :: &
           rhs(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      real(ar),       intent(in   ) :: &
           u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer             :: i, j, k
      real(ar)            :: odtdx, odtdy, odtdz, dv
      real(ar), parameter :: two = one / half 
      
      odtdx = one / ( dt * dx(1) )
      odtdy = one / ( dt * dx(2) )
      odtdz = one / ( dt * dx(3) )

      dv    = dx(1) * dx(2) * dx(3)

      sum_rhsdv = zero
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhs(i,j,k) =  -  ( ( u_g(i+1,j,k) - u_g(i,j,k) ) * odtdx + &
                                  ( v_g(i,j+1,k) - v_g(i,j,k) ) * odtdy + &
                                  ( w_g(i,j,k+1) - w_g(i,j,k) ) * odtdz )
               
               sum_rhsdv = sum_rhsdv + rhs(i,j,k) * dv
            end do
         end do
      end do

   end subroutine compute_ppe_rhs




   
   !
   ! Apply the pressure correction u = u^* - dt (ep_g/rop_g) grad(p)
   ! Note that the scalar ep_g/rop_g is 1/ro_g, hence the name oro_g 
   ! 
   subroutine apply_pressure_correction ( utlo, uthi, ug, ulo, uhi, oro_x,  &
        & vtlo, vthi, vg, vlo, vhi, oro_y, wtlo, wthi, wg, wlo, whi, oro_z, &
        & pg, slo, shi, dx, dt )  bind(C, name="apply_pressure_correction")

      ! Array bounds
      integer(c_int),   intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),   intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),   intent(in   ) :: wlo(3), whi(3)
      integer(c_int),   intent(in   ) :: slo(3), shi(3)
      
      ! Tile bounds
      integer(c_int),   intent(in   ) :: utlo(3), uthi(3)
      integer(c_int),   intent(in   ) :: vtlo(3), vthi(3)
      integer(c_int),   intent(in   ) :: wtlo(3), wthi(3)
      
      ! Grid and time step
      real(ar),         intent(in   ) :: dt
      real(ar),         intent(in   ) :: dx(3)
           
      ! Arrays
      real(ar),         intent(inout) ::                    &
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), & 
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(ar),         intent(in   ) ::                    &
           & oro_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), & 
           & oro_y(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           & oro_z(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      
      real(ar),         intent(in   ) ::                    &
           & pg(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) 

      ! Local variables
      integer      :: i, j, k
      real(ar)     :: dtodx, dtody, dtodz

      dtodx = dt / dx(1) 
      dtody = dt / dx(2) 
      dtodz = dt / dx(3) 
      
      ! x component
      do k = utlo(3), uthi(3)
         do j = utlo(2), uthi(2)
            do i = utlo(1), uthi(1)
               ug(i,j,k) = ug(i,j,k) - dtodx * oro_x(i,j,k) * &
                    ( pg(i,j,k) - pg(i-1,j,k) )   
            end do
         end do
      end do


      ! y component
      do k = vtlo(3), vthi(3)
         do j = vtlo(2), vthi(2)
            do i = vtlo(1), vthi(1)
               vg(i,j,k) = vg(i,j,k) - dtody * oro_y(i,j,k) * &
                    ( pg(i,j,k) - pg(i,j-1,k) )   
            end do
         end do
      end do

      
      ! z component
      do k = wtlo(3), wthi(3)
         do j = wtlo(2), wthi(2)
            do i = wtlo(1), wthi(1)
               wg(i,j,k) = wg(i,j,k) - dtodz * oro_z(i,j,k) * &
                    ( pg(i,j,k) - pg(i,j,k-1) )   
            end do
         end do
      end do

      
   end subroutine apply_pressure_correction

   
   !
   ! Compute the coefficients of the PPE, i.e. 1 / ro_g = eps_g/rho_g,
   ! at the faces of the pressure cells along the x-axis
   ! (u-velocity locations).
   ! 
   subroutine compute_oro_g_x ( lo, hi, oro_g_x, ulo, uhi, &
        rop_g, slo, shi, ep_g)  bind(C, name="compute_oro_g_x")

      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      
      real(ar),       intent(in   ) :: &
           rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
            ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
           oro_g_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      
      integer      :: i, j, k


      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               oro_g_x(i,j,k) = half * (              &
                    ep_g(i,j,k) / rop_g(i,j,k)     + &
                    ep_g(i-1,j,k) / rop_g(i-1,j,k) )
            end do
         end do
      end do

   end subroutine compute_oro_g_x

   
   !
   ! Compute the coefficients of the PPE, i.e. 1 / ro_g = eps_g/rho_g,
   ! at the faces of the pressure cells along the y-axis
   ! (v-velocity locations).
   ! 
   subroutine compute_oro_g_y ( lo, hi, oro_g_y, vlo, vhi, &
        rop_g, slo, shi, ep_g)  bind(C, name="compute_oro_g_y")
      
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      
      real(ar),       intent(in   ) :: &
           rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      real(ar),       intent(  out) :: &
           oro_g_y(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      
      integer      :: i, j, k
      
      
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               oro_g_y(i,j,k) = half * (              &
                    ep_g(i,j,k) / rop_g(i,j,k)     + &
                    ep_g(i,j-1,k) / rop_g(i,j-1,k) )
            end do
         end do
      end do
      
   end subroutine compute_oro_g_y


   
   !
   ! Compute the coefficients of the PPE, i.e. 1 / ro_g = eps_g/rho_g,
   ! at the faces of the pressure cells along the z-axis
   ! (w-velocity locations).
   ! 
   subroutine compute_oro_g_z ( lo, hi, oro_g_z, wlo, whi, &
        rop_g, slo, shi, ep_g)  bind(C, name="compute_oro_g_z")
      
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: wlo(3),whi(3)
      
      real(ar),       intent(in   ) :: &
           rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      real(ar),       intent(  out) :: &
           oro_g_z(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      
      integer      :: i, j, k
      
      
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               oro_g_z(i,j,k) = half * (              &
                    ep_g(i,j,k) / rop_g(i,j,k)     + &
                    ep_g(i,j,k-1) / rop_g(i,j,k-1) )
            end do
         end do
      end do
      
   end subroutine compute_oro_g_z

   
end module projection_mod
