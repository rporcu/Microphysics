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

contains
   
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
           
      integer(c_int)                 :: i, j, k
      real(ar)                       :: rhs, irho, f 
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               irho = half / rop(i,j,k) + half / rop(i-1,j,k) 

               ! Remeber the minus sign between the two parts
               f    = f_gds_i(i,j,k) * u_i(i,j,k) - drag_i(i,j,k)
               
               rhs  = - ugradu_i(i,j,k) + gravity(dir) + &
                    & irho * ( divtau_i(i,j,k) + f )

               u_i_star(i,j,k) = u_i(i,j,k) + dt * rhs      
              
               
            end do
         end do
      end do

   end subroutine compute_intermediate_velocity



   !
   ! Compute the RHS of the pressure poisson equation: rhs = - div(u*)/dt
   ! 
   subroutine compute_ppe_rhs ( lo, hi, rhs, slo, shi, u_g, ulo, uhi, v_g, vlo, vhi, &
        & w_g, wlo, whi, dx, dt )  bind(C, name="compute_ppe_rhs")

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      real(ar),       intent(in   ) :: dx(3), dt
      
      real(ar),       intent(  out) :: &
           rhs(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      real(ar),       intent(in   ) :: &
           u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer      :: i, j, k
      real(ar)     :: odx, ody, odz, odt

      odx = 1.d0 / dx(1)
      ody = 1.d0 / dx(2)
      odz = 1.d0 / dx(3)
      odt = 1.d0 / dt

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rhs(i,j,k) = - odt * ( &
                    (u_g(i+1,j,k)-u_g(i,j,k))*odx + &
                    (v_g(i,j+1,k)-v_g(i,j,k))*ody + &
                    (w_g(i,j,k+1)-w_g(i,j,k))*odz )
            end do
         end do
      end do

   end subroutine compute_ppe_rhs


   !
   ! Apply the pressure correction u = u^* - dt (ep_g/rop_g) dp/dx
   ! Note that the scalar ep_g/rop_g is 1/ro_g, hence the name oro_g 
   ! 
   subroutine apply_pressure_correction_x ( lo, hi, p_g, slo, shi, &
     &  u_g, ulo, uhi, oro_g, dx, dt )  bind(C, name="apply_pressure_correction_x")

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)
      real(ar),       intent(in   ) :: dx(3), dt
      
      real(ar),       intent(in   ) :: &
           p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      
      real(ar),       intent(inout) :: &
           u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           oro_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
           
      integer      :: i, j, k
      real(ar)     :: dtodx

      dtodx = dt / dx(1)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               u_g(i,j,k) = u_g(i,j,k) -  dtodx *  oro_g(i,j,k) * &
                    ( p_g(i,j,k) - p_g(i-1,j,k) )   
            end do
         end do
      end do

   end subroutine apply_pressure_correction_x


   
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
