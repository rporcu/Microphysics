
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_particle_beta                                      !
!                                                                      !
!  Purpose: This routine is called before the FLUID solve.             !
!  It calculates the source terms for the center coefficients and RHS  !
!  for the momentum equations. It also saves the drag coefficient for  !
!  each particle.                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine calc_particle_beta ( ep_g, slo, shi, &
 &                              ro_g,  mu_g,    &
 &                              vel_g, ulo, uhi,& 
 &                              np, particles, x0, dx ) bind(C)

   use amrex_fort_module,               only: rt => amrex_real
   use iso_c_binding,                   only: c_int
   use des_drag_gp_module,              only: des_drag_gp
   use particle_mod,                    only: particle_t
   use param,                           only: zero, half, one
   use interpolation_m,                 only: trilinear_interp
   
   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)
   integer(c_int), intent(in   ) :: ulo(3), uhi(3)
   
   ! Array
   real(rt), intent(in   ) :: &
    &  ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    &  ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    &  mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    & vel_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

   ! Particles
   integer(c_int),   intent(in   ) :: np
   type(particle_t), intent(inout) :: particles(np)

   ! Coordinates of lower corner of domain
   real(rt),         intent(in   ) :: x0(3)  

   ! Grid
   real(rt),         intent(in   ) :: dx(3)


   ! Local variables
   !---------------------------------------------------------------------//
   ! Loop counters: Particle, fluid cell, neighbor cells
   integer  :: p, i, j, k
   real(rt) :: velfp(3), velp(3), beta
   real(rt) :: odx, ody, odz, ovol

   odx  = one/dx(1)
   ody  = one/dx(2)
   odz  = one/dx(3)
   ovol = odx*ody*odz

   ! Calculate the gas phase forces acting on each particle.
   do p = 1, np

      associate ( ppos => particles(p) % pos )
         
         velfp(:)  = trilinear_interp(vel_g, ulo, uhi, 3, ppos, x0, dx) 

         ! Indeces of cell where particle is located
         i = floor((ppos(1) - x0(1))*odx)
         j = floor((ppos(2) - x0(2))*ody)
         k = floor((ppos(3) - x0(3))*odz)

         ! Calculate drag coefficient, beta
         call des_drag_gp(slo, shi, p, particles(p) % vel, velfp, &
          ep_g(i,j,k), ro_g, mu_g, beta, i, j, k,             &
          particles(p) % radius,  particles(p) % volume,      &
          particles(p) % density, particles(p) % phase )

         particles(p) % drag(1) = beta
         
      end associate
      
   end do

end subroutine calc_particle_beta

