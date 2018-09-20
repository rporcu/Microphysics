!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_particle_beta                                         !
!                                                                      !
!  Purpose: This routine is called before the FLUID solve.             !
!  It calculates the source terms for the center coefficients and RHS  !
!  for the momentum equations. It also saves the drag coefficient for  !
!  each particle.                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine calc_particle_beta ( slo, shi, ep_g, ro_g, vel_g, mu_g, &
                                np, particles, dx, dy, dz) &
   bind(C, name="calc_particle_beta")

   use amrex_fort_module,  only : rt => amrex_real
   use iso_c_binding ,     only: c_int
   use des_drag_gp_module, only: des_drag_gp
   use particle_mod,      only: particle_t

   implicit none

   integer(c_int), intent(in   ) :: slo(3),shi(3)
   integer(c_int), intent(in   ) :: np

   real(rt), intent(in   ) :: &
         ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        vel_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

   real(rt), intent(in   ) :: dx, dy, dz

   type(particle_t), intent(inout) :: particles(np)

   ! Local variables
   !---------------------------------------------------------------------//
   ! Loop counters: Particle, fluid cell, neighbor cells
   integer      :: p, i, j, k
   real(rt) :: velfp(3), velp(3), beta
   real(rt) :: odx, ody, odz, ovol
   real(rt) :: lx, ly, lz
   real(rt) :: sx_lo, sy_lo, sz_lo
   real(rt) :: sx_hi, sy_hi, sz_hi
   real(rt) :: plo(3)
   !......................................................................!

   odx = 1.0d0/dx
   ody = 1.0d0/dy
   odz = 1.0d0/dz

   ovol = 1.0d0/(dx*dy*dz)

   plo = 0.0d0 ! HACK -- This should get passed into routine

   ! Calculate the gas phase forces acting on each particle.

   do p = 1, np

      lx = (particles(p) % pos(1) - plo(1))*odx + 0.5d0
      ly = (particles(p) % pos(2) - plo(2))*ody + 0.5d0
      lz = (particles(p) % pos(3) - plo(3))*odz + 0.5d0

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      sx_hi = lx - i;  sx_lo = 1.0d0 - sx_hi
      sy_hi = ly - j;  sy_lo = 1.0d0 - sy_hi
      sz_hi = lz - k;  sz_lo = 1.0d0 - sz_hi

      ! Fluid velocity at the particle's position.
      velfp(1:3) = sx_lo*sy_lo*sz_lo*vel_g(i-1, j-1, k-1,1:3) + &
                   sx_lo*sy_lo*sz_hi*vel_g(i-1, j-1, k  ,1:3) + &
                   sx_lo*sy_hi*sz_lo*vel_g(i-1, j  , k-1,1:3) + &
                   sx_lo*sy_hi*sz_hi*vel_g(i-1, j  , k  ,1:3) + &
                   sx_hi*sy_lo*sz_lo*vel_g(i  , j-1, k-1,1:3) + &
                   sx_hi*sy_lo*sz_hi*vel_g(i  , j-1, k  ,1:3) + &
                   sx_hi*sy_hi*sz_lo*vel_g(i  , j  , k-1,1:3) + &
                   sx_hi*sy_hi*sz_hi*vel_g(i  , j  , k  ,1:3)

      velp(:) = particles(p) % vel

      ! Calculate drag coefficient, beta
      call des_drag_gp(slo, shi, p, velp, velfp, ep_g(i,j,k),  &
           ro_g, mu_g, beta, i, j, k, particles(p) % radius,   &
           particles(p) % volume, particles(p) % density,      &
           particles(p) % phase )

      particles(p) % drag(1) = beta

   end do

end subroutine calc_particle_beta

