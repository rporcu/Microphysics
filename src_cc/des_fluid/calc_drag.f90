!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_drag_particle                                      !
!                                                                      !
!  Purpose: This routine must be called after calc_particle_beta       !
!  because it assumes you have already computed beta and stored it in  !
!  particles(p)%drag(1).   Here we compute the actual drag on the      !
!  particle.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine calc_drag_particle( gp   , gplo, gphi, &
                               gp0  , gp0lo, gp0hi, &
                               vel_g, ulo, uhi, &
                               np, particles, dx, dy, dz) &
     bind(C, name="calc_drag_particle")

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int
   use des_drag_gp_module, only: des_drag_gp
   use particle_mod, only: particle_t

   implicit none

   integer(c_int), intent(in   ) ::  ulo(3), uhi(3)
   integer(c_int), intent(in   ) ::  gplo(3), gphi(3)
   integer(c_int), intent(in   ) ::  gp0lo(3), gp0hi(3)
   integer(c_int), intent(in   ) ::  np

   real(rt)  , intent(in   ) :: &
        vel_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3), &
           gp(gplo(1):gphi(1),gplo(2):gphi(2),gplo(3):gphi(3),3),  &
          gp0(gp0lo(1):gp0hi(1),gp0lo(2):gp0hi(2),gp0lo(3):gp0hi(3),3)

   real(rt),     intent(in   ) :: dx, dy, dz
   type(particle_t), intent(inout) :: particles(np)


   ! Local variables
   !---------------------------------------------------------------------//
   ! Loop counters: Particle, fluid cell, neighbor cells
   integer      :: p, i, j, k
   real(rt) :: velfp(3), gradpg(3)
   real(rt) :: beta(np)
   real(rt) :: odx, ody, odz
   real(rt) :: lx, ly, lz
   real(rt) :: sx_lo, sy_lo, sz_lo
   real(rt) :: sx_hi, sy_hi, sz_hi
   real(rt) :: ux_lo, vy_lo, wz_lo
   real(rt) :: ux_hi, vy_hi, wz_hi
   real(rt) :: plo(3)

   plo = 0.0d0 ! HACK -- This should get passed into routine

   odx = 1.0d0/dx
   ody = 1.0d0/dy
   odz = 1.0d0/dz

   do p = 1, np
      beta(p) = particles(p) % drag(1)
   end do

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

      gradpg(1:3) = sx_lo*sy_lo*sz_lo*(gp(i-1, j-1, k-1,1:3) + gp0(i-1, j-1, k-1,1:3)) + &
                    sx_lo*sy_lo*sz_hi*(gp(i-1, j-1, k  ,1:3) + gp0(i-1, j-1, k  ,1:3)) + &
                    sx_lo*sy_hi*sz_lo*(gp(i-1, j  , k-1,1:3) + gp0(i-1, j  , k-1,1:3)) + &
                    sx_lo*sy_hi*sz_hi*(gp(i-1, j  , k  ,1:3) + gp0(i-1, j  , k  ,1:3)) + &
                    sx_hi*sy_lo*sz_lo*(gp(i  , j-1, k-1,1:3) + gp0(i  , j-1, k-1,1:3)) + &
                    sx_hi*sy_lo*sz_hi*(gp(i  , j-1, k  ,1:3) + gp0(i  , j-1, k  ,1:3)) + &
                    sx_hi*sy_hi*sz_lo*(gp(i  , j  , k-1,1:3) + gp0(i  , j  , k-1,1:3)) + &
                    sx_hi*sy_hi*sz_hi*(gp(i  , j  , k  ,1:3) + gp0(i  , j  , k  ,1:3))

      particles(p) % drag = beta(p)*(velfp - particles(p) % vel) - &
                            gradpg(:) * particles(p) % volume

   enddo

end subroutine calc_drag_particle
