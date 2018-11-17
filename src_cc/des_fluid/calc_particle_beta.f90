
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
subroutine calc_particle_beta ( slo, shi, ep_g, ro_g, vel_g, mu_g,   &
 &                              np, particles, x0, dx ) bind(C)

   use amrex_fort_module,               only: rt => amrex_real
   use iso_c_binding,                   only: c_int
   use des_drag_gp_module,              only: des_drag_gp
   use particle_mod,                    only: particle_t
   use param,                           only: zero, half, one

   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)

   ! Array
   real(rt), intent(in   ) :: &
    &  ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    &  ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    &  mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    & vel_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

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
   real(rt) :: lx, ly, lz
   real(rt) :: sx_lo, sy_lo, sz_lo
   real(rt) :: sx_hi, sy_hi, sz_hi
   real(rt) :: wx, wy, wz

   odx  = one/dx(1)
   ody  = one/dx(2)
   odz  = one/dx(3)
   ovol = odx*ody*odz

   ! Calculate the gas phase forces acting on each particle.
   do p = 1, np


      ! Thir-linear interpolation
      ! By adding half we always pick the upper cell in the
      ! stencil
      lx = (particles(p) % pos(1) - x0(1))*odx + half
      ly = (particles(p) % pos(2) - x0(2))*ody + half
      lz = (particles(p) % pos(3) - x0(3))*odz + half

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      sx_hi = lx - i;  sx_lo = one - sx_hi
      sy_hi = ly - j;  sy_lo = one - sy_hi
      sz_hi = lz - k;  sz_lo = one - sz_hi

      velfp(1:3) = sx_lo*sy_lo*sz_lo*vel_g(i-1, j-1, k-1,1:3) + &
       &       sx_lo*sy_lo*sz_hi*vel_g(i-1, j-1, k  ,1:3) + &
       &       sx_lo*sy_hi*sz_lo*vel_g(i-1, j  , k-1,1:3) + &
       &       sx_lo*sy_hi*sz_hi*vel_g(i-1, j  , k  ,1:3) + &
       &       sx_hi*sy_lo*sz_lo*vel_g(i  , j-1, k-1,1:3) + &
       &       sx_hi*sy_lo*sz_hi*vel_g(i  , j-1, k  ,1:3) + &
       &       sx_hi*sy_hi*sz_lo*vel_g(i  , j  , k-1,1:3) + &
       &       sx_hi*sy_hi*sz_hi*vel_g(i  , j  , k  ,1:3)

      ! Calculate drag coefficient, beta
      call des_drag_gp(slo, shi, p, particles(p) % vel, velfp, &
       ep_g(i,j,k), ro_g, mu_g, beta, i, j, k,             &
       particles(p) % radius,  particles(p) % volume,      &
       particles(p) % density, particles(p) % phase )

      particles(p) % drag(1) = beta

   end do

end subroutine calc_particle_beta




subroutine calc_particle_beta_eb ( slo, shi, ep_g,  ro_g, vel_g, mu_g, &
 &                                 flags, flo, fhi, np, particles, x0, &
 &                                 dx ) bind(C)

   use amrex_ebcellflag_module,  only: is_covered_cell
   use amrex_fort_module,        only: rt => amrex_real
   use iso_c_binding,            only: c_int
   use des_drag_gp_module,       only: des_drag_gp
   use particle_mod,             only: particle_t
   use param,                    only: zero, half, one

   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)
   integer(c_int), intent(in   ) :: flo(3), fhi(3)

   ! Array
   real(rt), intent(in   ) :: &
    &  ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    &  ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    &  mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),  &
    & vel_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

   integer(c_int), intent(in   ) ::  &
    & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

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
   real(rt) :: lx, ly, lz
   real(rt) :: sx_lo, sy_lo, sz_lo
   real(rt) :: sx_hi, sy_hi, sz_hi
   real(rt) :: wx, wy, wz

   odx  = one/dx(1)
   ody  = one/dx(2)
   odz  = one/dx(3)
   ovol = odx*ody*odz

   ! Calculate the gas phase forces acting on each particle.
   do p = 1, np

      lx = (particles(p) % pos(1) - x0(1))*odx + half
      ly = (particles(p) % pos(2) - x0(2))*ody + half
      lz = (particles(p) % pos(3) - x0(3))*odz + half

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      ! Thir-linear interpolation if stencil is valid
      if ( ( .not. is_covered_cell(flags(i-1,j-1,k-1)) ) .and. &
       &   ( .not. is_covered_cell(flags(i-1,j-1,k  )) ) .and. &
       &   ( .not. is_covered_cell(flags(i-1,j  ,k-1)) ) .and. &
       &   ( .not. is_covered_cell(flags(i-1,j  ,k  )) ) .and. &
       &   ( .not. is_covered_cell(flags(i  ,j-1,k-1)) ) .and. &
       &   ( .not. is_covered_cell(flags(i  ,j-1,k  )) ) .and. &
       &   ( .not. is_covered_cell(flags(i  ,j  ,k-1)) ) .and. &
       &   ( .not. is_covered_cell(flags(i  ,j  ,k  )) ) ) then

         sx_hi = lx - i;  sx_lo = one - sx_hi
         sy_hi = ly - j;  sy_lo = one - sy_hi
         sz_hi = lz - k;  sz_lo = one - sz_hi

         velfp(1:3) = sx_lo*sy_lo*sz_lo*vel_g(i-1, j-1, k-1,1:3) + &
          &       sx_lo*sy_lo*sz_hi*vel_g(i-1, j-1, k  ,1:3) + &
          &       sx_lo*sy_hi*sz_lo*vel_g(i-1, j  , k-1,1:3) + &
          &       sx_lo*sy_hi*sz_hi*vel_g(i-1, j  , k  ,1:3) + &
          &       sx_hi*sy_lo*sz_lo*vel_g(i  , j-1, k-1,1:3) + &
          &       sx_hi*sy_lo*sz_hi*vel_g(i  , j-1, k  ,1:3) + &
          &       sx_hi*sy_hi*sz_lo*vel_g(i  , j  , k-1,1:3) + &
          &       sx_hi*sy_hi*sz_hi*vel_g(i  , j  , k  ,1:3)

      else

         i = floor((particles(p) % pos(1) - x0(1))*odx)
         j = floor((particles(p) % pos(2) - x0(2))*ody)
         k = floor((particles(p) % pos(3) - x0(3))*odz)

         velfp(1:3) = vel_g(i,j,k,1:3)

      endif

      ! Calculate drag coefficient, beta
      call des_drag_gp(slo, shi, p, particles(p) % vel, velfp, &
       ep_g(i,j,k), ro_g, mu_g, beta, i, j, k,             &
       particles(p) % radius,  particles(p) % volume,      &
       particles(p) % density, particles(p) % phase )

      particles(p) % drag(1) = beta

   end do

end subroutine calc_particle_beta_eb
