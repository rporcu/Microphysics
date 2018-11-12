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
subroutine calc_drag_particle( gp,  gplo,  gphi,  &
 &                            gp0, gp0lo, gp0hi,  &
 &                            vel,   ulo,   uhi,  &
 &                             np, particles, dx, &
 &                             x0 )  bind(C)

   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding,     only: c_int
   use particle_mod,      only: particle_t
   use param,             only: half, zero, one

   implicit none

   ! Array bounds
   integer(c_int), intent(in   )        ::   gplo(3),  gphi(3)
   integer(c_int), intent(in   )        ::  gp0lo(3), gp0hi(3)
   integer(c_int), intent(in   )        ::    ulo(3),   uhi(3)   

   ! Arrays
   real(rt),       intent(in   )        ::                           &
    & vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3),              &
    &  gp(gplo(1):gphi(1),gplo(2):gphi(2),gplo(3):gphi(3),3),        &
    & gp0(gp0lo(1):gp0hi(1),gp0lo(2):gp0hi(2),gp0lo(3):gp0hi(3),3)  

   ! Particles 
   integer(c_int),   intent(in   )    :: np
   type(particle_t), intent(inout)    :: particles(np)

   ! Grid 
   real(rt),         intent(in   )    :: dx(3)

   ! Coordinates of domain lower corner
   real(rt),         intent(in   )    :: x0(3)   


   ! Local variables
   !---------------------------------------------------------------------//
   ! Loop counters: Particle, fluid cell, neighbor cells
   integer  :: p, i, j, k
   real(rt) :: velfp(3), gradpg(3)
   real(rt) :: beta(np)
   real(rt) :: odx, ody, odz
   real(rt) :: lx, ly, lz
   real(rt) :: sx_lo, sy_lo, sz_lo
   real(rt) :: sx_hi, sy_hi, sz_hi
   real(rt) :: wx, wy, wz 

   odx = one/dx(1)
   ody = one/dx(2)
   odz = one/dx(3)

   ! Calculate the gas phase forces acting on each particle.
   do p = 1, np

      beta(p) = particles(p) % drag(1)

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

      gradpg(1:3) = sx_lo*sy_lo*sz_lo*(gp(i-1, j-1, k-1,1:3) + gp0(i-1, j-1, k-1,1:3)) + &
       &            sx_lo*sy_lo*sz_hi*(gp(i-1, j-1, k  ,1:3) + gp0(i-1, j-1, k  ,1:3)) + &
       &            sx_lo*sy_hi*sz_lo*(gp(i-1, j  , k-1,1:3) + gp0(i-1, j  , k-1,1:3)) + &
       &            sx_lo*sy_hi*sz_hi*(gp(i-1, j  , k  ,1:3) + gp0(i-1, j  , k  ,1:3)) + &
       &            sx_hi*sy_lo*sz_lo*(gp(i  , j-1, k-1,1:3) + gp0(i  , j-1, k-1,1:3)) + &
       &            sx_hi*sy_lo*sz_hi*(gp(i  , j-1, k  ,1:3) + gp0(i  , j-1, k  ,1:3)) + &
       &            sx_hi*sy_hi*sz_lo*(gp(i  , j  , k-1,1:3) + gp0(i  , j  , k-1,1:3)) + &
       &            sx_hi*sy_hi*sz_hi*(gp(i  , j  , k  ,1:3) + gp0(i  , j  , k  ,1:3))

      velfp(1:3) = sx_lo*sy_lo*sz_lo*vel(i-1, j-1, k-1,1:3) + &
       &           sx_lo*sy_lo*sz_hi*vel(i-1, j-1, k  ,1:3) + &
       &           sx_lo*sy_hi*sz_lo*vel(i-1, j  , k-1,1:3) + &
       &           sx_lo*sy_hi*sz_hi*vel(i-1, j  , k  ,1:3) + &
       &           sx_hi*sy_lo*sz_lo*vel(i  , j-1, k-1,1:3) + &
       &           sx_hi*sy_lo*sz_hi*vel(i  , j-1, k  ,1:3) + &
       &           sx_hi*sy_hi*sz_lo*vel(i  , j  , k-1,1:3) + &
       &           sx_hi*sy_hi*sz_hi*vel(i  , j  , k  ,1:3)          


      ! Particle drag calculation
      particles(p) % drag = beta(p)*(velfp - particles(p) % vel) - &
       &                    gradpg(:) * particles(p) % volume

   end do

end subroutine calc_drag_particle


subroutine calc_drag_particle_eb( gp,  gplo,   gphi,  &
 &                               gp0, gp0lo,  gp0hi,  &
 &                               vel,   ulo,    uhi,  &
 &                             flags,   flo,    fhi,  &
 &                                np, particles, dx,  &
 &                                x0 )  bind(C)

   use amrex_ebcellflag_module, only: is_covered_cell
   use amrex_fort_module,       only: rt => amrex_real
   use iso_c_binding,           only: c_int
   use particle_mod,            only: particle_t
   use param,                   only: half, zero, one

   implicit none

   ! Array bounds
   integer(c_int), intent(in   )        ::   gplo(3),  gphi(3)
   integer(c_int), intent(in   )        ::  gp0lo(3), gp0hi(3)
   integer(c_int), intent(in   )        ::    ulo(3),   uhi(3)   
   integer(c_int), intent(in   )        ::    flo(3),   fhi(3)

   ! Arrays
   real(rt),       intent(in   )        ::                           &
    & vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3),              &
    &  gp(gplo(1):gphi(1),gplo(2):gphi(2),gplo(3):gphi(3),3),        &
    & gp0(gp0lo(1):gp0hi(1),gp0lo(2):gp0hi(2),gp0lo(3):gp0hi(3),3)

   integer(c_int), intent(in   ) ::  &
    & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

   ! Particles 
   integer(c_int),   intent(in   )    :: np
   type(particle_t), intent(inout)    :: particles(np)

   ! Grid 
   real(rt),         intent(in   )    :: dx(3)

   ! Coordinates of domain lower corner
   real(rt),         intent(in   )    :: x0(3)   


   ! Local variables
   !---------------------------------------------------------------------//
   ! Loop counters: Particle, fluid cell, neighbor cells
   integer  :: p, i, j, k, ic, jc, kc
   real(rt) :: velfp(3), gradpg(3)
   real(rt) :: beta(np)
   real(rt) :: odx, ody, odz
   real(rt) :: lx, ly, lz
   real(rt) :: sx_lo, sy_lo, sz_lo
   real(rt) :: sx_hi, sy_hi, sz_hi
   real(rt) :: wx, wy, wz 


   odx = one/dx(1)
   ody = one/dx(2)
   odz = one/dx(3)

   ! Calculate the gas phase forces acting on each particle.
   do p = 1, np

      beta(p) = particles(p) % drag(1)


      ! Pick upper cell in the stencil
      lx = (particles(p) % pos(1) - x0(1))*odx + half
      ly = (particles(p) % pos(2) - x0(2))*ody + half
      lz = (particles(p) % pos(3) - x0(3))*odz + half

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      sx_hi = lx - i;  sx_lo = one - sx_hi
      sy_hi = ly - j;  sy_lo = one - sy_hi
      sz_hi = lz - k;  sz_lo = one - sz_hi      

      ! For now velocity is always interpolated via tri-linear interpolation,
      ! even near/at EB cells
      velfp(1:3) = sx_lo*sy_lo*sz_lo*vel(i-1, j-1, k-1,1:3) + &
       &           sx_lo*sy_lo*sz_hi*vel(i-1, j-1, k  ,1:3) + &
       &           sx_lo*sy_hi*sz_lo*vel(i-1, j  , k-1,1:3) + &
       &           sx_lo*sy_hi*sz_hi*vel(i-1, j  , k  ,1:3) + &
       &           sx_hi*sy_lo*sz_lo*vel(i  , j-1, k-1,1:3) + &
       &           sx_hi*sy_lo*sz_hi*vel(i  , j-1, k  ,1:3) + &
       &           sx_hi*sy_hi*sz_lo*vel(i  , j  , k-1,1:3) + &
       &           sx_hi*sy_hi*sz_hi*vel(i  , j  , k  ,1:3)

      !
      ! gradp is interpolated only if there are not covered cells
      ! in the stencil. If any covered cell is present in the stencil,
      ! we use the cell gradp
      !
      if ( ( .not. is_covered_cell(flags(i-1,j-1,k-1)) ) .and. &
       &   ( .not. is_covered_cell(flags(i-1,j-1,k  )) ) .and. &
       &   ( .not. is_covered_cell(flags(i-1,j  ,k-1)) ) .and. &
       &   ( .not. is_covered_cell(flags(i-1,j  ,k  )) ) .and. &
       &   ( .not. is_covered_cell(flags(i  ,j-1,k-1)) ) .and. &
       &   ( .not. is_covered_cell(flags(i  ,j-1,k  )) ) .and. &
       &   ( .not. is_covered_cell(flags(i  ,j  ,k-1)) ) .and. &
       &   ( .not. is_covered_cell(flags(i  ,j  ,k  )) ) ) then

         gradpg(1:3) = sx_lo*sy_lo*sz_lo*(gp(i-1, j-1, k-1,1:3) + gp0(i-1, j-1, k-1,1:3)) + &
          &            sx_lo*sy_lo*sz_hi*(gp(i-1, j-1, k  ,1:3) + gp0(i-1, j-1, k  ,1:3)) + &
          &            sx_lo*sy_hi*sz_lo*(gp(i-1, j  , k-1,1:3) + gp0(i-1, j  , k-1,1:3)) + &
          &            sx_lo*sy_hi*sz_hi*(gp(i-1, j  , k  ,1:3) + gp0(i-1, j  , k  ,1:3)) + &
          &            sx_hi*sy_lo*sz_lo*(gp(i  , j-1, k-1,1:3) + gp0(i  , j-1, k-1,1:3)) + &
          &            sx_hi*sy_lo*sz_hi*(gp(i  , j-1, k  ,1:3) + gp0(i  , j-1, k  ,1:3)) + &
          &            sx_hi*sy_hi*sz_lo*(gp(i  , j  , k-1,1:3) + gp0(i  , j  , k-1,1:3)) + &
          &            sx_hi*sy_hi*sz_hi*(gp(i  , j  , k  ,1:3) + gp0(i  , j  , k  ,1:3))


         velfp(1:3) = sx_lo*sy_lo*sz_lo*vel(i-1, j-1, k-1,1:3) + &
          &           sx_lo*sy_lo*sz_hi*vel(i-1, j-1, k  ,1:3) + &
          &           sx_lo*sy_hi*sz_lo*vel(i-1, j  , k-1,1:3) + &
          &           sx_lo*sy_hi*sz_hi*vel(i-1, j  , k  ,1:3) + &
          &           sx_hi*sy_lo*sz_lo*vel(i  , j-1, k-1,1:3) + &
          &           sx_hi*sy_lo*sz_hi*vel(i  , j-1, k  ,1:3) + &
          &           sx_hi*sy_hi*sz_lo*vel(i  , j  , k-1,1:3) + &
          &           sx_hi*sy_hi*sz_hi*vel(i  , j  , k  ,1:3)

      else

         ic  = floor((particles(p) % pos(1) - x0(1))*odx)
         jc  = floor((particles(p) % pos(2) - x0(2))*ody)
         kc  = floor((particles(p) % pos(3) - x0(3))*odz)

         gradpg(:) = gp(ic,jc,kc,:) + gp0(ic,jc,kc,:)
         velfp(:)  = vel(ic,jc,kc,1:3)

      end if

      ! Particle drag calculation
      particles(p) % drag = beta(p)*(velfp - particles(p) % vel) - &
       &                gradpg(:) * particles(p) % volume

   end do

end subroutine calc_drag_particle_eb
