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
   use interpolation_m,   only: trilinear_interp
   
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
   real(rt) :: pbeta
   real(rt) :: odx, ody, odz

   odx = one/dx(1)
   ody = one/dx(2)
   odz = one/dx(3)

   ! Calculate the gas phase forces acting on each particle.
   do p = 1, np

      associate( ppos => particles(p) % pos, pvel => particles(p) % vel)

         pbeta = particles(p) % drag(1)

         ! Pick upper cell in the stencil
         i = floor((ppos(1) - x0(1))*odx + half)
         j = floor((ppos(2) - x0(2))*ody + half)
         k = floor((ppos(3) - x0(3))*odz + half)
         
         velfp(:)  = trilinear_interp(vel, ulo, uhi, 3, ppos, x0, dx)           
         gradpg(:) = trilinear_interp(gp, gp0, gplo, gphi, 3, ppos, x0, dx )
            

         ! Particle drag calculation
         particles(p) % drag = pbeta*(velfp - pvel) - &
          &                gradpg(:) * particles(p) % volume

      end associate

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
   use interpolation_m,         only: trilinear_interp

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
   real(rt) :: pbeta
   real(rt) :: odx, ody, odz


   odx = one/dx(1)
   ody = one/dx(2)
   odz = one/dx(3)

   ! Calculate the gas phase forces acting on each particle.
   ! The velocity field is assumed to be reconstructed in a narrow band
   ! surrouding the walls so that we can interpolate with a standard trilinear
   ! interpolation
   do p = 1, np

      associate( ppos => particles(p) % pos, pvel => particles(p) % vel)
         
         pbeta = particles(p) % drag(1)

         ! Pick upper cell in the stencil
         i = floor((ppos(1) - x0(1))*odx + half)
         j = floor((ppos(2) - x0(2))*ody + half)
         k = floor((ppos(3) - x0(3))*odz + half)
         
         velfp(:) = trilinear_interp(vel, ulo, uhi, 3, ppos, x0, dx)

         !
         ! gradp is interpolated only if there are not covered cells
         ! in the stencil. If any covered cell is present in the stencil,
         ! we use the cell gradp
         !
         if ( any( is_covered_cell(flags(i-1:i,j-1:j,k-1:k) )) ) then

            ic  = floor((particles(p) % pos(1) - x0(1))*odx)
            jc  = floor((particles(p) % pos(2) - x0(2))*ody)
            kc  = floor((particles(p) % pos(3) - x0(3))*odz)

            gradpg(:) = gp(ic,jc,kc,:) + gp0(ic,jc,kc,:)
            
         else
            
            gradpg(:) = trilinear_interp(gp, gp0, gplo, gphi, 3, ppos, x0, dx )
            
         end if

         ! Particle drag calculation
         particles(p) % drag = pbeta*(velfp - pvel) - &
          &                gradpg(:) * particles(p) % volume

      end associate

   end do

end subroutine calc_drag_particle_eb
