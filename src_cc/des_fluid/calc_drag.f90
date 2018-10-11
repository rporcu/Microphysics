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
     &                        gp0, gp0lo, gp0hi,  &
     &                        vel,   ulo,   uhi,  &
     &                        xsl,   ysl,         &
     &                        zsl,   slo,   shi,  &
     &                        np, particles, dx,  &
     &                        x0,  interp_type )  &
     &                        bind(C, name="calc_drag_particle")

   use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
   use amrex_error_module,              only: amrex_abort
   use amrex_fort_module,               only: rt => amrex_real
   use iso_c_binding,                   only: c_int
   use particle_mod,                    only: particle_t
   use param,                           only: half, zero, one

   implicit none

   ! Array bounds
   integer(c_int), intent(in   )        ::   gplo(3),  gphi(3)
   integer(c_int), intent(in   )        ::  gp0lo(3), gp0hi(3)
   integer(c_int), intent(in   )        ::    ulo(3),   uhi(3)   
   integer(c_int), intent(in   )        ::    slo(3),   shi(3)

   ! Arrays
   real(rt),       intent(in   )        ::                               &
        & vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3),              &
        &  gp(gplo(1):gphi(1),gplo(2):gphi(2),gplo(3):gphi(3),3),        &
        & gp0(gp0lo(1):gp0hi(1),gp0lo(2):gp0hi(2),gp0lo(3):gp0hi(3),3),  &
        & xsl(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3),              &
        & ysl(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3),              &
        & zsl(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)


   ! Particles 
   integer(c_int),   intent(in   )    :: np
   type(particle_t), intent(inout)    :: particles(np)
   
   ! Grid 
   real(rt),         intent(in   )    :: dx(3)

   ! Coordinates of domain lower corner
   real(rt),         intent(in   )    :: x0(3)   

   ! Interpolation type
   integer(c_int),   intent(in   ), value :: interp_type



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

   if (amrex_pd_ioprocessor()) then
      if (interp_type==0) then 
         write(*,*) 'WARNING: No interpolation of particle drag force'
      else if (interp_type==1) then
         write(*,*) 'WARNING: Standard interpolation of particle drag force'
      else if (interp_type==2) then
         write(*,*) 'WARNING: Slope-interpolation of particle drag force'
      else
         call amrex_abort("calc_drag_particle(): &
              & argument 'interp_type' has wrong value. Must be set to 0,1, or 2")
      end if
   end if

   odx = one/dx(1)
   ody = one/dx(2)
   odz = one/dx(3)

   ! Calculate the gas phase forces acting on each particle.
   do p = 1, np

      beta(p) = particles(p) % drag(1)
      
      if ( interp_type == 0 ) then   !No interpolation

         i  = floor((particles(p) % pos(1) - x0(1))*odx)
         j  = floor((particles(p) % pos(2) - x0(2))*ody)
         k  = floor((particles(p) % pos(3) - x0(3))*odz)

         velfp(:)  = vel(i,j,k,:)
         gradpg(:) = gp(i,j,k,:) + gp0(i,j,k,:)

      else                           !Standard and slope-based interpolation

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
              &        sx_lo*sy_lo*sz_hi*(gp(i-1, j-1, k  ,1:3) + gp0(i-1, j-1, k  ,1:3)) + &
              &        sx_lo*sy_hi*sz_lo*(gp(i-1, j  , k-1,1:3) + gp0(i-1, j  , k-1,1:3)) + &
              &        sx_lo*sy_hi*sz_hi*(gp(i-1, j  , k  ,1:3) + gp0(i-1, j  , k  ,1:3)) + &
              &        sx_hi*sy_lo*sz_lo*(gp(i  , j-1, k-1,1:3) + gp0(i  , j-1, k-1,1:3)) + &
              &        sx_hi*sy_lo*sz_hi*(gp(i  , j-1, k  ,1:3) + gp0(i  , j-1, k  ,1:3)) + &
              &        sx_hi*sy_hi*sz_lo*(gp(i  , j  , k-1,1:3) + gp0(i  , j  , k-1,1:3)) + &
              &        sx_hi*sy_hi*sz_hi*(gp(i  , j  , k  ,1:3) + gp0(i  , j  , k  ,1:3))

         ! Fluid velocity at the particle's position.
         if ( interp_type == 2 ) then
            
            velfp(1:3) = sx_lo*sy_lo*sz_lo*vel(i-1, j-1, k-1,1:3) + &
                 &       sx_lo*sy_lo*sz_hi*vel(i-1, j-1, k  ,1:3) + &
                 &       sx_lo*sy_hi*sz_lo*vel(i-1, j  , k-1,1:3) + &
                 &       sx_lo*sy_hi*sz_hi*vel(i-1, j  , k  ,1:3) + &
                 &       sx_hi*sy_lo*sz_lo*vel(i  , j-1, k-1,1:3) + &
                 &       sx_hi*sy_lo*sz_hi*vel(i  , j-1, k  ,1:3) + &
                 &       sx_hi*sy_hi*sz_lo*vel(i  , j  , k-1,1:3) + &
                 &       sx_hi*sy_hi*sz_hi*vel(i  , j  , k  ,1:3)
            
         else

            lx = (particles(p) % pos(1) - x0(1))*odx 
            ly = (particles(p) % pos(2) - x0(2))*ody 
            lz = (particles(p) % pos(3) - x0(3))*odz 
            
            i = floor(lx)
            j = floor(ly)
            k = floor(lz)

            wx = lx - i - half
            wy = ly - j - half
            wz = lz - k - half

            velfp(1:3) = vel(i,j,k,1:3)    &
                 &     + xsl(i,j,k,1:3)*wx &
                 &     + ysl(i,j,k,1:3)*wy &
                 &     + zsl(i,j,k,1:3)*wz

         end if
         
      end if

      ! Particle drag calculation
      particles(p) % drag = beta(p)*(velfp - particles(p) % vel) - &
           &                gradpg(:) * particles(p) % volume

   end do

end subroutine calc_drag_particle
