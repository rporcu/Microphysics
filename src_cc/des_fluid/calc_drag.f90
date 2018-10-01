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
                               xslopes, yslopes,&
                               zslopes, slo, shi, &
                               np, particles, dx, dy, dz, use_slopes) &
     bind(C, name="calc_drag_particle")

   use amrex_paralleldescriptor_module, only : amrex_pd_ioprocessor
   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int
   use des_drag_gp_module, only: des_drag_gp
   use particle_mod, only: particle_t
   use param,        only: half, zero, one

   implicit none

   integer(c_int), intent(in   ) ::    ulo(3), uhi(3)
   integer(c_int), intent(in   ) ::   gplo(3), gphi(3)
   integer(c_int), intent(in   ) ::  gp0lo(3), gp0hi(3)
   integer(c_int), intent(in   ) ::    slo(3), shi(3)
   integer(c_int), intent(in   ) ::  np
   integer(c_int), intent(in   ), value :: use_slopes

   real(rt)  , intent(in   ) :: &
        vel_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3), &
        gp(gplo(1):gphi(1),gplo(2):gphi(2),gplo(3):gphi(3),3),  &
        gp0(gp0lo(1):gp0hi(1),gp0lo(2):gp0hi(2),gp0lo(3):gp0hi(3),3), &
        xslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
        yslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
        zslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

   real(rt),     intent(in   ) :: dx, dy, dz
   type(particle_t), intent(inout) :: particles(np)


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
   real(rt) :: plo(3)

   logical, parameter :: no_interpolation = .true.

   if(no_interpolation .and. amrex_pd_ioprocessor()) &
        write(*,*) 'WARNING: No interpolation of particle drag force'

   plo = zero ! HACK -- This should get passed into routine

   odx = one/dx
   ody = one/dy
   odz = one/dz

   do p = 1, np
      beta(p) = particles(p) % drag(1)
   end do

   ! Calculate the gas phase forces acting on each particle.

   do p = 1, np

      lx = (particles(p) % pos(1) - plo(1))*odx + half
      ly = (particles(p) % pos(2) - plo(2))*ody + half
      lz = (particles(p) % pos(3) - plo(3))*odz + half

      i = floor(lx)
      j = floor(ly)
      k = floor(lz)

      sx_hi = lx - i;  sx_lo = one - sx_hi
      sy_hi = ly - j;  sy_lo = one - sy_hi
      sz_hi = lz - k;  sz_lo = one - sz_hi

      ! Fluid velocity at the particle's position.
      if (use_slopes == 0 ) then
         velfp(1:3) = sx_lo*sy_lo*sz_lo*vel_g(i-1, j-1, k-1,1:3) + &
             &        sx_lo*sy_lo*sz_hi*vel_g(i-1, j-1, k  ,1:3) + &
             &        sx_lo*sy_hi*sz_lo*vel_g(i-1, j  , k-1,1:3) + &
             &        sx_lo*sy_hi*sz_hi*vel_g(i-1, j  , k  ,1:3) + &
             &        sx_hi*sy_lo*sz_lo*vel_g(i  , j-1, k-1,1:3) + &
             &        sx_hi*sy_lo*sz_hi*vel_g(i  , j-1, k  ,1:3) + &
             &        sx_hi*sy_hi*sz_lo*vel_g(i  , j  , k-1,1:3) + &
             &        sx_hi*sy_hi*sz_hi*vel_g(i  , j  , k  ,1:3)

      else

         block
            integer  :: ic, jc, kc  ! Indeces of the cell containing the particle
            real(rt) :: xc(3)       ! Coordinates of cell center
            real(rt) :: wx, wy, wz  ! Interpolation weights

            ic = floor((particles(p) % pos(1) - plo(1))*odx)
            jc = floor((particles(p) % pos(2) - plo(2))*ody)
            kc = floor((particles(p) % pos(3) - plo(3))*odz)

            xc = plo + [ic,jc,kc]*dx + half*dx ! Assuming that plo is the origin

            wx = (particles(p) % pos(1) - xc(1)) * odx
            wy = (particles(p) % pos(2) - xc(2)) * ody
            wz = (particles(p) % pos(3) - xc(3)) * odz

            velfp(1:3) = vel_g(ic,jc,kc,1:3)     &
                 &     + xslopes(i,j,k,1:3)*wx   &
                 &     + yslopes(i,j,k,1:3)*wy   &
                 &     + zslopes(i,j,k,1:3)*wz

         end block

      end if

      gradpg(1:3) = sx_lo*sy_lo*sz_lo*(gp(i-1, j-1, k-1,1:3) + gp0(i-1, j-1, k-1,1:3)) + &
                    sx_lo*sy_lo*sz_hi*(gp(i-1, j-1, k  ,1:3) + gp0(i-1, j-1, k  ,1:3)) + &
                    sx_lo*sy_hi*sz_lo*(gp(i-1, j  , k-1,1:3) + gp0(i-1, j  , k-1,1:3)) + &
                    sx_lo*sy_hi*sz_hi*(gp(i-1, j  , k  ,1:3) + gp0(i-1, j  , k  ,1:3)) + &
                    sx_hi*sy_lo*sz_lo*(gp(i  , j-1, k-1,1:3) + gp0(i  , j-1, k-1,1:3)) + &
                    sx_hi*sy_lo*sz_hi*(gp(i  , j-1, k  ,1:3) + gp0(i  , j-1, k  ,1:3)) + &
                    sx_hi*sy_hi*sz_lo*(gp(i  , j  , k-1,1:3) + gp0(i  , j  , k-1,1:3)) + &
                    sx_hi*sy_hi*sz_hi*(gp(i  , j  , k  ,1:3) + gp0(i  , j  , k  ,1:3))


      ! Use cell center values
      if(no_interpolation) then

         i = floor((particles(p) % pos(1) - plo(1))*odx)
         j = floor((particles(p) % pos(2) - plo(2))*ody)
         k = floor((particles(p) % pos(3) - plo(3))*odz)

         velfp(:) = vel_g(i,j,k,:)

         gradpg(1:3) = gp(i,j,k,:) + gp0(i,j,k,:)
      endif


      particles(p) % drag = beta(p)*(velfp - particles(p) % vel) - &
                            gradpg(:) * particles(p) % volume

   enddo

end subroutine calc_drag_particle
