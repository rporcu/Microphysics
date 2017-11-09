!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine usr2_des( np, particles )

   use amrex_fort_module, only: c_real => amrex_real
   use usr,               only: bounce_count, max_bounce, yvelo, yposo, max_height
   use particle_mod,      only: particle_t

   implicit none

   integer,          intent(in   ) :: np
   type(particle_t), intent(in   ) :: particles(np)

   real(c_real) :: pPos, pVel

   ! Check if particle reached the peak of the bounce.

   pVel = minval(particles(1)%vel(1:2))
   pPos = sqrt(particles(1)%pos(1)**2 + particles(1)%pos(2)**2) - 0.5d0

   if ( pVel < 0.0 ) then
      if(yvelo > 0.0) then
         if(bounce_count < max_bounce) then
            bounce_count = bounce_count + 1
            max_height(bounce_count) = max(yPoso, pPos)
         endif
      endif
   endif

   yVelo = pVel
   yPoso = pPos

   RETURN
END SUBROUTINE USR2_DES
