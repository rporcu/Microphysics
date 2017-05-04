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
subroutine USR2_DES( np, particles )

   use amrex_fort_module, only: c_real => amrex_real
   use particle_mod,      only: particle_t
   
   implicit none

   integer,          intent(in   ) :: np
   type(particle_t), intent(inout) :: particles(np)
   integer                         :: p

   ! Move particles 63-93 below particles 32-62 to fake a wall.
   do p = 63, np
      particles(p) % vel(:)   = 0.0d0
      particles(p) % pos(2)   = 0.0475d0
      particles(p) % pos(1)   = particles(p-31) % pos(1)
      particles(p) % pos(3)   = particles(p-31) % pos(3)
      particles(p) % omega(:) = 0.0d0
   end do
   
end subroutine USR2_DES
