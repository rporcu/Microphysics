module CFNEWVALUES_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none
   private


   public  cfnewvalues



contains

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  Module name: CFNEWVALUES                                            !
   !                                                                      !
   !  Purpose: DES - Calculate the new values of particle velocity,       !
   !           position, angular velocity etc                             !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine cfnewvalues( particles, fc, tow, dtsolid )

      use discretelement, only: normal_particle, exiting_particle
      use param,          only: zero
      use constant,       only: gravity
      use particle_mod,   only: particle_t
      
      implicit none

      type(particle_t), intent(inout) :: particles(:)
      real(c_real),     intent(inout) :: fc(:,:), tow(:,:)
      real(c_real),     intent(in)    :: dtsolid
      integer                         :: L
      real(c_real)                    :: DD(3), lVELo(3), lPOSo(3)
      logical, save                   :: first_pass = .true.

      ! Adams-Bashforth defaults to Euler for the first time step.
      if (first_pass) then
         
         first_pass = .false.
         
         do l = 1, size(particles)
            if ( particles(l) % state == normal_particle .or. &
                 particles(l) % state == exiting_particle ) then
               particles(l) % acc   = ( fc(l,:) + particles(l) % drag ) / particles(l) % mass  +  gravity(:)
               particles(l) % alpha = tow(l,:)
            end if
         end do
         
      end if

      do l = 1, size(particles)

         if ( particles(l) % state == normal_particle .or. &
              particles(l) % state == exiting_particle ) then

            fc(l,:) = ( fc(l,:) + particles(l) % drag )/ particles(l) % mass + gravity(:)

            ! Advance particle position, velocity
            lvelo = particles(l) % vel
            lposo = particles(l) % pos

            ! Second-order Adams-Bashforth/Trapezoidal scheme
            particles(l) % vel   = lvelo(:) + 0.5d0 * &
                 & ( 3.d0 * fc(l,:) - particles(l) % acc ) * dtsolid

            particles(l) % omega = particles(l) % omega + 0.5d0 * &
                 ( 3.d0 * tow(l,:) * particles(l) % omoi - particles(l) % alpha ) * dtsolid

            dd(:) = 0.5d0*( lvelo(:) + particles(l) % vel ) * dtsolid

            particles(l) % pos    = lposo(:) + dd(:)
            particles(l) % acc   = fc(l,:)
            particles(l) % alpha = tow(l,:) * particles(l) % omoi

            ! Reset total contact force and torque
            fc(l,:)  = zero
            tow(l,:) = zero
            
         end if
         
      end do

   
   end subroutine cfnewvalues






   
end module cfnewvalues_module
