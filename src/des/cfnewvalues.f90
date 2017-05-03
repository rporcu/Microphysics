module CFNEWVALUES_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none
   private

   public  cfnewvalues
   public  cfnewvalues_aos



contains

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  Module name: CFNEWVALUES                                            !
   !                                                                      !
   !  Purpose: DES - Calculate the new values of particle velocity,       !
   !           position, angular velocity etc                             !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine CFNEWVALUES(max_pip, particle_state, pmass,&
        omoi, drag, des_pos_new, des_vel_new, omega_new, fc, tow, &
        des_acc_old, rot_acc_old)

      use discretelement, only: dtsolid
      use discretelement, only: normal_particle, exiting_particle
      use param, only: zero
      use constant, only: gravity

      implicit none

      integer     , intent(in   ) :: max_pip
      integer     , intent(in   ) :: particle_state(max_pip)
      real(c_real), intent(in   ) :: omoi(max_pip)
      real(c_real), intent(in   ) :: pmass(max_pip)
      real(c_real), intent(in   ) :: drag(max_pip,3)
      real(c_real), intent(inout) :: des_pos_new(max_pip,3)
      real(c_real), intent(inout) :: des_vel_new(max_pip,3)
      real(c_real), intent(inout) :: omega_new(max_pip,3)
      real(c_real), intent(inout) :: fc(max_pip,3)
      real(c_real), intent(inout) :: tow(max_pip,3)
      real(c_real), intent(inout) :: des_acc_old(max_pip,3)
      real(c_real), intent(inout) :: rot_acc_old(max_pip,3)

      !-----------------------------------------------
      ! Local Variables
      !-----------------------------------------------
      integer :: L
      real(c_real) :: DD(3)
      logical, save :: first_pass = .true.

      real(c_real) :: lVELo(3), lPOSo(3)

      !-----------------------------------------------

      ! Adams-Bashforth defaults to Euler for the first time step.
      if(first_pass) then
         first_pass = .false.
         do l =1, max_pip
            if(particle_state(l) == normal_particle .or.      &
                 particle_state(l) == exiting_particle) then
               des_acc_old(l,:) = (fc(l,:) + drag(l,:))/pmass(l) + gravity(:)
               rot_acc_old(l,:) = tow(l,:)
            endif
         enddo
      endif

      do l = 1, max_pip
         if(particle_state(l) == normal_particle .or.      &
              particle_state(l) == exiting_particle) then

            fc(l,:) = (fc(l,:) + drag(l,:))/pmass(l) + gravity(:)

            ! Advance particle position, velocity
            lvelo = des_vel_new(l,:)
            lposo = des_pos_new(l,:)

            ! Second-order Adams-Bashforth/Trapezoidal scheme
            des_vel_new(l,:) = lvelo(:) + 0.5d0*&
                 ( 3.d0*fc(l,:)-des_acc_old(l,:) )*dtsolid

            omega_new(l,:)   =  omega_new(l,:) + 0.5d0*&
                 ( 3.d0*tow(l,:)*omoi(l)-rot_acc_old(l,:) )*dtsolid

            dd(:) = 0.5d0*( lvelo(:)+des_vel_new(l,:) )*dtsolid

            des_pos_new(l,:) = lposo(:) + dd(:)
            des_acc_old(l,:) = fc(l,:)
            rot_acc_old(l,:) = tow(l,:)*omoi(l)

            ! Reset total contact force and torque
            fc(l,:) = zero
            tow(l,:) = zero
         endif
      enddo

      return
   end subroutine cfnewvalues



   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  Module name: CFNEWVALUES                                            !
   !                                                                      !
   !  Purpose: DES - Calculate the new values of particle velocity,       !
   !           position, angular velocity etc                             !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine cfnewvalues_aos( particles, fc, tow )

      use discretelement, only: dtsolid
      use discretelement, only: normal_particle, exiting_particle
      use param,          only: zero
      use constant,       only: gravity
      use particle_mod,   only: particle_t
      
      implicit none

      type(particle_t), intent(inout) :: particles(:)
      real(c_real),     intent(inout) :: fc(:,:), tow(:,:)
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

   
   end subroutine cfnewvalues_aos






   
end module cfnewvalues_module
