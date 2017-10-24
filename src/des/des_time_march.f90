module des_time_march_module

   use amrex_fort_module, only: c_real => amrex_real
   use iso_c_binding ,    only: c_int

   implicit none

   private

   ! Define interface for external functions so there will be no
   ! warning at compile time
   interface

      subroutine usr0_des()
      end subroutine usr0_des

      subroutine usr1_des()
      end subroutine usr1_des

   end interface


contains


   function des_is_continuum_coupled () result( is_coupled ) &
        bind(C, name="des_continuum_coupled")

      use discretelement, only: des_continuum_coupled

      integer  :: is_coupled

      is_coupled = 0

      if ( des_continuum_coupled ) is_coupled = 1

   end function des_is_continuum_coupled

   subroutine des_init_time_loop ( tstart, dt, nsubsteps, subdt ) &
        bind(C, name="des_init_time_loop")

      use discretelement,  only: dtsolid, des_continuum_coupled
      use run,             only: tstop

      real(c_real),   intent(in   ) :: tstart, dt
      integer(c_int), intent(  out) :: nsubsteps
      real(c_real),   intent(  out) :: subdt

      ! Initialize time stepping variables for
      ! coupled gas/solids simulations.
      if ( des_continuum_coupled ) then

         if ( dt >= dtsolid ) then
            nsubsteps = ceiling ( real ( dt / dtsolid ) )
            subdt     =  dt / nsubsteps
         else
            nsubsteps = 1
            subdt     = dtsolid
         end if

         ! Initialize time stepping variable for pure granular simulations.
      else

         nsubsteps = ceiling ( real ( (tstop - tstart) / dtsolid ) )
         subdt     = ( tstop - tstart ) / nsubsteps

      end if

   end subroutine des_init_time_loop

   subroutine call_usr3_des( np, particles ) &
        bind(c, name="call_usr3_des")

      use run,            only: call_usr
      use particle_mod,   only: particle_t

      integer(c_int),   intent(in)    :: np
      type(particle_t), intent(inout) :: particles(np)

      if ( call_usr ) call usr3_des(np, particles)

   end subroutine call_usr3_des

   subroutine des_time_loop_ops_nl ( n_do, nrp, rparticles, ngp, gparticles, size_nl, &
        nbor_list, subdt, dx, xlength, ylength, zlength, ncoll, stime, &
        flag, fglo, fghi, bcent, blo, bhi, apx, axlo, axhi, apy, aylo, ayhi, &
        apz, azlo, azhi)  &
          bind(C, name="des_time_loop_ops_nl")

     use particle_mod
     use calc_collision_wall,     only: calc_dem_force_with_wall
     use calc_force_dem_module,   only: calc_force_dem_nl
     use discretelement       ,   only: des_continuum_coupled
     use output_manager_module,   only: output_manager
     use run,                     only: call_usr

     integer(c_int),   intent(in   )     :: n_do, nrp, ngp, size_nl
     real(c_real),     intent(in   )     :: subdt, dx(3)
     real(c_real),     intent(in   )     :: xlength, ylength, zlength
     type(particle_t), intent(inout)     :: rparticles(nrp), gparticles(ngp)
     integer(c_int),   intent(in   )     :: nbor_list(size_nl)
     integer(c_int),   intent(inout)     :: ncoll
     real(c_real),     intent(inout)     :: stime

     integer, dimension(3), intent(in) :: axlo,axhi
     integer, dimension(3), intent(in) :: aylo,ayhi
     integer, dimension(3), intent(in) :: azlo,azhi
     integer, dimension(3), intent(in) :: fglo,fghi
     integer, dimension(3), intent(in) :: blo,bhi

      integer,      intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
      real(c_real), intent(in) :: bcent  (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
      real(c_real), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
      real(c_real), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
      real(c_real), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

      type(particle_t), allocatable       :: particles(:)

      integer(c_int)            :: nstep
      real(c_real), allocatable :: tow(:,:), fc(:,:)

      allocate(tow(nrp+ngp,3))
      allocate( fc(nrp+ngp,3))
      allocate(particles(nrp+ngp))

      do nstep = 1, n_do

         tow  = 0
         fc   = 0

         particles(    1:nrp) = rparticles
         particles(nrp+1:   ) = gparticles

         ! calculate forces from particle-wall collisions
         call calc_dem_force_with_wall (particles, nrp+ngp, nrp, fc, tow, subdt, &
              flag, fglo, fghi, bcent, blo, bhi, apx, axlo, axhi, apy, aylo, ayhi, &
              apz, azlo, azhi, dx)

         ! calculate forces from particle-particle collisions
         call calc_force_dem_nl ( particles, nrp, nbor_list, size_nl, fc, tow, subdt, ncoll )

         ! call user functions.
         if ( call_usr ) call usr1_des

         ! update position and velocities
         call des_euler_update ( rparticles, nrp+ngp, nrp, fc, tow, subdt )

         if ( call_usr ) call usr2_des(nrp, rparticles );

         if ( .not.des_continuum_coupled ) call output_manager(nrp+ngp,  &
              stime, subdt, xlength, ylength, zlength, nstep, particles, 0)

         stime = stime + subdt

      end do

      deallocate(tow, fc, particles)

   end subroutine des_time_loop_ops_nl

   subroutine des_euler_update ( particles, np, nrp, fc, tow, dt )

      use constant,       only: gravity
      use particle_mod,   only: particle_t

      type(particle_t), intent(inout)  :: particles(np)
      real(c_real),     intent(inout)  :: fc(np,3), tow(np,3)
      real(c_real),     intent(in   )  :: dt
      integer,          intent(in   )  :: np, nrp
      integer                          :: p

      do p = 1, nrp

            particles(p) % vel     = particles(p) % vel   + dt * &
                ( ( fc(p,:) +  particles(p) % drag ) / particles(p) % mass + gravity )
            particles(p) % pos     = particles(p) % pos   + dt * particles(p) % vel
            particles(p) % omega   = particles(p) % omega + dt * tow(p,:) * particles(p) % omoi

      end do

   end subroutine des_euler_update

end module des_time_march_module
