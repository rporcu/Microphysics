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

   subroutine call_usr2_des( np, particles ) &
        bind(c, name="call_usr2_des")

      use run,            only: call_usr
      use particle_mod,   only: particle_t

      integer(c_int),   intent(in)    :: np
      type(particle_t), intent(inout) :: particles(np)

      if ( call_usr ) call usr2_des(np, particles)

   end subroutine call_usr2_des

   subroutine des_time_loop_ops ( nrp, rparticles, ngp, gparticles, &
        & subdt, dx, dy, dz, xlength, ylength, zlength, nstep, ncoll)  &
        bind(C, name="des_time_loop_ops")

      use particle_mod
      use calc_collision_wall,     only: calc_dem_force_with_wall_stl
      use calc_force_dem_module,   only: calc_force_dem, calc_force_dem_nl
      use output_manager_module,   only: output_manager
      use run,                     only: call_usr

      integer(c_int),   intent(in   )     :: nrp, ngp
      real(c_real),     intent(in   )     :: subdt, dx, dy, dz
      real(c_real),     intent(in   )     :: xlength, ylength, zlength
      type(particle_t), intent(inout)     :: rparticles(nrp), gparticles(ngp)
      integer(c_int),   intent(in   )     :: nstep
      integer(c_int),   intent(inout)     :: ncoll

      real(c_real)                        :: tow(nrp+ngp,3), fc(nrp+ngp,3)
      type(particle_t)                    :: particles(nrp+ngp)

      tow  = 0
      fc   = 0

      particles(    1:nrp) = rparticles
      particles(nrp+1:   ) = gparticles

      ! calculate forces from particle-wall collisions
      call calc_dem_force_with_wall_stl ( particles, fc, tow, &
                                          xlength, ylength, zlength, subdt )

      ! calculate forces from particle-particle collisions
      call calc_force_dem ( particles, fc, tow, subdt, nstep, ncoll )

      rparticles = particles(    1: nrp)
      gparticles = particles(nrp+1:    )

      ! call user functions.
      if ( call_usr ) call usr1_des

      ! update position and velocities
      call des_euler_update ( rparticles, gparticles, fc, tow, subdt )

   end subroutine des_time_loop_ops

   subroutine des_time_loop_ops_nl ( nrp, rparticles, ngp, gparticles, size_nl, nbor_list, &
        & subdt, dx, dy, dz, xlength, ylength, zlength, nstep, ncoll )  &
        bind(C, name="des_time_loop_ops_nl")

      use particle_mod
      use calc_collision_wall,     only: calc_dem_force_with_wall_stl
      use calc_force_dem_module,   only: calc_force_dem_nl
      use output_manager_module,   only: output_manager
      use run,                     only: call_usr

      integer(c_int),   intent(in   )     :: nrp, ngp, size_nl
      real(c_real),     intent(in   )     :: subdt, dx, dy, dz
      real(c_real),     intent(in   )     :: xlength, ylength, zlength
      type(particle_t), intent(inout)     :: rparticles(nrp), gparticles(ngp)
      integer(c_int),   intent(in   )     :: nbor_list(size_nl)
      integer(c_int),   intent(in   )     :: nstep
      integer(c_int),   intent(inout)     :: ncoll

      type(particle_t), allocatable       :: particles(:)

      real(c_real), allocatable :: tow(:,:), fc(:,:)

      allocate(tow(nrp+ngp,3))
      allocate( fc(nrp+ngp,3))
      allocate(particles(nrp+ngp))

      tow  = 0
      fc   = 0

      particles(    1:nrp) = rparticles
      particles(nrp+1:   ) = gparticles

      ! calculate forces from particle-wall collisions
      call calc_dem_force_with_wall_stl ( particles, fc, tow, &
                                          xlength, ylength, zlength, subdt )

      ! calculate forces from particle-particle collisions
      call calc_force_dem_nl ( particles, nbor_list, size_nl, fc, tow, subdt, nstep, ncoll )

      rparticles = particles(    1:nrp)
      gparticles = particles(nrp+1:   )

      ! call user functions.
      if ( call_usr ) call usr1_des

      ! update position and velocities
      call des_euler_update ( rparticles, gparticles, fc, tow, subdt )

      deallocate(tow, fc, particles)

   end subroutine des_time_loop_ops_nl

   subroutine des_euler_update ( particles, grid_nbors, fc, tow, dt )

      use constant,       only: gravity
      use particle_mod,   only: particle_t

      type(particle_t), intent(inout)  :: particles(:)
      type(particle_t), intent(inout)  :: grid_nbors(:)
      real(c_real),     intent(inout)  :: fc(:,:), tow(:,:)
      real(c_real),     intent(in   )  :: dt
      integer                          :: p,np,ng

      np = size (particles)
      ng = size (grid_nbors)

      do p = 1, np

            particles(p) % vel     = particles(p) % vel   + dt * &
                ( ( fc(p,:) +  particles(p) % drag ) / particles(p) % mass + gravity )
            particles(p) % pos     = particles(p) % pos   + dt * particles(p) % vel
            particles(p) % omega   = particles(p) % omega + dt * tow(p,:) * particles(p) % omoi

      end do

      do p = 1, ng

            grid_nbors(p) % vel     = grid_nbors(p) % vel   + dt * &
                ( ( fc(np+p,:) +  grid_nbors(p) % drag ) / grid_nbors(p) % mass + gravity )
            grid_nbors(p) % pos     = grid_nbors(p) % pos   + dt * grid_nbors(p) % vel
            grid_nbors(p) % omega   = grid_nbors(p) % omega + dt * tow(np+p,:) * grid_nbors(p) % omoi

      end do

   end subroutine des_euler_update

end module des_time_march_module
