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

   subroutine des_init_time_loop ( tstart, dt, nsubsteps, subdt, subdt_io) &
        bind(C, name="des_init_time_loop")

      use discretelement,  only: dtsolid
      use run,             only: glob_subdt_io => subdt_io, des_tstart, des_dt

      real(c_real),   intent(in   ) :: tstart, dt
      integer(c_int), intent(in   ) :: subdt_io
      integer(c_int), intent(  out) :: nsubsteps
      real(c_real),   intent(  out) :: subdt

      ! set the global subdt_io (in run module) to toggle sub-dt I/O
      glob_subdt_io = .true.
      if ( subdt_io == 0 )  glob_subdt_io = .false.

      ! update the global des_tstart, and des_dt (in run module) corresponding to this 
      ! des run: this enables usr[2,3]_des to know the time
      des_tstart = tstart
      des_dt     = dt

      ! Initialize time stepping variables 
      if ( dt >= dtsolid ) then
          nsubsteps = ceiling ( real ( dt / dtsolid ) )
          subdt     =  dt / nsubsteps
      else
          nsubsteps = 1
          subdt     = dt
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

   subroutine des_time_loop ( np, particles, &
                              nf , tow, fc, subdt, &
                              xlength, ylength, zlength, stime, nstep) & 
          bind(C, name="des_time_loop")

      use particle_mod
      use constant                       , only: gravity
      use output_manager_module          , only: output_manager
      use run                            , only: call_usr, subdt_io
      use bc                             , only: BC_shaker_A, BC_shaker_F
 
      integer(c_int),   intent(in   )     :: nf, np
      real(c_real),     intent(in   )     :: subdt, xlength, ylength, zlength
      type(particle_t), intent(inout)     :: particles(np)
      real(c_real),     intent(inout)     :: tow(nf,3)
      real(c_real),     intent(inout)     :: fc(nf,3)
      real(c_real),     intent(inout)     :: stime
      integer(c_int),   intent(in   )     :: nstep

      integer :: p
      
      ! call user functions.
      if ( call_usr ) call usr1_des

      if ( subdt_io ) call output_manager(np,  &
           stime, subdt, xlength, ylength, zlength, nstep, particles, 0)

      ! Update the time for the purpose of printing
      stime = stime + subdt

      ! Update position and velocities
      do p = 1, np
            ! forces acting on particles => acceleration
            !   |---> fc: particle-particle collisions and particle-wall collisions
            !   |---> drag: drag force (particle-fluid coupling)
            !   `---> gravity
            particles(p) % vel     = particles(p) % vel   + subdt * &
                        ( ( fc(p,:) +  particles(p) % drag ) / particles(p) % mass + gravity )

            ! in case of shaking, the simulation is in the co-shaken frame => non-interial frame => fictious forces
            ! fictitions force due to frame tranlation = - m * x_frame'' (without rotation)
            if ( ( BC_shaker_F .gt. 0.d0) .and. any( BC_shaker_A .gt. 0.d0 ) ) then
                particles(p) % vel = particles(p) % vel    - subdt * &
                        BC_shaker_A * BC_shaker_F**2 * sin(BC_shaker_F * stime)
            end if

            ! tranlate particle (Euler-step using updated velocity)
            particles(p) % pos     = particles(p) % pos   + subdt * particles(p) % vel

            ! update particle rotational velocity 
            particles(p) % omega   = particles(p) % omega + subdt * tow(p,:) * particles(p) % omoi

            ! note: particle omega should also couple to particle vel
      end do

      if ( call_usr ) call usr2_des(np, particles );

   end subroutine des_time_loop

end module des_time_march_module
