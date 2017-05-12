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
        bind(C, name="mfix_des_continuum_coupled")

      use discretelement, only: des_continuum_coupled
      
      integer  :: is_coupled

      is_coupled = 0

      if ( des_continuum_coupled ) is_coupled = 1

   end function des_is_continuum_coupled


   subroutine des_init_time_loop ( tstart, dt, nsubsteps, subdt ) &
        bind(C, name="mfix_des_init_time_loop")

      use discretelement,  only: dtsolid, des_continuum_coupled
      use run,             only: tstop
      use param,           only: zero
      
      real(c_real),   intent(in   ) :: tstart, dt
      integer(c_int), intent(  out) :: nsubsteps
      real(c_real),   intent(  out) :: subdt
     
      ! Initialize time stepping variables for coupled gas/solids simulations.
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
        bind(c, name="mfix_call_usr3_des")

      use run,            only: call_usr
      use particle_mod,   only: particle_t
      
      integer(c_int),   intent(in)    :: np
      type(particle_t), intent(inout) :: particles(np)
            
      if ( call_usr ) call usr3_des(np, particles)

   end subroutine call_usr3_des

   subroutine call_usr2_des( np, particles ) &
        bind(c, name="mfix_call_usr2_des")

      use run,            only: call_usr
      use particle_mod,   only: particle_t
      
      integer(c_int),   intent(in)    :: np
      type(particle_t), intent(inout) :: particles(np)
                  
      if ( call_usr ) call usr2_des(np, particles)

   end subroutine call_usr2_des



   
   subroutine des_time_loop_ops ( np, particles, subdt, dx, dy, dz, &
        & xlength, ylength, zlength, nstep )  &
        bind(C, name="mfix_des_time_loop_ops")
   
      use particle_mod
      use calc_collision_wall,     only: calc_dem_force_with_wall_stl
      use calc_force_dem_module,   only: calc_force_dem
      use cfnewvalues_module,      only: cfnewvalues
      use discretelement,          only: dtsolid, s_time, des_continuum_coupled
      use output_manager_module,   only: output_manager
      use run,                     only: call_usr

      integer(c_int),   intent(in   ) :: np
      real(c_real),     intent(in   ) :: subdt, xlength, ylength, zlength, dx, dy, dz
      type(particle_t), intent(inout) :: particles(np)
      integer(c_int),   intent(inout) :: nstep
      real(c_real)                    :: tow(np,3), fc(np,3)

      tow  = 0
      fc   = 0

      ! calculate forces from particle-wall collisions
      call calc_dem_force_with_wall_stl( particles, fc, tow, xlength, ylength, zlength, subdt )

      ! calculate forces from particle-particle collisions
      call calc_force_dem( particles, fc, tow, subdt )

      ! call user functions.
      if ( call_usr ) call usr1_des

      ! update position and velocities
      call cfnewvalues( particles, fc, tow, subdt )

      ! update time to reflect changes
      s_time = s_time + subdt

   end subroutine des_time_loop_ops

end module des_time_march_module
