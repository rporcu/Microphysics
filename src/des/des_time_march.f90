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

   
   subroutine des_set_dt( time, dt, nsteps, dtsolid_tmp, tmp_dts) &
        bind(C, name="mfix_des_set_dt")

      use discretelement,  only: dtsolid, des_continuum_coupled
      use run,             only: tstop
      use param,           only: zero
      
      real(c_real),     intent(inout) :: time, dt
      integer(c_int),   intent(inout) :: nsteps
      real(c_real),     intent(out)   :: dtsolid_tmp, tmp_dts
      
      ! In case of restarts assign S_TIME from MFIX TIME
      TMP_DTS = ZERO


      ! ! In case of restarts assign S_TIME from MFIX TIME
      ! s_time = time
      ! tmp_dts = zero

      ! ! Initialize time stepping variables for coupled gas/solids simulations.
      ! if ( dt >= dtsolid ) then
      !    nsteps = ceiling(real(dt/dtsolid))
      ! else
      !    nsteps = 1
      !    dtsolid_tmp = dtsolid
      !    dtsolid = dt
      ! end if
      
      ! Initialize time stepping variables for coupled gas/solids simulations.
      if ( des_continuum_coupled ) then

         if ( DT >= DTSOLID ) then
            nsteps = ceiling(real(DT/DTSOLID))
         else
            nsteps = 1
            DTSOLID_TMP = DTSOLID
            DT = DTSOLID
         end if

         ! Initialize time stepping variable for pure granular simulations.
      else

         nsteps = ceiling(real((TSTOP-TIME)/DTSOLID))
         DT = DTSOLID
         
      end if
      
   end subroutine des_set_dt

 

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

   
   subroutine des_finalize_time_loop( dt, dtsolid_tmp, tmp_dts ) &
        bind(C, name="mfix_des_finalize_time_loop")

      use discretelement, only: dtsolid
      use param,          only: zero
      
      real(c_real),     intent(inout) :: dt
      real(c_real),     intent(inout) :: dtsolid_tmp, tmp_dts
      
      ! When coupled, and if needed, reset the discrete time step accordingly
      if ( DT < DTSOLID_TMP ) then
         DTSOLID = DTSOLID_TMP
      endif

      if ( abs(TMP_DTS) > ZERO ) then
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      endif

   end subroutine des_finalize_time_loop



   subroutine des_check_dt ( quit, dt, time, tmp_dts )  &
        bind(C, name="mfix_des_check_dt") 

      use discretelement,  only: dtsolid, s_time, des_continuum_coupled

      integer(c_int),   intent(out  ) :: quit
      real(c_real),     intent(in   ) :: dt, time
      real(c_real),     intent(inout) :: tmp_dts
      
      quit = 0
      
      if ( des_continuum_coupled ) then
         ! If the current time in the discrete loop exceeds the current time in
         ! the continuum simulation, exit the discrete loop
         if ( S_TIME > (TIME+DT) ) then
            quit = 1
            return
         end if
         ! If next time step in the discrete loop will exceed the current time
         ! in the continuum simulation, modify the discrete time step so final
         ! time will match
         if ( (S_TIME+DTSOLID) > (TIME+DT) ) then
            TMP_DTS = DTSOLID
            DTSOLID = TIME + DT - S_TIME
         end if
      end if

   end subroutine des_check_dt

 
   subroutine des_time_loop_ops ( np, particles, time, dt, dx, dy, dz, &
        & xlength, ylength, zlength, nstep, dtsolid_tmp, tmp_dts )  &
        bind(C, name="mfix_des_time_loop_ops")

      use particle_mod
      use calc_collision_wall,     only: calc_dem_force_with_wall_stl
      use calc_force_dem_module,   only: calc_force_dem
      use cfnewvalues_module,      only: cfnewvalues
      use discretelement,          only: dtsolid, s_time, des_continuum_coupled
      use output_manager_module,   only: output_manager
      use run,                     only: call_usr

      integer(c_int),   intent(in   ) :: np
      real(c_real),     intent(in   ) :: xlength, ylength, zlength, dx, dy, dz
      real(c_real),     intent(inout) :: time, dt
      type(particle_t), intent(inout) :: particles(np)
      integer(c_int),   intent(inout) :: nstep
      real(c_real)                    :: tow(np,3), fc(np,3)
      real(c_real),     intent(inout) :: dtsolid_tmp, tmp_dts

      tow  = 0
      fc   = 0

      ! calculate forces from particle-wall collisions
      call calc_dem_force_with_wall_stl( particles, fc, tow, xlength, ylength, zlength )

      ! calculate forces from particle-particle collisions
      call calc_force_dem( particles, fc, tow )

      ! call user functions.
      if ( call_usr ) call usr1_des

      ! update position and velocities
      call cfnewvalues( particles, fc, tow, dtsolid )

      ! update time to reflect changes
      s_time = s_time + dtsolid

      ! the following section targets data writes for dem only cases:
      if ( .not. des_continuum_coupled ) then
         ! keep track of time and number of steps for dem simulations
         time  = time + dtsolid 
         nstep = nstep + 1
         
      end if 

   end subroutine des_time_loop_ops

end module des_time_march_module
