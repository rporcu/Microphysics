module des_time_march_module

   use amrex_fort_module, only: c_real => amrex_real
   use iso_c_binding ,    only: c_int

   implicit none

   private
   ! public des_time_march

   ! Let's make some variables module variables to make the transition
   ! to C++ easier ( even though it is really a nasty thing to do ).
   real(c_real),   save  :: dtsolid_tmp

   ! Temporary variables when des_continuum_coupled is T to track
   ! changes in solid time step
   real(c_real),   save  :: TMP_DTS

   real(c_real),   allocatable :: fc(:,:), tow(:,:)

   ! Define interface for external functions so there will be no
   ! warning at compile time
   interface

      subroutine usr0_des()
      end subroutine usr0_des

      subroutine usr1_des()
      end subroutine usr1_des

      ! subroutine usr2_des_aos(particles)
      !    use particle_mod, only: particle_t
      !    type(particle_t), intent(inout) :: particles(:)
      ! end subroutine usr2_des_aos
      
      ! subroutine usr2_des(np, pos, vel, omega)
      !    import  c_real
      !    integer     , intent(in)    :: np
      !    real(c_real), intent(in)    :: pos(np,3), vel(np,3)
      !    real(c_real), intent(inout) :: omega(np,3)
      ! end subroutine usr2_des

      subroutine usr3_des(np, pos, vel, omega)
         import  c_real
         integer     , intent(in)    :: np
         real(c_real), intent(in)    :: pos(np,3), vel(np,3), omega(np,3)
      end subroutine usr3_des

   end interface


contains


   subroutine des_init_time_loop( np, state, phase, radius, vol, pos, vel, omega, drag,&
        & time, dt, dx, dy, dz, xlength, ylength, zlength, nstep, nsubsteps) &
        bind(C, name="mfix_des_init_time_loop")

      use discretelement,          only: dtsolid, s_time
      use discretelement,          only: des_continuum_coupled
      use output_manager_module,   only: output_manager
      use run,                     only: call_usr, tstop
      use param,                   only: zero
      use error_manager,           only: err_msg, ival, flush_err_msg

      integer(c_int), intent(in)    :: np
      real(c_real),   intent(in)    :: xlength, ylength, zlength, dx, dy, dz
      real(c_real),   intent(inout) :: time, dt

      real(c_real), intent(inout) :: vol(np), radius(np)
      real(c_real), intent(inout) :: drag(np,3), omega(np,3)
      real(c_real), intent(inout) :: pos(np,3), vel(np,3)

      integer(c_int), intent(inout) :: nstep, state(np), phase(np)
      integer(c_int), intent(out)   :: nsubsteps

      allocate ( tow(np,3), fc(np,3) )

      tow        = zero
      fc         = zero

      ! In case of restarts assign S_TIME from MFIX TIME
      S_TIME = TIME
      TMP_DTS = ZERO

      ! Initialize time stepping variables for coupled gas/solids simulations.
      if ( des_continuum_coupled ) then

         if ( DT >= DTSOLID ) then
            nsubsteps = ceiling(real(DT/DTSOLID))
         else
            nsubsteps = 1
            DTSOLID_TMP = DTSOLID
            DT = DTSOLID
         end if

         ! Initialize time stepping variable for pure granular simulations.
      else

         nsubsteps = ceiling(real((TSTOP-TIME)/DTSOLID))
         DT = DTSOLID
         call output_manager(np, time, dt, xlength, ylength, zlength, nstep, &
              state, radius, pos, vel, omega, 0)

      end if

      if (des_continuum_coupled) then
         write(ERR_MSG, 1000) trim(iVal(nsubsteps))
         call FLUSH_ERR_MSG(HEADER=.false., FOOTER=.false., LOG=.false.)
      else
         write(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(nsubsteps))
         call FLUSH_ERR_MSG(HEADER=.false., FOOTER=.false., LOG=.false.)
      endif

1000  FORMAT(/'DEM NITs: ',A)
1100  FORMAT(/'Time: ',g12.5,3x,'DT: ',g12.5,3x,'DEM NITs: ',A)

      IF(call_usr) call USR0_DES

   end subroutine des_init_time_loop


   subroutine des_init_time_loop_aos( np, particles, time, dt, dx, dy, dz, xlength, ylength, zlength, &
        & nstep, nsubsteps) &
        bind(C, name="mfix_des_init_time_loop_aos")

      use discretelement,          only: dtsolid, s_time
      use discretelement,          only: des_continuum_coupled
      use output_manager_module,   only: output_manager_aos
      use run,                     only: call_usr, tstop
      use param,                   only: zero
      use error_manager,           only: err_msg, ival, flush_err_msg
      use particle_mod,            only: particle_t
      
      integer(c_int),   intent(in)    :: np
      real(c_real),     intent(in)    :: xlength, ylength, zlength, dx, dy, dz
      real(c_real),     intent(inout) :: time, dt
      type(particle_t), intent(inout) :: particles(np)
      integer(c_int),   intent(inout) :: nstep
      integer(c_int),   intent(out)   :: nsubsteps

      allocate ( tow(np,3), fc(np,3) )
     
      tow        = zero
      fc         = zero

      ! In case of restarts assign S_TIME from MFIX TIME
      S_TIME = TIME
      TMP_DTS = ZERO

      ! Initialize time stepping variables for coupled gas/solids simulations.
      if ( des_continuum_coupled ) then

         if ( DT >= DTSOLID ) then
            nsubsteps = ceiling(real(DT/DTSOLID))
         else
            nsubsteps = 1
            DTSOLID_TMP = DTSOLID
            DT = DTSOLID
         end if

         ! Initialize time stepping variable for pure granular simulations.
      else

         nsubsteps = ceiling(real((TSTOP-TIME)/DTSOLID))
         DT = DTSOLID
         call output_manager_aos(np, time, dt, xlength, ylength, zlength, nstep, &
              particles, 0)

      end if

      if (des_continuum_coupled) then
         write(ERR_MSG, 1000) trim(iVal(nsubsteps))
         call FLUSH_ERR_MSG(HEADER=.false., FOOTER=.false., LOG=.false.)
      else
         write(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(nsubsteps))
         call FLUSH_ERR_MSG(HEADER=.false., FOOTER=.false., LOG=.false.)
      endif

1000  FORMAT(/'DEM NITs: ',A)
1100  FORMAT(/'Time: ',g12.5,3x,'DT: ',g12.5,3x,'DEM NITs: ',A)

      IF(call_usr) call USR0_DES

   end subroutine des_init_time_loop_aos


   

   subroutine des_finalize_time_loop( np, dt, pos, vel, omega ) &
        bind(C, name="mfix_des_finalize_time_loop")

      use discretelement, only: dtsolid
      use param,          only: zero
      use run,            only: call_usr

      integer(c_int), intent(in)    :: np
      real(c_real),   intent(inout) :: dt, pos(np,3), vel(np,3), omega(np,3)

      if ( call_usr ) call usr3_des(np, pos, vel, omega)

      ! When coupled, and if needed, reset the discrete time step accordingly
      if ( DT < DTSOLID_TMP ) then
         DTSOLID = DTSOLID_TMP
      endif

      if ( abs(TMP_DTS) > ZERO ) then
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      endif

      deallocate( tow, fc )

   end subroutine des_finalize_time_loop


   subroutine des_finalize_time_loop_aos( np, dt, particles ) &
        bind(C, name="mfix_des_finalize_time_loop_aos")

      use discretelement, only: dtsolid
      use param,          only: zero
      use run,            only: call_usr
      use particle_mod,   only: particle_t
      
      integer(c_int),   intent(in)    :: np
      real(c_real),     intent(inout) :: dt
      type(particle_t), intent(inout) :: particles(np)
      
      !if ( call_usr ) call usr3_des(np, pos, vel, omega)

      ! When coupled, and if needed, reset the discrete time step accordingly
      if ( DT < DTSOLID_TMP ) then
         DTSOLID = DTSOLID_TMP
      endif

      if ( abs(TMP_DTS) > ZERO ) then
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      endif

      deallocate( tow, fc )

   end subroutine des_finalize_time_loop_aos


   

   subroutine des_time_loop_ops( np, &
        & state, phase, radius, vol, mass, omoi, pos, vel, omega, acc, alpha, drag,&
        & time, dt, dx, dy, dz, xlength, ylength, zlength, nstep, quit )  &
        bind(C, name="mfix_des_time_loop_ops")

      use calc_collision_wall,     only: calc_dem_force_with_wall_stl
      use calc_force_dem_module,   only: calc_force_dem
      use cfnewvalues_module,      only: cfnewvalues
      use discretelement,          only: dtsolid, s_time, des_continuum_coupled
      use output_manager_module,   only: output_manager
      use run,                     only: call_usr

      integer(c_int), intent(in)    :: np
      real(c_real),   intent(in)    :: xlength, ylength, zlength, dx, dy, dz
      real(c_real),   intent(inout) :: time, dt

      real(c_real), intent(inout) :: vol(np), mass(np), radius(np), omoi(np)
      real(c_real), intent(inout) :: alpha(np,3), drag(np,3), omega(np,3)
      real(c_real), intent(inout) :: pos(np,3), vel(np,3), acc(np,3)

      integer(c_int), intent(inout) :: nstep, state(np), phase(np)
      integer(c_int), intent(out)   :: quit

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

      ! calculate forces from particle-wall collisions
      call calc_dem_force_with_wall_stl( phase, state, radius, &
           pos, vel, omega, fc, tow, xlength, ylength, zlength )

      ! calculate forces from particle-particle collisions
      call calc_force_dem(phase, radius, pos, vel, omega, state, &
           fc, tow)

      ! call user functions.
      if ( call_usr ) call usr1_des

      ! update position and velocities
      call cfnewvalues(np, state, mass, omoi, drag, &
           pos, vel, omega, fc, tow, &
           acc, alpha)

      ! update time to reflect changes
      s_time = s_time + dtsolid

      ! the following section targets data writes for dem only cases:
      if(.not.des_continuum_coupled) then
         ! keep track of time and number of steps for dem simulations
         time = s_time
         nstep = nstep + 1

         ! call the output manager to write res data.
         call output_manager(np, time, dt, &
              xlength, ylength, zlength, nstep, &
              state, radius, &
              pos, vel, omega, 0)
      endif  ! end if (.not.des_continuum_coupled)

      if(call_usr) call usr2_des(np, pos, vel, omega)

   end subroutine des_time_loop_ops


   subroutine des_time_loop_ops_aos( np, particles, time, dt, dx, dy, dz, &
        & xlength, ylength, zlength, nstep, quit )  &
        bind(C, name="mfix_des_time_loop_ops_aos")

      use particle_mod
      use calc_collision_wall,     only: calc_dem_force_with_wall_stl_aos
      use calc_force_dem_module,   only: calc_force_dem_aos
      use cfnewvalues_module,      only: cfnewvalues_aos
      use discretelement,          only: dtsolid, s_time, des_continuum_coupled
      use output_manager_module,   only: output_manager_aos
      use run,                     only: call_usr

      integer(c_int),   intent(in   ) :: np
      real(c_real),     intent(in   ) :: xlength, ylength, zlength, dx, dy, dz
      real(c_real),     intent(inout) :: time, dt
      type(particle_t), intent(inout) :: particles(np)
      integer(c_int),   intent(inout) :: nstep
      integer(c_int),   intent(out  ) :: quit

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

      ! calculate forces from particle-wall collisions
      call calc_dem_force_with_wall_stl_aos( particles, fc, tow, xlength, ylength, zlength )

      ! calculate forces from particle-particle collisions
      call calc_force_dem_aos( particles, fc, tow )

      ! call user functions.
      if ( call_usr ) call usr1_des

      ! update position and velocities
      call cfnewvalues_aos( particles, fc, tow )

      ! update time to reflect changes
      s_time = s_time + dtsolid

      ! the following section targets data writes for dem only cases:
      if ( .not. des_continuum_coupled ) then
         ! keep track of time and number of steps for dem simulations
         time  = s_time
         nstep = nstep + 1

         ! call the output manager to write res data.
         call output_manager_aos(np, time, dt,  xlength, ylength, zlength, &
              nstep, particles, 0 )
         
      end if 

      if (call_usr)  call usr2_des( size(particles), particles ) 

   end subroutine des_time_loop_ops_aos

end module des_time_march_module
