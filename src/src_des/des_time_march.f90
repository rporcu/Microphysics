module des_time_march_module

   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding ,    only: c_int

   implicit none

   private

   ! Define interface for external functions so there will be no
   ! warning at compile time
   interface

      subroutine usr0_des ()
      end subroutine usr0_des

      subroutine usr1_des ()
      end subroutine usr1_des

   end interface

contains

   subroutine des_init_time_loop (tstart, dt, nsubsteps, subdt) &
        bind(C, name="des_init_time_loop")

      use discretelement,  only: dtsolid
      use run,             only: des_tstart, des_dt

      real(rt),   intent(in   ) :: tstart, dt
      integer(c_int), intent(  out) :: nsubsteps
      real(rt),   intent(  out) :: subdt

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

   subroutine call_usr2_des (np, particles) &
        bind(c, name="call_usr2_des")

      use particle_mod,   only: particle_t

      integer(c_int),   intent(in)    :: np
      type(particle_t), intent(inout) :: particles(np)

      call usr2_des(np, particles)

   end subroutine call_usr2_des

   subroutine call_usr3_des (np, particles) &
        bind(c, name="call_usr3_des")

      use particle_mod,   only: particle_t

      integer(c_int),   intent(in)    :: np
      type(particle_t), intent(inout) :: particles(np)

      call usr3_des(np, particles)

   end subroutine call_usr3_des

end module des_time_march_module
