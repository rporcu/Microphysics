module des_time_march_module

   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding ,    only: c_int

   implicit none

   private

contains

   subroutine des_init_time_loop (dt, nsubsteps, subdt) &
        bind(C, name="des_init_time_loop")

      use discretelement,  only: dtsolid

      real(rt),       intent(in   ) :: dt
      integer(c_int), intent(  out) :: nsubsteps
      real(rt),       intent(  out) :: subdt

      ! Initialize time stepping variables
      if ( dt >= dtsolid ) then
          nsubsteps = ceiling ( real ( dt / dtsolid ) )
          subdt     =  dt / nsubsteps
      else
          nsubsteps = 1
          subdt     = dt
      end if

   end subroutine des_init_time_loop

end module des_time_march_module
