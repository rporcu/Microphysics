subroutine set_particle_properties(pstate, pradius, pdensity, pvol, pmass, omoi, omega) &
     bind(C, name="mfix_set_particle_properties")

  use amrex_fort_module, only: c_real => amrex_real
  use iso_c_binding ,    only: c_int
  use discretelement,    only: normal_particle
  use constant,          only: pi
  use param,            only: zero

  implicit none

  integer(c_int), intent(in)  :: pstate
  real(c_real),   intent(in)  :: pradius, pdensity
  real(c_real),   intent(out) :: pvol, pmass, omoi, omega

  if ( pstate == normal_particle ) then
     pvol  = (4.0d0/3.0d0)*pi*pradius**3
     pmass = pvol * pdensity
     omoi  = 2.5d0/(pmass * pradius**2)
     omega = 0.d0
  endif

end subroutine set_particle_properties
