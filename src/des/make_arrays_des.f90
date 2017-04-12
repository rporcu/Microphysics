!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine make_arrays_des(max_pip, particle_state, &
  des_radius, ro_sol, pvol, pmass, omoi) &
  bind(C, name="mfix_make_arrays_des")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

! Global data .......................................................//
  use discretelement, only: normal_particle
  use constant, only: pi
  use param1, only: zero

  implicit none

  integer(c_int), intent(in   ) :: max_pip
  integer(c_int), intent(in   ) :: particle_state(max_pip)

  real(c_real),   intent(in   ) :: des_radius(max_pip)
  real(c_real),   intent(in   ) :: ro_sol(max_pip)
  real(c_real),   intent(  out) :: pvol(max_pip)
  real(c_real),   intent(  out) :: pmass(max_pip)
  real(c_real),   intent(  out) :: omoi(max_pip)


! Local variables ...................................................//
  integer :: lc

  ! Populate the particle property arrays.
  do lc = 1, max_pip
     if(particle_state(lc) == normal_particle) then
        pvol(lc) = (4.0d0/3.0d0)*pi*des_radius(lc)**3
        pmass(lc) = pvol(lc)*ro_sol(lc)
        omoi(lc) = 2.5d0/(pmass(lc)*des_radius(lc)**2)
     endif
  enddo

end subroutine make_arrays_des
