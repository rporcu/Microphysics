!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: check_domain                                            !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine check_domain(dx,dy,dz,xlength,ylength,zlength,domlo,domhi) &
  bind(C, name="check_domain")

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use check_initial_conditions_module, only: check_initial_conditions
  use check_boundary_conditions_module, only: check_boundary_conditions
  use check_point_sources_module, only: check_point_sources

  integer(c_int), intent(in) :: domlo(3),domhi(3)
  real(c_real)  , intent(in) :: dx, dy ,dz
  real(c_real)  , intent(in) :: xlength, ylength, zlength

  call check_initial_conditions(dx,dy,dz,domlo,domhi)
  call check_boundary_conditions(dx,dy,dz,xlength,ylength,zlength,domlo,domhi)
  call check_point_sources(dx,dy,dz)

end subroutine check_domain
