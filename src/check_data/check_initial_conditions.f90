module check_initial_conditions_module

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public check_initial_conditions

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS                                !
!                                                                      !
!  Purpose: check the initial conditions input section                 !
!     - check geometry of any specified IC region                      !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_initial_conditions(dx,dy,dz,domlo,domhi) &
      bind(C, name="check_initial_conditions")

    integer(c_int), intent(in) :: domlo(3),domhi(3)
    real(rt)      , intent(in) :: dx, dy, dz


 end subroutine check_initial_conditions

end module check_initial_conditions_module
