!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Purpose: This routine is called from the time loop and is           !
!           user-definable.  The user may insert code in this routine  !
!           or call appropriate user defined subroutines.  This        !
!           can be used for setting or checking errors in quantities   !
!           that vary with time.  This routine is not called from an   !
!           IJK loop, hence all indices are undefined.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine usr1(time)

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int, c_char

  implicit none

  real(rt),   intent(in ) :: time

  return

end subroutine usr1
