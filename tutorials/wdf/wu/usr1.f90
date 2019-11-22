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
  use bc, only: bc_u_g

  implicit none

  real(rt),   intent(in ) :: time
  real(rt) :: usr_pi, usr_umf

  usr_pi = 4.0d0*DATAN(1.0d0)
  
  usr_umf = MIN(0.02d0*time, 0.041d0)
  bc_u_g(1) = usr_umf*(2.38d0 + 2.03d0*DSIN(10.0d0*usr_pi*time))

  return

end subroutine usr1
