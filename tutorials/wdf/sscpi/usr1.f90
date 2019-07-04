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
  real(rt) :: usr_U

!-----------------------------------------------------!
! BCs          case :  2     3     4       (*Umf)     !
!      xlo.velocity =  2.19  3.28  4.38    (m/s)      !
!-----------------------------------------------------!

! pull final U from BC ID 101 
  usr_U = bc_u_g(101)

! ramp up over 2s 
  bc_u_g(1) = min(1.0d0, max(usr_U, 0.5d0*usr_U*time))

  return

end subroutine usr1
