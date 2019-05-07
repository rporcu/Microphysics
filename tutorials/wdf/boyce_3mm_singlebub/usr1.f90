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

  real(rt) :: inj_dt

   inj_dt = 154.2d-3
!  inj_dt = 101.4d-3
!  inj_dt =  66.7d-3
!  inj_dt =  51.4d-3
!  inj_dt =  25.0d-3

  if ((time .GE. 1.0d0 - inj_dt) .AND. (time .LE. 1.0d0)) then
    bc_u_g(10) = 50.0d0
  else
    bc_u_g(10) = 0.0d0
  endif 

!
! umf bc: uncomment to measure U_mf 
!  
!  bc_u_g(1) = max(0.0d0, min(1.4d0*time, 1.5d0 - time*0.1d0))

  return

end subroutine usr1
