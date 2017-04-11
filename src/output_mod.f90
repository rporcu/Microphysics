!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: output                                                      !
!                                                                      !
!  Purpose: Contain data for output control.                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module output

  use param, only: dimension_usr
  use amrex_fort_module, only : c_real => amrex_real

  ! Interval at which user-defined output files are updated.
  real(c_real) :: usr_dt (dimension_usr), usr_time(dimension_usr)
  ! Interval in number of time steps at which LOG file is written
  integer :: nlog
  ! Flag to display messages and residuals on the screen
  logical :: full_log

end module output
