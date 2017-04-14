!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: output                                                      !
!                                                                      !
!  Purpose: Contain data for output control.                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module output

  use param, only: dim_usr
  use amrex_fort_module, only : c_real => amrex_real

  ! Interval at which user-defined output files are updated.
  real(c_real) :: usr_dt (dim_usr), usr_time(dim_usr)
  ! Interval in number of time steps at which LOG file is written

end module output
