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
  real(c_real) :: usr_dt (dim_usr)
  real(c_real) :: usr_time(dim_usr) = 0.0d0

end module output
