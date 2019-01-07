!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: output                                                      !
!                                                                      !
!  Purpose: Contain data for output control.                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module output

  use param, only: dim_usr
  use amrex_fort_module, only : rt => amrex_real

  ! Interval at which user-defined output files are updated.
  real(rt) :: usr_dt (dim_usr)
  real(rt) :: usr_time(dim_usr) = 0.0d0

  real(rt) :: USR_X_w(dim_usr), USR_X_e(dim_usr)
  real(rt) :: USR_Y_s(dim_usr), USR_Y_n(dim_usr)
  real(rt) :: USR_Z_b(dim_usr), USR_Z_t(dim_usr)

contains

  !......................................................................!
  !                                                                      !
  !......................................................................!
  integer function newUnit(punit)

    integer, intent(in   ), optional :: punit

    logical :: is_open
    integer :: lc1

    newUnit = -1

    if(present(punit)) then
       inquire(unit=punit, opened=is_open)
       if(.not.is_open) newUnit = punit
    endif

    if(newUnit <= 6 .or. newUnit == -1) then
       do lc1=100, 9999
          inquire(unit=lc1, opened=is_open)
          if(.not.is_open) then
             newUnit = lc1
             return
          endif
       enddo
    endif

  end function newUnit

end module output
