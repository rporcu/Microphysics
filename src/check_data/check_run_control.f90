module check_run_control_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use error_manager,  only: init_err_msg, flush_err_msg, finl_err_msg, &
                            err_msg, ivar, ival

    use param, only: undefined_c, zero
    use param, only: is_defined, is_undefined

  implicit none
  private

  public check_run_control

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_RUN_CONTROL                                       !
!  Purpose: Check the run control namelist section                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_run_control(dt)

    use run,    only: nlog, discretize
    use param, only: dim_eqs


    real(rt), intent(in) :: dt

    integer  :: lc

    CALL init_err_msg("CHECK_RUN_CONTROL")

    if(is_defined(dt)) then

       if (dt < zero) then
          write(err_msg,1002) 'DT', dt
          call flush_err_msg(abort=.true.)
       endif

    endif

    do lc = 1,dim_eqs
       if(discretize(lc) /= 0 .and. discretize(lc) /= 2) then
          write(err_msg,2002) trim(ivar('DISCRETIZE',lc)),&
            trim(ival(discretize(lc)))
          call flush_err_msg(abort=.true.)
       endif
    enddo

    if(nlog <= 0) then
       write(err_msg,1001) 'NLOG', ival(nlog)
       call flush_err_msg(abort=.true.)
    endif

2002 format('Error 2002: Invalid option ', A,' = ', A, '.',/  &
       'Please correct the input deck.')

    ! Clear the error manager
    call finl_err_msg

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
          'Please correct the input deck.')

1002 format('Error 1002: Illegal or unknown input: ',A,' = ',G14.4,/  &
          'Please correct the input deck.')

  end subroutine check_run_control

end module check_run_control_module
