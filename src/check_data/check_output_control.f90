module check_output_control_module

  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg, &
                          & ivar, ival, err_msg

  implicit none
  private

  public check_output_control

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_OUTPUT_CONTROL                                    !
!  Purpose: Check the output control namelist section                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_output_control

    use output, only: NLOG

    ! Initialize the error manager.
    call init_err_msg("CHECK_OUTPUT_CONTROL")

    ! Verify that the LOG frequency is valid.
    if(NLOG <= 0) then
       write(ERR_MSG,1003) 'NLOG', NLOG
       call flush_err_msg(ABORT=.true.)
    endif

    ! Finalize the error manager.
    call finl_err_msg

1003 format('Error 1003: Illegal input: ',A,' = ',I4,/     &
          'Please correct the input deck.')

  end subroutine check_output_control

end module check_output_control_module
