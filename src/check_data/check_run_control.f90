module check_run_control_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg, &
                          & ivar, ival, err_msg

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
  subroutine check_run_control(time, dt)

    use run,    only: RUN_TYPE, DESCRIPTION, TSTOP
    use param1, only: UNDEFINED_C, IS_UNDEFINED, ZERO
! Flag to adjust dt when residuals diverge
      use run, only: DETECT_STALL

    real(c_real), intent(inout) :: time, dt


    ! Initialize the error manager.
    CALL init_err_msg("CHECK_RUN_CONTROL")

    ! Clear out the run description if not specified.
    if (DESCRIPTION == UNDEFINED_C) DESCRIPTION = ' '

    ! Verify that DT is valid.
    if (DT < ZERO) then
       write(ERR_MSG,1002) 'DT', DT
       call flush_err_msg(ABORT=.true.)

       ! Steady-state simulation.
    elseif(IS_UNDEFINED(DT) .or. IS_UNDEFINED(DT)) then
       DETECT_STALL = .FALSE.
       TIME = ZERO

       ! Transient simulation.
    else
       ! Verify the remaining time settings.
       if (IS_UNDEFINED(TIME)) then
          write(ERR_MSG,1000) 'TIME'
          call flush_err_msg(ABORT=.true.)

       elseif (IS_UNDEFINED(TSTOP)) then
          write(ERR_MSG,1000) 'TSTOP'
          call flush_err_msg(ABORT=.true.)

       elseif (TIME < ZERO) then
          write(ERR_MSG,1002)'TIME', TIME
          call flush_err_msg(ABORT=.true.)

       elseif (TSTOP < ZERO) then
          write(ERR_MSG,1002) 'TSTOP', TSTOP
          call flush_err_msg(ABORT=.true.)
       endif
    endif

    ! Verify the run type.
    if(.not.(RUN_TYPE=='NEW')) then
       write(ERR_MSG,1001) 'RUN_TYPE', RUN_TYPE
       call flush_err_msg(ABORT=.true.)
    endif


    ! Clear the error manager
    call finl_err_msg


1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
          'correct the input deck.')

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
          'Please correct the input deck.')

1002 format('Error 1002: Illegal or unknown input: ',A,' = ',G14.4,/  &
          'Please correct the input deck.')

  end subroutine check_run_control

end module check_run_control_module
