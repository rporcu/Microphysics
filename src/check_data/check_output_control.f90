module check_output_control_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use run,            only: IFILE_NAME
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
  !  Author: P. Nicoletti                               Date: 27-NOV-91  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_output_control

    use output, only: RES_DT, NLOG
    use param1, only: ZERO, IS_UNDEFINED

    ! Initialize the error manager.
    call init_err_msg("CHECK_OUTPUT_CONTROL")

    ! Check the values specified for the RES file.
    if (IS_UNDEFINED(RES_DT))then
       write(ERR_MSG,1000) 'RES_DT', trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    elseif(RES_DT <= ZERO) then
       write(ERR_MSG,1002) 'RES_DT', RES_DT, trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
    
    ! Verify that the LOG frequency is valid.
    if(NLOG <= 0) then
       write(ERR_MSG,1003) 'NLOG', NLOG, trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
    
    ! Finalize the error manager.
    call finl_err_msg

    
1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
          'correct the ',A,' file.')
    
1002 format('Error 1002: Illegal or unknown input: ',A,' = ',E14.6,/  &
          'Please correct the ',A,'  file.')
    
1003 format('Error 1003: Illegal or unknown input: ',A,' = ',I4,/     &
          'Please correct the ',A,'  file.')
    
  end subroutine check_output_control

end module check_output_control_module
