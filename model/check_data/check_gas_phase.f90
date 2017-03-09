module check_gas_phase_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  implicit none
  private

  public check_gas_phase

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Module name: CHECK_GAS_PHASE                                        !
  !  Purpose: Check the gas phase input section                          !
  !                                                                      !
  !  Author: P.Nicoletti                                Date: 02-DEC-91  !
  !          J.Musser                                   Date: 01-FEB-14  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_gas_phase

    use fld_const,     only: mu_g0, ro_g0, mw_avg
    use param1,        only: is_undefined, is_defined, zero
    use run,           only: IFILE_NAME
    use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, &
                           & ivar, ival, err_msg


    ! Initialize the error manager.
    call init_err_msg("CHECK_GAS_PHASE")


    ! CHECK MU_g0
    if (mu_g0 <= zero) then
       write(err_msg,1001) 'MU_G0', ival(mu_g0), trim(IFILE_NAME)
       call flush_err_msg(abort=.true.)
    endif


    ! CHECK MW_AVG
    ! When the species equations are not solved and the gas phase is
    ! compressible, verify that the user provided average molecular weight
    ! has a physical value. (This does not include the case where MW_AVG
    ! is UNDEFINED.)
    if (is_undefined(ro_g0)) then
       if (is_undefined(mw_avg)) then
          write(err_msg, 1001) 'MW_AVG', ival(mw_avg), trim(IFILE_NAME)
          call flush_err_msg(abort=.true.)
       endif
    else
       ! Gas density for incompressible flows must be positive.
       if (ro_g0 < zero) then
          write(err_msg, 1001) 'RO_g0', ival(ro_g0), trim(IFILE_NAME)
          call flush_err_msg(abort=.true.)
       endif
       ! Incompressible simulations do not need MW_AVG. Notify the user that
       ! the provided data is ignored.
       if (is_defined(mw_avg))then
          write(err_msg, 1100) 'RO_g0 is specified'
          call flush_err_msg
       endif

    endif


    ! Finalize the error manager
    call finl_err_msg

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the ',A,' file.')


1100 format('Message 2000: MW_AVG is not needed when ',A,'.')

  end subroutine check_gas_phase
end module check_gas_phase_module
