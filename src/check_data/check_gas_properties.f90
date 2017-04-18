module check_gas_prop_module

  use error_manager, only: init_err_msg, finl_err_msg, flush_err_msg, &
                           err_msg, ivar, ival

  use param, only: is_undefined, is_defined, zero

  implicit none
  private

  public check_gas_properties

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_GAS_PHASE                                        !
!  Purpose: Check the gas phase input section                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_gas_properties

    use fld_const, only: mu_g0, ro_g0, mw_avg

    call init_err_msg("CHECK_GAS_PHASE")

    if(ro_g0 > zero) then
       if (is_defined(ro_g0)) then
          if (ro_g0 < zero) then ! incompressible
             write(err_msg, 1001) 'RO_g0', ival(ro_g0)
             call flush_err_msg(abort=.true.)
          endif
          if (is_defined(mw_avg))then
             write(err_msg, 1100) 'RO_g0 is specified'
             call flush_err_msg
          endif
1100 format('Message 2000: MW_AVG is not needed when ',A,'.')

       else ! compressible

          ! Verify that the user provided average molecular weight
          if (is_undefined(mw_avg)) then
             write(err_msg, 1001) 'MW_AVG', ival(mw_avg)
             call flush_err_msg(abort=.true.)
          endif
       endif

       if (mu_g0 < zero) then
          write(err_msg,1001) 'MU_G0', ival(mu_g0)
          call flush_err_msg(abort=.true.)
       endif
    endif

    call finl_err_msg

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the input deck.')

  end subroutine check_gas_properties
end module check_gas_prop_module
