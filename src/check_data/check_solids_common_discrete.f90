module check_solids_common_discrete_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg,    &
       &  ivar,  ival, err_msg

  implicit none
  private

  public check_solids_common_discrete

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_SOLIDS_COMMON_DISCRETE                            !
  !  Author: J.Musser                                   Date: 02-FEB-14  !
  !                                                                      !
  !  Purpose:                                                            !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_solids_common_discrete

    use constant,       only: mmax,  d_p0
    use param1,         only: undefined
    use discretelement, only: des_intg_method, intg_adams_bashforth, &
                            & intg_euler

    integer :: lM  ! Solids phase Index

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_COMMON_DISCRETE")

    ! Check for valid integration method
    select case(trim(DES_INTG_METHOD))
    case ('EULER')
       intg_euler = .true.
       intg_adams_bashforth = .false.
    case ('ADAMS_BASHFORTH')
       intg_euler = .false.
       intg_adams_bashforth = .true.
    case DEFAULT
       write(err_msg,2020) trim(des_intg_method)
       call flush_err_msg(abort=.true.)

2020   format('Error 2020:Invalid DES_INGT_METHOD: ',A,/'Please ',      &
            'correct the input deck.')

    end select

    call finl_err_msg

  end subroutine check_solids_common_discrete

end module check_solids_common_discrete_module
