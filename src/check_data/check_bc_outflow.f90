module check_bc_outflow_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use param1,         only: one, undefined, zero, is_undefined, is_defined, equal
  use error_manager,  only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival


  implicit none
  private

  public check_bc_outflow
  public check_bc_p_outflow

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_OUTFLOW                                         !
!                                                                      !
! Purpose: Provided a detailed error message concerning specification  !
! of bc_ep_g + bc_ep_s at a outflow boundary (and pressure inflow)     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_outflow(M_TOT, BCV)

     use bc,     only: bc_ep_g, bc_ep_s
     use param1, only: one, equal

    integer, intent(in) :: BCV, M_TOT
    integer             :: M
    real(c_real)        :: SUM_EP
    logical             :: FLAG_WARNING

    FLAG_WARNING = .true.

    CALL init_err_msg("CHECK_BC_OUTFLOW")

    if (.not.equal(bc_ep_g(bcv),one)) then
       sum_ep = bc_ep_g(bcv)
       do m = 1, m_tot
          if(flag_warning) then
             write(err_msg, 1101) trim(iVar('BC_EP_g',BCV))
             call flush_err_msg
             flag_warning = .false.
          endif
          sum_ep = sum_ep + bc_ep_s(bcv,m)
       enddo
       if(.not.equal(sum_ep,one)) then
          write(err_msg,1103) bcv, trim(ival(sum_ep))
          call flush_err_msg(abort=.true.)
       endif

   ! bc_ep_g is not defined but check if any bc_ep_s are defined
    else

       sum_ep = zero
       do m = 1, m_tot
          if(bc_ep_s(bcv,m) > zero) then
             write(err_msg, 1101) trim(ivar('BC_EP_s',bcv,m))
             call flush_err_msg
             sum_ep = sum_ep + bc_ep_s(bcv,m)
          endif
       enddo
       if(sum_ep > one) then
          write(err_msg,1103) bcv, trim(ival(sum_ep))
          call flush_err_msg(abort=.true.)
       endif

    endif

    call finl_err_msg

1101 format('Warning 1101: ',A,' should not be specified for ', &
          'outflow BCs',/'with DEM/PIC runs except for a mass outflow ',&
          'boundary with specified ',/ 'flow rate(s). In this case ',&
          'volume fraction data is used for ',/ 'conversion to ',&
          'velocity(s). However, the solids volume fraction data ',/&
          'is effectively disregarded and it is the solids velocity ',&
          'that is ',/'used to direct any solids at the boundary.')

1103 format('Error 1103: Illegal boundary condition region: ',I3,'. ',&
          'Sum of volume',/'fractions does NOT equal ONE. (SUM = ',A,   &
          ')',/'Please correct the input deck.')

  end subroutine check_bc_outflow


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: CHECK_BC_P_OUTFLOW                                       !
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Provided a detailed error message on bc                     !
  !                                                                      !
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_p_outflow(BCV)

    use fld_const, only: RO_g0
    use bc       , only: BC_P_g

    integer, intent(in) :: BCV

    call init_err_msg("CHECK_BC_P_OUTFLOW")

    if (is_undefined(bc_p_g(bcv))) then
       write(err_msg,1000) trim(ivar('BC_P_g',bcv))
       call flush_err_msg(abort=.true.)

    elseif (bc_p_g(bcv)<=zero .and. is_undefined(ro_g0)) then
       write(err_msg, 1100) bcv, trim(ival(bc_p_g(bcv)))
       call flush_err_msg(abort=.true.)
    endif

1000 format('Error 1000: Required input not specified: ',A,/&
        'Please correct the input deck.')

1100 format('Error 1100: Pressure must be greater than zero for ',    &
        'compressible flow',/3x,'BC_P_g(',I3,') = ',A,/&
        'Please correct the input deck.')

    ! Clean up and return.
    call finl_err_msg

  end subroutine check_bc_p_outflow

end module check_bc_outflow_module
