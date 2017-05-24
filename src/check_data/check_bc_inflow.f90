module check_bc_inflow_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use param,         only: undefined, one, zero, is_undefined, is_defined
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg,    &
                         &  ivar,  ival, err_msg

  implicit none
  private

  public check_bc_mass_inflow
  public check_bc_p_inflow

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: check_bc_mass_inflow                                     !
! Purpose: Provided a detailed error message on BC                     !
!                                                                      !
! Comments:                                                            !
!     The velocities at the inflow face are fixed and the momentum     !
!     equations are not solved in the inflow cells. Since the flow is  !
!     into the domain all other scalars that are used need to be       !
!     specified (e.g., mass fractions, void fraction, etc.,)           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_bc_mass_inflow(m_tot, skip, bcv)

    use bc,        only: bc_ep_g, bc_p_g, bc_ep_s, bc_massflow_g
    use param    , only: dim_m
    use fld_const, only: ro_g0
    use param,    only: equal

    integer, intent(in) :: BCV, M_TOT
    logical, intent(in) :: SKIP(DIM_M)
    integer             :: M
    real(c_real)        :: SUM_EP


    CALL init_err_msg("CHECK_BC_MASS_INFLOW")


    ! Verify compressible boundary condition variables.
    if(is_undefined(ro_g0)) then
       if(is_undefined(bc_p_g(bcv))) then
          if(is_undefined(bc_massflow_g(bcv)) .and.                   &
               abs(bc_massflow_g(bcv)) > zero) then
             write(err_msg, 1100) trim(ivar('BC_P_g',bcv))
             call flush_err_msg(abort=.true.)
          endif

1100 format('Error 1100: ',A,' must be specified for compressible ',  &
        'flows',/'when specifying BC_MASSFLOW_g to make the ',        &
        'conversion to velocity.',/'Please correct the input deck.')

       elseif(BC_P_G(BCV) <= ZERO) then
          write(ERR_MSG, 1101) BCV, trim(iVal(BC_P_G(BCV)))
          call flush_err_msg(ABORT=.true.)
       endif
1101   format('Error 1101: Pressure must be greater than zero for ',    &
          'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/&
          'Please correct the input deck.')
    endif

    ! Initialize the sum of the total volume fraction.
    sum_ep = bc_ep_g(bcv)
    do m = 1, m_tot
       sum_ep = sum_ep + bc_ep_s(bcv,m)
    enddo

    ! verify that the volume fractions sum to one.
    if(.not.equal(sum_ep,one)) then
       write(err_msg,1215) bcv, trim(ival(sum_ep))
       call flush_err_msg(abort=.true.)
    endif
1215 format('Error 1215: Illegal boundary condition region: ',I3,'. ',&
          'Sum of volume',/'fractions does NOT equal ONE. (SUM = ',A,   &
          ')',/'Please correct the input deck.')

    call finl_err_msg

  end subroutine check_bc_mass_inflow


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: check_bc_p_inflow                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided detailed error message on bc                       !
!                                                                      !
! Comments:                                                            !
!     Unlike the MI boundary, for the PI boundary the velocities at    !
!     the inflow face are calculated by solving the momentum eqns      !
!     and are not fixed. In this way, the PI is similar to the PO      !
!     except that the flow is into the domain and hence all other      !
!     scalars (e.g., mass fractions, void fraction, temperature,       !
!     etc.,) at the inflow cells need to be specified. To satisfy      !
!     the error routines at the start of the simulation, both the      !
!     tangential and normal components at the inflow also need to      !
!     be specified. The velocities values essentially serve as IC.     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_p_inflow(m_tot, skip, bcv)

    use param    , only: dim_m
    use fld_const, only: ro_g0
    use bc,        only: bc_p_g

    integer, intent(in) :: bcv, m_tot
    logical, intent(in) :: skip(dim_m)

    call init_err_msg("CHECK_BC_P_INFLOW")

    ! Remove checks on bc_ep_g; using routine check_bc_outflow
    if (is_undefined(bc_p_g(bcv))) then
       write(err_msg,1000) 'BC_P_g', bcv
       call flush_err_msg(abort=.true.)
1000  format('Error 1000: Required input not specified: ',A,/&
         'Please correct the input deck.')

    elseif (bc_p_g(bcv)<=zero .and. is_undefined(ro_g0)) then
       write(err_msg, 1101) bcv, trim(ival(bc_p_g(bcv)))
       call flush_err_msg(abort=.true.)
1101   format('Error 1101: Pressure must be greater than zero for ',    &
          'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/&
          'Please correct the input deck.')
    endif

    call finl_err_msg

  end subroutine check_bc_p_inflow

end module check_bc_inflow_module
