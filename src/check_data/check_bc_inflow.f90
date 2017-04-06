module check_bc_inflow_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use param1,         only: undefined, one, zero, is_undefined, is_defined
  use run,            only: IFILE_NAME
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
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Provided a detailed error message on BC                     !
  !                                                                      !
  ! Comments:                                                            !
  !     The velocities at the inflow face are fixed and the momentum     !
  !     equations are not solved in the inflow cells. Since the flow is  !
  !     into the domain all other scalars that are used need to be       !
  !     specified (e.g., mass fractions, void fraction, etc.,)           !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_bc_mass_inflow(M_TOT, SKIP, BCV)

    use bc,        only: bc_ep_g, bc_p_g, bc_rop_s, bc_ep_s, bc_massflow_g
    use param    , only: dim_m
    use fld_const, only: ro_g0
    use constant , only: ro_s0
    use toleranc , only: compare

    integer, intent(in) :: BCV, M_TOT
    logical, intent(in) :: SKIP(DIM_M)
    integer             :: M
    real(c_real)        :: SUM_EP
    real(c_real)        :: BC_ROs(DIM_M) ! Solids phase density in BC region.


    CALL init_err_msg("CHECK_BC_MASS_INFLOW")

    ! Check gas phase volume fraction.
    if(IS_UNDEFINED(BC_EP_G(BCV))) then
       write(ERR_MSG, 1000) trim(iVar('BC_EP_g',BCV)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif

    ! Verify compressible boundary condition variables.
    if(IS_UNDEFINED(RO_G0)) then
       if(IS_UNDEFINED(BC_P_G(BCV))) then
          if(IS_UNDEFINED(BC_MASSFLOW_G(BCV)) .and.                   &
               abs(BC_MASSFLOW_G(BCV)) > ZERO) then
             write(ERR_MSG, 1100) trim(iVar('BC_P_g',BCV)), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
1100      format('Error 1100: ',A,' must be specified for compressible ',  &
               'flows',/'when specifying BC_MASSFLOW_g to make the ',        &
               'conversion to velocity.',/'Please correct the ',A,' file.')
          
       elseif(BC_P_G(BCV) <= ZERO) then
          write(ERR_MSG, 1101) BCV, trim(iVal(BC_P_G(BCV))), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
1101   format('Error 1101: Pressure must be greater than zero for ',    &
            'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/'Please ',  &
            'correct the ',A,' file.')
    endif
    
    ! Calculate the solids volume fraction from the gas phase if there is
    ! only one solids phase.
    if(M_TOT == 1 .and. IS_UNDEFINED(BC_EP_S(BCV,1))) then
       BC_EP_S(BCV,1) = ONE - BC_EP_g(BCV)
    endif
    
    ! Bulk density or solids volume fraction must be explicitly defined
    ! if there are more than one solids phase.
    if(M_TOT > 1 .and. .not.COMPARE(BC_EP_g(BCV),ONE)) then
       do M = 1, M_TOT
          if(IS_UNDEFINED(BC_ROP_S(BCV,M)) .and. &
               IS_UNDEFINED(BC_EP_S(BCV,M))) then
             write(ERR_MSG, 1200) M, BCV, 'BC_ROP_s and BC_EP_s', trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
       enddo
    endif
1200 format('Error 1200: Insufficient solids phase ',I2,' data ',     &
          'for BC',I3,'. ',/A,' not specified.',/'Please correct the ', &
          A,' file.')
    
    ! Initialize the sum of the total volume fraction.
    SUM_EP = BC_EP_G(BCV)
    do M = 1, M_TOT
       
       ! If this phase is not present, clear out EPs and ROPs for the BC and
       ! cycle the solids loop. No need to continue checks.
       if(SKIP(M)) then
          BC_EP_S(BCV,M)  = ZERO
          BC_ROP_S(BCV,M) = ZERO
          cycle
       endif
       
       ! Set the solids density for the BC region.
       BC_ROs(M) = RO_s0(M)
       
       ! If both input parameters are defined. Make sure they are equivalent.
       if(IS_DEFINED(BC_ROP_S(BCV,M)) .and.                         &
            IS_DEFINED(BC_EP_S(BCV,M))) then
          
          if(.not.COMPARE(BC_EP_S(BCV,M)*BC_ROs(M),                  &
               BC_ROP_S(BCV,M))) then
             write(ERR_MSG,1214) BCV, trim(IFILE_NAME), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
1214      format('Error 1214: Illegal initial condition region : ',I3,/    &
               'BC_EP_s and BC_ROP_s are inconsistent. Please correct the ',/&
               A, ' file.')
          
          ! Compute BC_EP_s from BC_ROP_s
       elseif(IS_UNDEFINED(BC_EP_S(BCV,M))) then
          BC_EP_S(BCV,M) = BC_ROP_S(BCV,M) / BC_ROs(M)
          
          ! Compute BC_ROP_s from BC_EP_s and BC_ROs
       elseif(IS_UNDEFINED(BC_ROP_S(BCV,M))) then
          BC_ROP_S(BCV,M) = BC_EP_S(BCV,M) * BC_ROs(M)
          
       endif
       ! Add this phase to the total volume fraction.
       SUM_EP = SUM_EP + BC_EP_S(BCV,M)
    enddo
    
    ! Verify that the volume fractions sum to one.
    if(.not.COMPARE(SUM_EP,ONE)) then
       write(ERR_MSG,1215) BCV, trim(iVal(SUM_EP)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
1215 format('Error 1215: Illegal boundary condition region: ',I3,'. ',&
          'Sum of volume',/'fractions does NOT equal ONE. (SUM = ',A,   &
          ')',/'Please correct the ',A,' file.')
    
    call finl_err_msg

1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the ',A,' file.')

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
  subroutine check_bc_p_inflow(M_TOT, SKIP, BCV)

    use param    , only: dim_m
    use fld_const, only: ro_g0
    use bc,        only: bc_p_g, bc_rop_s, bc_u_g, bc_v_g, bc_w_g, &
                      &  bc_u_s, bc_v_s, bc_w_s

    integer, intent(in) :: BCV, M_TOT
    logical, intent(in) :: SKIP(DIM_M)
    integer             :: M

    CALL init_err_msg("CHECK_BC_P_INFLOW")

    ! Remove checks on bc_ep_g/bc_rop_s; using routine check_bc_outflow
    if (IS_UNDEFINED(BC_P_G(BCV))) then
       write(ERR_MSG,1000) 'BC_P_g', BCV, trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
1000   format('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the ',A,' file.')
       
    elseif (BC_P_G(BCV)<=ZERO .and. IS_UNDEFINED(RO_G0)) then
       write(ERR_MSG, 1101) BCV, trim(iVal(BC_P_G(BCV))), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
1101   format('Error 1101: Pressure must be greater than zero for ',    &
            'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/'Please ',  &
            'correct the ',A,' file.')
    endif
    
    ! Check that velocities are also specified. These are essentially used
    ! as initial conditions for the boundary region. If they are not
    ! specified then a default value is set here.
    if(IS_UNDEFINED(BC_U_G(BCV))) then
       BC_U_G(BCV) = ZERO
       write(ERR_MSG, 1300) trim(iVar('BC_U_g',BCV))
       call flush_err_msg(ABORT=.false.)
    endif
    
    if(IS_UNDEFINED(BC_V_G(BCV))) then
       BC_V_G(BCV) = ZERO
       write(ERR_MSG, 1300) trim(iVar('BC_V_g',BCV))
       call flush_err_msg(ABORT=.false.)
    endif
    
    if(IS_UNDEFINED(BC_W_G(BCV))) then
       BC_W_G(BCV) = ZERO
       write(ERR_MSG, 1300) trim(iVar('BC_W_g',BCV))
       call flush_err_msg(ABORT=.false.)
    endif
    
    do M = 1, M_TOT
       if (SKIP(M)) then
          BC_U_S(BCV,M) = ZERO
          BC_V_S(BCV,M) = ZERO
          BC_W_S(BCV,M) = ZERO
       else
          if(IS_UNDEFINED(BC_U_S(BCV,M))) then
             BC_U_S(BCV,M) = ZERO
             if(abs(BC_ROP_S(BCV,M)) > ZERO) then
                write(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M))
                call flush_err_msg(ABORT=.false.)
             endif
          endif
          
          if(IS_UNDEFINED(BC_V_S(BCV,M))) then
             BC_V_S(BCV,M) = ZERO
             if(abs(BC_ROP_S(BCV,M)) > ZERO) then
                write(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M))
                call flush_err_msg(ABORT=.false.)
             endif
          endif
          
          if(IS_UNDEFINED(BC_W_S(BCV,M))) then
             BC_W_S(BCV,M) = ZERO
             if(abs(BC_ROP_S(BCV,M)) > ZERO) then
                write(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M))
                call flush_err_msg(ABORT=.false.)
             endif
          endif
       endif
    enddo
    
1300 format('Warning 1300: ',A,' was undefined. This variable was ', &
          'set ',/ 'to zero to be used as the initial value in the BC ',&
          'region.')
    
    call finl_err_msg

  end subroutine check_bc_p_inflow

end module check_bc_inflow_module
