module check_bc_outflow_module

  use bl_fort_module, only: c_real
  use iso_c_binding , only: c_int
  use run,            only: IFILE_NAME
  use param1,         only: one, undefined, zero, is_undefined, is_defined, equal
  use error_manager,  only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival


  implicit none
  private 

  public check_bc_outflow
  public check_bc_p_outflow
  public check_bc_mass_outflow

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: CHECK_BC_OUTFLOW                                         !
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Provided a detailed error message concerning specification  !
  ! of bc_ep_g + bc_rop_s at a outflow boundary (and pressure inflow)    !
  !                                                                      !
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_outflow(M_TOT, BCV)

    use bc,       only: bc_ep_g, bc_rop_s
    use constant, only: ro_s0
    use run,      only: solids_model
    use toleranc, only: compare

    integer, intent(in) :: BCV, M_TOT
    integer             :: M
    real(c_real)        :: SUM_EP
    logical             :: FLAG_WARNING

    FLAG_WARNING = .true.

    CALL init_err_msg("CHECK_BC_OUTFLOW")

    ! if bc_ep_g is defined at the outflow boundary, then the sum of ep_g
    ! and ep_s at the boundary may not equal one given the code in the
    ! subroutine set_outflow (see code for details).
    ! therefore if bc_ep_g and/or bc_rop_s are defined, perform possible
    ! data consistency checks and, when appropriate, provide the user with
    ! a warning about their chosen settings.
    if (IS_DEFINED(BC_EP_G(BCV))) then

       SUM_EP = BC_EP_G(BCV)
       do M = 1, M_TOT
          
          if(SOLIDS_MODEL(M) /= 'TFM' .and. FLAG_WARNING) then
             write(ERR_MSG, 1101) trim(iVar('BC_EP_g',BCV))
             call flush_err_msg
             FLAG_WARNING = .false.
          endif
          
          if(IS_UNDEFINED(BC_ROP_S(BCV,M))) then
             
             if(EQUAL(BC_EP_G(BCV), ONE)) then
                ! what does it mean to force the bulk density to zero at the
                ! boundary? (does this value matter anyway?)
                BC_ROP_S(BCV,M) = ZERO
             elseif(M_TOT == 1 ) then
                BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S0(M)
             else
                ! bc_ep_g is defined but some bc_rop_s(m) are undefined.
                ! in this case, ep_p in the outflow boundary will be based on the user
                ! defined value of bc_ep_g, while rop_s would become based on the
                ! value in the adjacent fluid cell. consequently, no check ensures
                ! the result is consistent with a requirement for ep_g+ep_s=1.
                write(ERR_MSG, 1102) trim(iVar('BC_EP_g',BCV))
                call flush_err_msg(ABORT=.true.)
1102            format('Warning 1102: Volume fraction may not sum to one when ',/&
                     A,' is defined.')
             endif
          endif  ! end if(bc_rop_s(bcv,m) == undefined)
          
          ! by this point bc_rop_s should either be defined or mfix exited
          ! therefore we can check that sum of void fraction and solids volume
          ! fractions
          SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S0(M)
       enddo
       
       ! now verify that the volume fractions sum to one.
       if(.not.COMPARE(SUM_EP,ONE)) then
          write(ERR_MSG,1103) BCV, trim(iVal(SUM_EP)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       
       ! bc_ep_g is not defined but check if any bc_rop_s are defined
    else
       
       SUM_EP = ZERO
       do M = 1, M_TOT
          if(IS_DEFINED(BC_ROP_S(BCV,M))) then
             if(SOLIDS_MODEL(M) /= 'TFM') then
                write(ERR_MSG, 1101) trim(iVar('BC_ROP_s',BCV,M))
                call flush_err_msg
             endif
             SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S0(M)
          endif
       enddo
       
       ! verify that the sum of any specified volume fractions is not greater
       ! than one
       if(SUM_EP > ONE) then
          write(ERR_MSG,1103) BCV, trim(iVal(SUM_EP)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
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
          ')',/'Please correct the ',A,' file.')
    
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

    if (IS_UNDEFINED(BC_P_G(BCV))) then
       write(ERR_MSG,1000) trim(iVar('BC_P_g',BCV)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
       
    elseif (BC_P_G(BCV)<=ZERO .and. IS_UNDEFINED(RO_G0)) then
       write(ERR_MSG, 1100) BCV, trim(iVal(BC_P_G(BCV))), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
    
1100 format('Error 1100: Pressure must be greater than zero for ',    &
          'compressible flow',/3x,'BC_P_g(',I3,') = ',A,/'Please ',     &
          'correct the ',A,' file.')
    
    ! Clean up and return.
    call finl_err_msg

1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the ',A,' file.')

  end subroutine check_bc_p_outflow
  

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: CHECK_BC_MASS_OUTFLOW                                    !
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Provided a detailed error message when the sum of volume    !
  !                                                                      !
  ! Comments:                                                            !
  !     The velocities at the outflow face are fixed and the momentum    !
  !     equations are not solved in the outflow cells. Since the flow    !
  !     is out of the domain none of the other scalars should need to    !
  !     be specified (e.g., mass fractions, void fraction, etc.,).       !
  !     Such values will become defined according to their adjacent      !
  !     fluid cell                                                       !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_bc_mass_outflow(M_TOT, BCV)

    use fld_const, only: ro_g0
    use bc,        only: bc_plane, bc_dt_0, bc_massflow_g, bc_volflow_g, &
                       & bc_massflow_s, bc_volflow_s, bc_ep_g, bc_rop_s, &
                       & bc_p_g, bc_t_g, bc_u_g, bc_v_g, bc_w_g

    integer, intent(in) :: BCV, M_TOT
    integer             :: M


    call init_err_msg("CHECK_BC_MASS_OUTFLOW")

    if(IS_UNDEFINED(BC_DT_0(BCV))) then
       write(ERR_MSG, 1000) trim(iVar('BC_DT_0',BCV)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
    
    if(IS_DEFINED(BC_MASSFLOW_G(BCV)) .or. &
         IS_DEFINED(BC_VOLFLOW_G(BCV))) then
       if (IS_UNDEFINED(BC_EP_G(BCV))) then
          write(ERR_MSG,1101) trim(iVar('BC_EP_G',BCV)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
1101   format('Error 1101: Invalid mass outflow boundary condition: ', /&
            'BC_MASSFLOW_G and/or BC_VOLFLOW_G are DEFINED but ',&
            A,' is not ',/'Please correct the ',A,' file.')
    endif
    
    do M = 1, M_TOT
       if(IS_DEFINED(BC_MASSFLOW_S(BCV,M)) .or. &
            IS_DEFINED(BC_VOLFLOW_S(BCV,M))) then
          write(ERR_MSG,1102) trim(iVar('BC_MASSFLOW_S',BCV,M)), &
               trim(iVar('BC_VOLFLOW_S',BCV,M))
1102      format('Warning 1102: ', A,' and/or ', A,' have been defined',/&
               'at a mass outflow boundary. A specified solids flow ',&
               'rate may not be ',/'physically achievable depending on the ',&
               'system and simulation ',/'setup.')
          
          if (IS_UNDEFINED(BC_ROP_S(BCV,M))) then
             write(ERR_MSG,1103) trim(iVar('BC_ROP_S',BCV,M)), trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
1103      format('Error 1103: Invalid mass outflow boundary condition: ', /&
               'BC_MASSFLOW_S and/or BC_VOLFLOW_S are DEFINED but ',&
               A,' is not ',/'Please correct the ',A,' file.')
          
       endif
    enddo
    
    ! This check probably needs changed.
    if(IS_UNDEFINED(RO_G0) .and. (IS_UNDEFINED(BC_P_G(BCV)) .or.       &
         IS_UNDEFINED(BC_T_G(BCV))) .and.abs(BC_MASSFLOW_G(BCV)) > ZERO) then
       
       if(BC_PLANE(BCV)=='W' .or. BC_PLANE(BCV)=='E') then
          if(abs(BC_U_G(BCV)) > ZERO) then
             write(ERR_MSG, 1100) BCV, 'BC_U_g', trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
       elseif(BC_PLANE(BCV)=='N' .or. BC_PLANE(BCV)=='S') then
          if(abs(BC_V_G(BCV)) > ZERO) then
             write(ERR_MSG, 1100) BCV, 'BC_V_g', trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
       elseif (BC_PLANE(BCV)=='T' .or. BC_PLANE(BCV)=='B') then
          if(abs(BC_W_G(BCV)) > ZERO) then
             write(ERR_MSG, 1100)  BCV, 'BC_W_g', trim(IFILE_NAME)
             call flush_err_msg(ABORT=.true.)
          endif
       endif
    endif   ! end if/else (ro_g0 /=undefined)
    
1100 format('Error 1100: Invalid mass outflow boundary condition: ',  &
          I3,/'RO_g0, BC_P_g, and BC_T_g are UNDEFINED and ',A,' is ',  &
          'non-zero',/'Please correct the ',A,' file.')    
    
    call finl_err_msg
    
1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
          'correct the ',A,' file.')
    
  end subroutine check_bc_mass_outflow


end module check_bc_outflow_module
