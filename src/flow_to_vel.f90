module flow_to_vel_new_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use exit_mod,       only: mfix_exit
  use toleranc,       only: compare
  use param,          only: DIM_M
  use run,            only: IFILE_NAME
  use param1,         only: UNDEFINED, ZERO, ONE, IS_DEFINED, IS_UNDEFINED
  use error_manager,  only: finl_err_msg, err_msg, flush_err_msg, &
                          &init_err_msg, ivar

  use bc, only: BC_MASSFLOW_G,  BC_MASSFLOW_S, BC_VOLFLOW_G, BC_VOLFLOW_S, &
              & bc_type, bc_plane, bc_u_s, bc_v_s, bc_w_s, bc_ep_s,        &
              & bc_ep_g, bc_area, bc_u_g, bc_v_g, bc_w_g

  implicit none
  private

  public flow_to_vel_new

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: FLOW_TO_VEL_NEW                                         !
  !  Author: M. Syamlal                                 Date: 28-JUL-92  !
  !                                                                      !
  !  Purpose: Convert volumetric and mass flow rates to velocities       !
  !     A specified mass flow rate is first converted to volumetric      !
  !     flow rate. The volumetric flow rate is then converted to a       !
  !     velocity.                                                        !
  !                                                                      !
  !    When both flow rates and velocities are specified, a consistency  !
  !    check is done. The first time flow_to_vel is called in by setting !
  !    the logical DO_VEL_CHECK to .TRUE.. If cut-cells are not used,    !
  !    flow_to_vel is only called once.  When cut-cells are used,        !
  !    flow_to_vel is called another time after the cut-cell pre-        !
  !    processing stage. During, the second call, the velocity check     !
  !    should not be performed, because the velocity assigned suring the !
  !    first call will not match the flow rate. Therfore, when called    !
  !    from cut_cell_preprocessing.f DO_VEL_CHECK is set to .FALSE..     !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine flow_to_vel_new(do_vel_check, m_tot, skip, bcv)

    logical, intent(in) :: do_vel_check, skip(dim_m)
    integer, intent(in) :: m_tot, bcv
    logical             :: converted = .false.
    integer             :: m

    call init_err_msg("FLOW_TO_VEL_NEW")

    ! mass flows rates are converted to volumetric flow rates.
    if(is_defined(bc_massflow_g(bcv))) &
         call gas_massflow_to_volflow(bcv)

    do m=1,m_tot
       if(is_defined(bc_massflow_s(bcv,m))) &
            call solids_massflow_to_volflow(bcv,m,skip(m))
    enddo

    ! volumetric flow rates are converted to velocities.
    if(is_defined(bc_volflow_g(bcv))) then
       call gas_volflow_to_velocity(do_vel_check, bcv)
       ! set the conversion flag.
       converted = .true.
    endif

    do m=1,m_tot
       if(is_defined(bc_volflow_s(bcv,m))) then
          call solids_volflow_to_velocity(do_vel_check,bcv,m,skip(m))
          ! set the conversion flag.
          converted = .true.
       endif
    enddo

    if(converted ) then
       write(err_msg, 1100)
       call flush_err_msg
    endif

    call finl_err_msg


1100 format('Warning 1100: Some volumetric or mass flow rates have ', &
         'been converted',/'velocity. Ensure that the third (unused) ',&
         'dimension in 2D simulations',/'is correctly specified.')

  end subroutine flow_to_vel_new



  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: GAS_MASSFLOW_TO_VOLFLOW                                 !
  !  Author: M. Syamlal                                 Date: 28-JUL-92  !
  !                                                                      !
  !  Purpose: Convert a gas phase BC input from a mass flow rate to      !
  !  a volumetric flow rate.                                             !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine gas_massflow_to_volflow(bcv)

    use bc,        only: bc_massflow_g, bc_p_g, bc_t_g, bc_volflow_g
    use eos,       only: eosg
    use fld_const, only: mw_avg, ro_g0
    use scales   , only: p_ref

    integer, intent(in) :: bcv
    real(c_real)        :: volflow
    real(c_real)        :: mw

    call init_err_msg("GAS_MASSFLOW_TO_VOLFLOW")

    ! No need to convert if the mass flow is zero.
    if(compare(bc_massflow_g(bcv),zero)) then
       volflow = zero

       ! incompressible gas bc.
    elseif(is_defined(ro_g0)) then
       volflow = bc_massflow_g(bcv)/ro_g0

       ! well-defined compresible gas bc.
    elseif(is_defined(bc_p_g(bcv)) .and. is_defined(bc_t_g(bcv))) then
       mw = mw_avg
       volflow = bc_massflow_g(bcv) / &
            eosg(mw,(bc_p_g(bcv)-p_ref),bc_t_g(bcv))

       ! fails. this shouldn't happen as previous checks should catch any
       ! errors leading to this routine.
    else
       write(err_msg, 1100) bcv
       call flush_err_msg(abort=.true.)

1100   format('Error 1100: Boundary condition ',I3,' failed sanity ',   &
            'check.',/'Please report this to the MFIX mailing list.')

    endif

    ! Check that a specified volumetric flow matches the calculated value.
    if(is_defined(bc_volflow_g(bcv))) then
       if(.not.compare(volflow,bc_volflow_g(bcv))) then
          write(err_msg,1101) trim(ivar('BC_MASSFLOW_g',bcv)), bcv,  &
               volflow, bc_volflow_g(bcv), trim(IFILE_NAME)
          call flush_err_msg(abort=.true.)
       endif
    else

       ! store the calculated volumetric flow rate.
       bc_volflow_g(bcv) = volflow
    endif

1101 FORMAT('Error 1101: Volumetric flow rate calculated from ',A,/   &
         'does NOT equal the specified volumetric flow rate for BC',I3,&
         /3x,'>>> Calculated: ',G14.7,/3x,'>>> Specified:  ',G14.7,/   &
         'Please correct the ',A,' file.')

    ! clean up and return.
    call finl_err_msg

  end subroutine gas_massflow_to_volflow



  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: SOLIDS_MASSFLOW_TO_VOLFLOW                              !
  !  Author: M. Syamlal                                 Date: 28-JUL-92  !
  !                                                                      !
  !  Purpose: Convert solids phase BC input from a mass flow rate to     !
  !  a volumetric flow rate.                                             !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine solids_massflow_to_volflow(bcv,m, skip_m)

    use bc,       only: bc_massflow_s, bc_volflow_s
    use constant, only: ro_s0

    integer, intent(in) :: bcv, m
    logical, intent(in) :: skip_m
    real(c_real)        :: volflow

    call init_err_msg("SOLIDS_MASSFLOW_TO_VOLFLOW")


    if(skip_m) then
       write(err_msg,1100) m, bcv, trim(ivar("BC_MASSFLOW_S",bcv,m)), &
            trim(IFILE_NAME)
       call flush_err_msg(abort=.true.)
    endif

1100 FORMAT('Error 1100: Solids phase ',I2,' has a specified mass ',  &
         'flow rate',/'at BC ',I3,', ',A,'. But, both BC_ROP_s and ',&
         'BC_EP_s are zero or undefined.',/'Please correct the ',&
         A,' file.')

    VOLFLOW = BC_MASSFLOW_S(BCV,M)/RO_S0(M)

    ! If volumetric flow is also specified compare both
    if(is_defined(bc_volflow_s(bcv,m))) then
       if(.not.compare(volflow,bc_volflow_s(bcv,m))) then
          write(err_msg,1101) trim(ivar('BC_MASSFLOW_S',bcv,m)), bcv, &
               volflow, bc_volflow_s(bcv,m),  trim(IFILE_NAME)
          call flush_err_msg(abort=.true.)
       endif
    else
       bc_volflow_s(bcv,m) = volflow
    endif

    call finl_err_msg

1101 format('Error 1101: Volumetric flow rate calculated from ',A,/   &
         'does NOT equal the specified volumetric flow rate for BC',I3,&
         /3x,'>>> Calculated: ',G14.7,/3x,'>>> Specified:  ',G14.7,/   &
         'Please correct the ',A,' file.')


  end subroutine solids_massflow_to_volflow




  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: GAS_VOLFLOW_TO_VELOCITY                                 !
  !  Author: M. Syamlal                                 Date: 28-JUL-92  !
  !                                                                      !
  !  Purpose: Convert gas phase volumetric rate to a velocity.           !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine gas_volflow_to_velocity(do_vel_check, bcv)

    use compar, only: mype

    integer, intent(in) :: bcv
    logical, intent(in) :: do_vel_check
    real(c_real)        :: sgn, off, vel

    call init_err_msg("GAS_VOLFLOW_TO_VELOCITY")

    select case (trim(bc_type(bcv)))
    case ('MASS_INFLOW');  SGN =  ONE; OFF = ZERO
    case ('MASS_OUTFLOW'); SGN = -ONE; OFF = ONE
    case DEFAULT
       write(*,*) 'error in GAS_VOLFLOW_TO_VELOCITY'
       call mfix_exit(myPE)
    end select

    select case (bc_plane(BCV))
    case ('W'); SGN = -SGN
    case ('S'); SGN = -SGN
    case ('B'); SGN = -SGN
    end select

    ! Calculate the velocity based on the volumetric flow rate,
    ! BC area and BC volume fraction.
    VEL = SGN*BC_VOLFLOW_G(BCV)/(BC_AREA(BCV)*BC_EP_G(BCV))

    ! if the user also defined the boundary velocity through the plane, then
    ! check that the calculated value agrees with the specified value. if
    ! the user did not define the boundary velocity through the plane, then
    ! if mass_inflow set the value of the boundary velocity to the
    ! calculated value. otherwise do nothing.
    if(bc_plane(BCV) == 'W' .or. bc_plane(BCV)== 'E') then

       if(is_defined(bc_u_g(bcv)) .and. do_vel_check) then
          if(.not.compare(vel,bc_u_g(bcv))) then
             write(err_msg,1100) bcv, vel, 'BC_U_g', bc_u_g(bcv)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_u_g(bcv) = vel
          bc_v_g(bcv) = off * bc_v_g(bcv)
          bc_w_g(bcv) = off * bc_w_g(bcv)
       endif

    elseif(bc_plane(BCV) == 'S' .or. bc_plane(BCV)== 'N') then
       if(is_defined(bc_v_g(bcv)) .and. do_vel_check) then
          if(.not.compare(vel,bc_v_g(bcv))) then
             write(err_msg, 1100) bcv, vel, 'BC_V_g', bc_v_g(bcv)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_v_g(bcv) = vel
          bc_u_g(bcv) = off * bc_u_g(bcv)
          bc_w_g(bcv) = off * bc_w_g(bcv)
       endif

    elseif(bc_plane(BCV) == 'B' .or. bc_plane(BCV)== 'T') then
       if(is_defined(bc_w_g(bcv)) .and. do_vel_check) then
          if(.not.compare(vel, bc_w_g(bcv))) then
             write(err_msg, 1100) bcv, vel, 'BC_W_g', bc_w_g(bcv)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_w_g(bcv) = vel
          bc_u_g(bcv) = off * bc_u_g(bcv)
          bc_v_g(bcv) = off * bc_v_g(bcv)
       endif

    endif

    call finl_err_msg


1100 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,') = ',G14.7,/1X,70('*')/)

  end subroutine gas_volflow_to_velocity

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: SOLIDS_VOLFLOW_TO_VELOCITY                              !
  !  Author: M. Syamlal                                 Date: 28-JUL-92  !
  !                                                                      !
  !  Purpose: Convert volumetric and mass flow rates to velocities       !
  !     A specified mass flow rate is first converted to volumetric      !
  !     flow rate. The volumetric flow rate is then converted to a       !
  !     velocity.                                                        !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine solids_volflow_to_velocity(do_vel_check, bcv, m, skip_m)

    use compar, only: myPE

    integer, intent(in) :: bcv, m
    logical, intent(in) :: do_vel_check, skip_m
    real(c_real)        :: vel, sgn, off


    call init_err_msg("SOLIDS_VOLFLOW_TO_VELOCITY")

    if(skip_m) then
       write(err_msg,1100) m, bcv, trim(ivar("BC_VOLFLOW_S",bcv,m)), &
            & trim(IFILE_NAME)
       call flush_err_msg(abort=.true.)
    endif

1100 format('Error 1100: Solids phase ',I2,' has a specified ',       &
         'volumetric flow rate',/'at BC ',I3,', ',A,'. But, both ',&
         'BC_ROP_s and BC_EP_s are zero or undefined.',/'Please ',&
         'the ',A,' file.')

    select case (trim(BC_TYPE(BCV)))
    case ('MASS_INFLOW');  SGN =  ONE; OFF = ZERO
    case ('MASS_OUTFLOW'); SGN = -ONE; OFF = ONE
    case DEFAULT
       write(*,*) 'error in SOLIDS_VOLFLOW_TO_VELOCITY'
       call mfix_exit(myPE)
    end select

    select case (bc_plane(BCV))
    CASE ('W'); SGN = -SGN
    case ('S'); SGN = -SGN
    case ('B'); SGN = -SGN
    end select

    if(abs(bc_ep_s(bcv,m)) > zero) then
       vel = sgn * bc_volflow_s(bcv,m)/(bc_area(bcv)*bc_ep_s(bcv,m))
    else
       if(abs(bc_volflow_s(bcv,m)) > zero) then
          vel = zero
       else
          call flush_err_msg(abort=.true.)
       endif
    endif


    if(bc_plane(BCV) == 'W' .or. bc_plane(BCV)== 'E') then
       if(is_defined(bc_u_s(bcv,m)) .and. do_vel_check) then
          if(.not.compare(vel, bc_u_s(bcv,m))) then
             write(err_msg, 1300) bcv, (-vel), 'BC_u_s', m, bc_u_s(bcv,m)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_u_s(bcv,m) = vel
          bc_v_s(bcv,m) = off * bc_v_s(bcv,m)
          bc_w_s(bcv,m) = off * bc_w_s(bcv,m)
       endif

    elseif(bc_plane(BCV) == 'S' .or. bc_plane(BCV)== 'N') then
       if(is_defined(bc_v_s(bcv,m)) .and. do_vel_check) then
          if(.not.compare(vel,bc_v_s(bcv,m))) then
             write(err_msg,1300) bcv, vel, 'BC_v_s', m, bc_v_s(bcv,m)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_v_s(bcv,m) = vel
          bc_u_s(bcv,m) = off * bc_u_s(bcv,m)
          bc_w_s(bcv,m) = off * bc_w_s(bcv,m)
       endif

    elseif(bc_plane(BCV) == 'B' .or. bc_plane(BCV)== 'T') then
       if(is_defined(bc_w_s(bcv,m)) .and. do_vel_check) then
          if(.not.compare(vel,bc_w_s(bcv,m))) then
             write(err_msg, 1300) bcv, vel, 'BC_w_s', m, bc_w_s(bcv,m)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_w_s(bcv,m) = vel
          bc_u_s(bcv,m) = off * bc_u_s(bcv,m)
          bc_v_s(bcv,m) = off * bc_v_s(bcv,m)
       endif
    endif


    call finl_err_msg



1300 format(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,I1,') = ',G14.7,/1X,70('*')/)

  end subroutine solids_volflow_to_velocity

end module flow_to_vel_new_module
