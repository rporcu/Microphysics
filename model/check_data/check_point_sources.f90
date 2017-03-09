module check_point_sources_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use get_ps_module,  only: get_ps
  use run,            only: IFILE_NAME
  use param1,         only: undefined, undefined_i, is_defined, is_undefined, zero
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg

  implicit none
  private 

  public check_point_sources


contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_POINT_SOURCES                                     !
  !  Author: J. Musser                                  Date: 10-JUN-13  !
  !                                                                      !
  !  Purpose: Check point source specifications.                         !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_point_sources(dx,dy,dz)

    use ps,    only: PS_DEFINED
    use param, only: DIMENSION_PS

    integer                  :: PSV
    real(c_real), intent(in) :: dx,dy,dz

    ! Initialize the error manager.
    call init_err_msg("CHECK_POINT_SOURCES")
    
    ! Determine which PSs are DEFINED
    call check_ps_geometry
    
    ! Loop over all PS arrays.
    do PSV = 1, DIMENSION_PS
       
       ! Verify user input for defined defined PS.
       if(PS_DEFINED(PSV)) then
          call get_ps(PSV,dx,dy,dz)
          call check_ps_gas_phase(PSV)
       else
          ! Verify that no data was defined for unspecified PS.
          call check_ps_overflow(PSV)
       endif
    enddo
    
    ! Clear the error manager.
    call finl_err_msg
    
    
  end subroutine check_point_sources


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_PS_GEOMETRY                                       !
  !  Author: J. Musser                                  Date: 10-JUN-13  !
  !                                                                      !
  !  Purpose: Check point source specifications.                         !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_ps_geometry


    use param, only: DIMENSION_PS
    use ps,    only: PS_DEFINED, POINT_SOURCE,       &
                   & PS_X_e, PS_X_w, PS_I_e, PS_I_w, &
                   & PS_Y_n, PS_Y_s, PS_J_n, PS_J_s, &
                   & PS_Z_t, PS_Z_b, PS_K_t, PS_K_b

    INTEGER :: PSV

    ! Initialize the error manager.
    CALL init_err_msg("CHECK_PS_GEOMETRY")

    ! Initialize the PS runtime flag.
    POINT_SOURCE = .false.
    
    ! Determine which point source indices have values.
    PSV_LP: do PSV = 1, DIMENSION_PS
       
       if (IS_DEFINED(PS_X_W(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_X_E(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_Y_S(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_Y_N(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_Z_B(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_Z_T(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_I_W(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_I_E(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_J_S(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_J_N(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_K_B(PSV)))   PS_DEFINED(PSV) = .true.
       if (IS_DEFINED(PS_K_T(PSV)))   PS_DEFINED(PSV) = .true.
       
       ! Skip consistency checks if nothing was defined.
       if (.not.PS_DEFINED(PSV)) cycle PSV_LP
       
       ! Flag that one or more point sources has been detected.
       POINT_SOURCE = .true.

       if(IS_UNDEFINED(PS_X_W(PSV)) .and. IS_UNDEFINED(PS_I_W(PSV))) then
          write(ERR_MSG,1101) PSV, 'PS_X_w and PS_I_w ', trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       if(IS_UNDEFINED(PS_X_E(PSV)) .and. IS_UNDEFINED(PS_I_E(PSV))) then
          write(ERR_MSG,1101) PSV, 'PS_X_e and PS_I_e ', trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       if(IS_UNDEFINED(PS_Y_S(PSV)) .and. IS_UNDEFINED(PS_J_S(PSV))) then
          write(ERR_MSG,1101) PSV, 'PS_Y_s and PS_J_s ', trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       if(IS_UNDEFINED(PS_Y_N(PSV)) .and. IS_UNDEFINED(PS_J_N(PSV))) then
          write(ERR_MSG,1101) PSV, 'PS_Y_n and PS_J_n ', trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       if(IS_UNDEFINED(PS_Z_B(PSV)) .and. IS_UNDEFINED(PS_K_B(PSV))) then
          write(ERR_MSG,1101) PSV, 'PS_Z_b and PS_K_b ', trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       if(IS_UNDEFINED(PS_Z_T(PSV)) .and. IS_UNDEFINED(PS_K_T(PSV))) then
          write(ERR_MSG,1101) PSV, 'PS_Z_t and PS_K_t ', trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       
1101   format('Error 1101: Point source ',I3,' is ill-defined.',/A,     &
            ' are not specified.',/'Please correct the ',A,' file.')

    enddo PSV_LP
    
    call finl_err_msg

  end subroutine check_ps_geometry


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_PS_GAS_PHASE                                      !
  !  Author: J. Musser                                  Date: 10-JUN-13  !
  !                                                                      !
  !  Purpose: Check point source specifications.                         !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_ps_gas_phase(PSV)

    use ps, only: PS_MASSFLOW_G, PS_U_g, PS_V_g, PS_W_g

    integer, intent(in) :: PSV

    ! Initialze the error manager.
    CALL init_err_msg("CHECK_PS_GAS_PHASE")

    ! Check mass flow and velocity
    if(IS_UNDEFINED(PS_MASSFLOW_G(PSV))) then
       if(IS_DEFINED(PS_U_g(PSV)) .or. &
            IS_DEFINED(PS_V_g(PSV)) .or. &
            IS_DEFINED(PS_W_g(PSV))) then
          
          write(ERR_MSG,1100) PSV, trim(iVar('PS_MASSFLOW_G',PSV)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
          
1100      format('Error 1100: Invalid specification for point source ',I3,&
               '.',/A,' is undefined but velocity is given.',/'Please ',    &
               'correct the ',A,' file.')
          
       else
          PS_MASSFLOW_G(PSV) = ZERO
          PS_U_g(PSV) = ZERO
          PS_V_g(PSV) = ZERO
          PS_W_g(PSV) = ZERO
       endif
       
    elseif(abs(PS_MASSFLOW_G(PSV)) < epsilon(ZERO)) then
       if(abs(PS_U_g(PSV)) < epsilon(ZERO) .or. &
            abs(PS_V_g(PSV)) < epsilon(ZERO) .or. &
            abs(PS_W_g(PSV)) < epsilon(ZERO)) then
          
          write(ERR_MSG,1101) PSV, trim(iVar('PS_MASSFLOW_G',PSV)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       
1101   format('Error 1101: Invalid specification for point source ',I3,&
            '.',/A,' is zero but velocity is given.',/'Please correct ', &
            'the ',A,' file.')
       
       ! Verify a physical mass flow
    elseif(PS_MASSFLOW_G(PSV) < ZERO) then
       write(ERR_MSG,1102) PSV, trim(iVar('PS_MASSFLOW_G',PSV)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
       
1102   format('Error 1102: Invalid specifications for point source ',I3,&
            '.',/A,' < 0.0. Point sources can only add mass to a system',/&
            'Please correct the ',A,' file.')
       
       
       ! Mass flow is specified:
    else
       
       ! Velocity does not have to be defined (no momentum source). If the
       ! components are UNDEFINED, zero them out.
       if(IS_UNDEFINED(PS_U_g(PSV))) PS_U_g(PSV) = ZERO
       if(IS_UNDEFINED(PS_V_g(PSV))) PS_V_g(PSV) = ZERO
       if(IS_UNDEFINED(PS_W_g(PSV))) PS_W_g(PSV) = ZERO
       
    endif
    
    call finl_err_msg

  end subroutine check_ps_gas_phase



  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_PS_OVERFLOW                                       !
  !  Author: J. Musser                                  Date: 10-JUN-13  !
  !                                                                      !
  !  Purpose: Check point source specifications.                         !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_ps_overflow(PSV)


    use ps, only: PS_MASSFLOW_G, PS_U_g, PS_V_g, PS_W_g

    integer, intent(in) :: PSV

    ! Initialize the error manager.
    call init_err_msg("CHECK_PS_OVERFLOW")
    
    
    if(IS_DEFINED(PS_MASSFLOW_G(PSV))) then
       write(ERR_MSG,1010) trim(iVar('PS_MASSFLOW_G',PSV))
       call flush_err_msg(ABORT=.true.)
    elseif(IS_DEFINED(PS_U_g(PSV))) then
       write(ERR_MSG,1010) trim(iVar('PS_U_g',PSV))
       call flush_err_msg(ABORT=.true.)
    elseif(IS_DEFINED(PS_V_g(PSV))) then
       write(ERR_MSG,1010) trim(iVar('PS_V_g',PSV))
       call flush_err_msg(ABORT=.true.)
    elseif(IS_DEFINED(PS_W_g(PSV))) then
       write(ERR_MSG,1010) trim(iVar('PS_W_g',PSV))
       call flush_err_msg(ABORT=.true.)
    endif
    
    call FINL_ERR_MSG
    
1010 format('Error 1010: ',A,' specified in an undefined PS region.')
    
  end subroutine check_ps_overflow

end module check_point_sources_module
