  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: check_bc_flow                                             !
  !  Author: P. Nicoletti                               Date: 10-DEC-91  !
  !                                                                      !
  !  Purpose: Check boundary condition specifications                    !
  !     - convert physical locations to i, j, k's                        !
  !     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
  !     - check specification of physical quantities                     !
  !                                                                      !
  !  Comments:                                                           !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

  subroutine check_bc_flow() &
    bind(C,name ="check_bc_flow")

    use constant,      only: mmax
    use param1,        only: zero, one, undefined, is_undefined, is_defined, equal
    use param,         only: dimension_bc, dim_m
    use run,           only: IFILE_NAME
    use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, err_msg
    use bc,            only: bc_defined, bc_type, bc_rop_s, bc_ep_s, bc_ep_g, &
                           & bc_u_g, bc_v_g, bc_w_g, bc_u_s, bc_v_s, bc_w_s,  &
                           & bc_plane

    implicit none

    integer :: bcv, i, mmax_tot
    logical :: skip(1:dim_m)     ! Flag to skip checks on indexed solid phase.

    ! Initialize the error manager.
    call init_err_msg("SET_BC_FLOW")

    ! Total number of solids.
    mmax_tot = mmax

    ! Loop over each defined BC and check the user data.
    do bcv = 1, dimension_bc

       if(.not.bc_defined(bcv)) cycle

       ! Determine which solids phases are present.
       skip = .false.
       do i = 1, dim_m
          if ((equal(bc_rop_s(bcv,i), undefined).or.equal(bc_rop_s(bcv,i), zero)) &
               .and.(equal(bc_ep_s(bcv,i), undefined).or.equal(bc_ep_s(bcv,i), zero))) then
             skip = .true.
          endif
       end do

       if(mmax_tot == 1 .and. .not.equal(bc_ep_g(bcv), one)) skip(1) = .false.

       select case (trim(BC_TYPE(BCV)))

       case ('MASS_INFLOW','MI')
          call check_bc_vel_inflow(mmax_tot, skip, bcv)

       case ('MASS_OUTFLOW','MO')
          call check_bc_vel_outflow(mmax_tot, skip, bcv)
       end select
    enddo

    ! Cleanup and exit.
    call finl_err_msg

    contains


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: check_bc_vel_inflow                                      !
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Provided a detailed error message when the sum of volume    !
  !                                                                      !
  ! Comments:                                                            !
  !     The velocities at the inflow face are fixed and the momentum     !
  !     equations are not solved in the inflow cells. Since the flow is  !
  !     into the domain all other scalars that are used need to be       !
  !     specified (e.g., mass fractions, void fraction, etc.,)           !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_bc_vel_inflow(m_tot, skip, bcv)

    use param, only: dim_m

    integer, intent(in)           :: bcv, m_tot
    logical, intent(in)           :: skip(dim_m)
    integer                       :: m
    character(len=:), allocatable :: fmt1, fmt2


    ! Define format for error messages
    fmt1 = "('Error 1300: Invalid flow direction.'," //&
          & "A,' should be ', A,' zero. ',/"               //&
          & "'Please correct the "//trim(IFILE_NAME)//" file.')"

    fmt2 = "('Error 1000: Required input not specified: ',A,/'Please ',"//&
         & "'correct the "//trim(IFILE_NAME)//" file.')"

    call init_err_msg("CHECK_BC_VEL_INFLOW")


    ! Check that gas phase velocities are defined.
    if(IS_UNDEFINED(BC_U_G(BCV))) then
       write(err_msg,fmt2) trim(iVar('BC_U_g',BCV))
       call flush_err_msg(ABORT=.true.)
    endif

    if (IS_UNDEFINED(BC_V_G(BCV))) then
       write(err_msg,fmt2) trim(iVar('BC_V_g',BCV))
       call flush_err_msg(ABORT=.true.)
    endif

    if(IS_UNDEFINED(BC_W_G(BCV))) then
       write(err_msg,fmt2) trim(iVar('BC_W_g',BCV))
       call flush_err_msg(ABORT=.true.)
    endif

    ! Check that solids phase velocities are defined.
    do M = 1, M_TOT
       if(IS_UNDEFINED(BC_U_S(BCV,M))) then
          if(SKIP(M)) then
             BC_U_S(BCV,M) = ZERO
          else
             write(err_msg,fmt2) trim(iVar('BC_U_s',BCV,M))
             call flush_err_msg(ABORT=.true.)
          endif
       endif

       if(IS_UNDEFINED(BC_V_S(BCV,M))) then
          if(SKIP(M)) then
             BC_V_S(BCV,M) = ZERO
          else
             WRITE(err_msg,fmt2) trim(iVar('BC_V_s',BCV,M))
             call flush_err_msg(ABORT=.true.)
          endif
       endif

       if(IS_UNDEFINED(BC_W_S(BCV,M))) then
          if(SKIP(M)) then
             BC_W_S(BCV,M) = ZERO
          else
             write(err_msg,fmt2) trim(iVar('BC_W_s',BCV,M))
             call flush_err_msg(ABORT=.true.)
          endif
       endif
    enddo

    ! Check that gas phase velocities are consistent.
    select case (bc_plane(bcv))

    case ('W')
       if(BC_U_G(BCV) > ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_U_g',BCV)), '<'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_U_S(BCV,M) > ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_U_s',BCV,M)), '<'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('E')
       if(BC_U_G(BCV) < ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_U_g',BCV)), '>'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_U_S(BCV,M) < ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_U_s',BCV,M)), '>'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('S')
       if(BC_V_G(BCV) > ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_V_g',BCV)), '<'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_V_S(BCV,M) > ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_V_s',BCV,M)), '<'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('N')
       if(BC_V_G(BCV) < ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_V_g',BCV)), '>'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_V_S(BCV,M) < ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_V_s',BCV,M)), '>'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('B')
       if(BC_W_G(BCV) > ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_W_g',BCV)), '<'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_W_S(BCV,M) > ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_W_s',BCV,M)), '<'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('T')
       if(BC_W_G(BCV) < ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_W_g',BCV)), '>'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_W_S(BCV,M) < ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_W_s',BCV,M)), '>'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    end select

    call finl_err_msg

    deallocate( fmt1, fmt2 )

  end subroutine check_bc_vel_inflow



  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: check_bc_vel_outflow                                     !
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
  subroutine check_bc_vel_outflow(m_tot, skip, bcv)

    use bc, only: dim_m

    ! loop/variable indices
    integer, intent(in)           :: bcv, m_tot
    logical, intent(in)           :: skip(dim_m)
    integer                       :: m
    character(len=:), allocatable :: fmt1, fmt2


    ! Set format for error messages
    fmt1 = "('Error 1300: Invalid flow direction. ',A,' should be ',"// &
         & " A,' zero. ',/'Please correct the "//trim(IFILE_NAME)//" file.')"
    fmt2 = " ('Error 1000: Required input not specified: ',A,/'Please ',"// &
         & " 'correct the "//trim(IFILE_NAME)//" file.')"

    call init_err_msg("CHECK_BC_VEL_OUTFLOW")

    ! Check that gas phase velocities are defined.
    if(IS_UNDEFINED(BC_U_G(BCV))) then
       write(err_msg,fmt2) trim(iVar('BC_U_g',BCV))
       call flush_err_msg(ABORT=.true.)
    endif

    if (IS_UNDEFINED(BC_V_G(BCV))) then
       write(err_msg,fmt2) trim(iVar('BC_V_g',BCV))
       call flush_err_msg(ABORT=.true.)
    endif

    if(IS_UNDEFINED(BC_W_G(BCV))) then
       write(err_msg,fmt2) trim(iVar('BC_W_g',BCV))
       call flush_err_msg(ABORT=.true.)
    endif

    ! Check that solids phase velocities are defined.
    do M = 1, M_TOT
       if(IS_UNDEFINED(BC_U_S(BCV,M))) then
          if(SKIP(M)) then
             BC_U_S(BCV,M) = ZERO
          else
             write(err_msg,fmt2) trim(iVar('BC_U_s',BCV,M))
             call flush_err_msg(ABORT=.true.)
          endif
       endif

       if(IS_UNDEFINED(BC_V_S(BCV,M))) then
          if(SKIP(M)) then
             BC_V_S(BCV,M) = ZERO
          else
             write(err_msg,fmt2) trim(iVar('BC_V_s',BCV,M))
             call flush_err_msg(ABORT=.true.)
          endif
       endif

       if(IS_UNDEFINED(BC_W_S(BCV,M))) then
          if(SKIP(M)) then
             BC_W_S(BCV,M) = ZERO
          else
             write(err_msg,fmt2) trim(iVar('BC_W_s',BCV,M))
             call flush_err_msg(ABORT=.true.)
          endif
       endif
    enddo


    ! Check that gas phase velocities are consistent.
    select case (bc_plane(BCV))

    case ('W')
       if(BC_U_G(BCV) < ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_U_g',BCV)), '>'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_U_S(BCV,M) < ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_U_s',BCV,M)), '>'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('E')
       if(BC_U_G(BCV) > ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_U_g',BCV)), '<'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_U_S(BCV,M) > ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_U_s',BCV,M)), '<'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('S')
       if(BC_V_G(BCV) < ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_V_g',BCV)), '>'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_V_S(BCV,M) < ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_V_s',BCV,M)), '>'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('N')
       if(BC_V_G(BCV) > ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_V_g',BCV)), '<'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_V_S(BCV,M) > ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_V_s',BCV,M)), '<'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('B')
       if(BC_W_G(BCV) < ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_W_g',BCV)), '>'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_W_S(BCV,M) < ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_W_s',BCV,M)), '>'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('T')
       if(BC_W_G(BCV) > ZERO) then
          write(err_msg,fmt1) trim(iVar('BC_W_g',BCV)), '<'
          call flush_err_msg
       ENDIF
       do M = 1, M_TOT
          if(BC_W_S(BCV,M) > ZERO) then
             write(err_msg, fmt1) trim(iVar('BC_W_s',BCV,M)), '<'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    end select


    call finl_err_msg

    deallocate(fmt1, fmt2)

  end subroutine check_bc_vel_outflow

  end subroutine check_bc_flow
