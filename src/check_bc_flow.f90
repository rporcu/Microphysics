  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: check_bc_flow                                           !
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
    use bc,            only: bc_defined, bc_type, bc_ep_s, bc_ep_g, &
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

       if(bc_defined(bcv)) then

         ! Determine which solids phases are present.
          do i = 1, dim_m
             skip(i) = (equal(bc_ep_s(bcv,i), zero))
          end do

          select case (trim(bc_type(bcv)))
          case ('MASS_INFLOW','MI')
             call check_bc_vel_inflow(mmax_tot, skip, bcv)

          case ('MASS_OUTFLOW','MO')
             call check_bc_vel_outflow(mmax_tot, skip, bcv)
          end select
       endif
    enddo

    ! Cleanup and exit.
    call finl_err_msg

    contains


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: check_bc_vel_inflow                                      !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_bc_vel_inflow(m_tot, skip, bcv)

    use param, only: dim_m

    integer, intent(in)           :: bcv, m_tot
    logical, intent(in)           :: skip(dim_m)
    integer                       :: m

    ! Define format for error messages

    call init_err_msg("CHECK_BC_VEL_INFLOW")


1000 format('Error 1000: Required input not specified: ',A,/&
        'Please correct the input deck.')

    ! Check that gas phase velocities are consistent.
    select case (bc_plane(bcv))

    case ('W')
       if(bc_u_g(bcv) > zero) then
          write(err_msg,1300) trim(ivar('BC_U_g',bcv)), '<'
          call flush_err_msg
       endif
       do m = 1, m_tot
          if(bc_u_s(bcv,m) > zero) then
             write(err_msg, 1300) trim(ivar('BC_U_s',bcv,m)), '<'
             call flush_err_msg(abort=.true.)
          endif
       enddo

    case('E')
       if(bc_u_g(bcv) < zero) then
          write(err_msg,1300) trim(ivar('BC_U_g',bcv)), '>'
          call flush_err_msg
       endif
       do m = 1, m_tot
          if(bc_u_s(bcv,m) < zero) then
             write(err_msg, 1300) trim(ivar('BC_U_s',bcv,m)), '>'
             call flush_err_msg(abort=.true.)
          endif
       enddo

    case('S')
       if(bc_v_g(bcv) > zero) then
          write(err_msg,1300) trim(ivar('BC_V_g',bcv)), '<'
          call flush_err_msg
       endif
       do m = 1, m_tot
          if(bc_v_s(bcv,m) > zero) then
             write(err_msg, 1300) trim(ivar('BC_V_s',bcv,m)), '<'
             call flush_err_msg(abort=.true.)
          endif
       enddo

    case('N')
       if(BC_V_G(BCV) < ZERO) then
          write(err_msg,1300) trim(iVar('BC_V_g',BCV)), '>'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_V_S(BCV,M) < ZERO) then
             write(err_msg, 1300) trim(iVar('BC_V_s',BCV,M)), '>'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('B')
       if(BC_W_G(BCV) > ZERO) then
          write(err_msg,1300) trim(iVar('BC_W_g',BCV)), '<'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_W_S(BCV,M) > ZERO) then
             write(err_msg, 1300) trim(iVar('BC_W_s',BCV,M)), '<'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    case('T')
       if(BC_W_G(BCV) < ZERO) then
          write(err_msg,1300) trim(iVar('BC_W_g',BCV)), '>'
          call flush_err_msg
       endif
       do M = 1, M_TOT
          if(BC_W_S(BCV,M) < ZERO) then
             write(err_msg, 1300) trim(iVar('BC_W_s',BCV,M)), '>'
             call flush_err_msg(ABORT=.true.)
          endif
       enddo

    end select

    call finl_err_msg
1300 format('Error 1300: Invalid flow direction.',/&
        A,' should be ', A,' zero. ',/&
        'Please correct the input deck.')

  end subroutine check_bc_vel_inflow



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: check_bc_vel_outflow                                     !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
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
