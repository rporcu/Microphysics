module check_point_sources_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use param,         only: undefined, undefined_i, is_defined, is_undefined, zero
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

    use ps,    only: ps_defined
    use param, only: dim_ps

    integer                  :: PSV
    real(c_real), intent(in) :: dx,dy,dz

    ! Initialize the error manager.
    call init_err_msg("CHECK_POINT_SOURCES")

    ! Determine which PSs are DEFINED
    call check_ps_geometry(dx,dy,dz)

    ! Loop over all PS arrays.
    do psv = 1, dim_ps

       ! Verify user input for defined defined PS.
       if(ps_defined(psv)) then
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
  subroutine check_ps_geometry(dx, dy, dz)


    use param, only: dim_ps
    use ps, only: ps_defined
    use ps, only: ps_x_e, ps_y_n, ps_z_t
    use ps, only: ps_x_w, ps_y_s, ps_z_b

    use calc_cell_module, only: calc_cell_ps

    real(c_real), intent(in) :: dx,dy,dz

    integer :: psv, ier
    integer :: i_w , i_e , j_s , j_n , k_b , k_t

    ! Initialize the error manager.
    call init_err_msg("CHECK_PS_GEOMETRY")

    ! Determine which point source indices have values.
    psv_lp: do psv = 1, dim_ps

       if(ps_defined(psv)) then
          if(is_undefined(ps_x_w(psv))) then
             write(err_msg,1101) psv, 'PS_X_w'
             call flush_err_msg(abort=.true.)
          endif
          if(is_undefined(ps_x_e(psv))) then
             write(err_msg,1101) psv, 'PS_X_e '
             call flush_err_msg(abort=.true.)
          endif
          if(is_undefined(ps_y_s(psv))) then
             write(err_msg,1101) psv, 'PS_Y_s'
             call flush_err_msg(abort=.true.)
          endif
          if(is_undefined(ps_y_n(psv))) then
             write(err_msg,1101) psv, 'PS_Y_n'
             call flush_err_msg(abort=.true.)
          endif
          if(is_undefined(ps_z_b(psv))) then
             write(err_msg,1101) psv, 'PS_Z_b'
             call flush_err_msg(abort=.true.)
          endif
          if(is_undefined(ps_z_t(psv))) then
             write(err_msg,1101) psv, 'PS_Z_t'
             call flush_err_msg(abort=.true.)
          endif

1101   format('Error 1101: Point source ',I3,' is ill-defined.',/A,     &
          ' are not specified.',/'Please correct the input deck.')

           call calc_cell_ps(psv, dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)

           ier = 0

           if(k_b > k_t) ier = 1
           if(j_s > j_n) ier = 1
           if(i_w > i_e) ier = 1

           if(ier /= 0)then
              write(err_msg,1103) psv,                         &
                 'X', ps_x_w(psv), ps_x_e(psv), 'I', i_w, i_e, &
                 'Y', ps_y_s(psv), ps_y_n(psv), 'J', j_s, j_n, &
                 'Z', ps_z_b(psv), ps_z_t(psv), 'K', k_b, k_t
              call flush_err_msg(abort=.true.)
           endif

 1103 FORMAT('Error 1101: Invalid location specified for PS ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the input deck.')

        endif
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

      use ps, only: ps_massflow_g, ps_u_g, ps_v_g, ps_w_g
      use param, only: zero, equal

      integer, intent(in) :: PSV

    ! Initialze the error manager.
      CALL init_err_msg("CHECK_PS_GAS_PHASE")

    ! Check mass flow and velocity
      if(equal(ps_massflow_g(psv),zero)) then
         if(.not.equal(ps_u_g(psv),zero) .or. &
            .not.equal(ps_v_g(psv),zero) .or. &
            .not.equal(ps_w_g(psv),zero)) then

            write(err_msg,1100) psv, trim(ivar('PS_MASSFLOW_G',Psv))
            call flush_err_msg(abort=.true.)

1100  format('Error 1100: Invalid specification for point source ',I3,&
         '.',/A,' is zero but velocity is given.',/&
         'Please correct the input deck.')

         endif

      elseif(ps_massflow_g(psv) < zero) then
       write(err_msg,1102) psv, trim(ivar('PS_MASSFLOW_g',psv))
       call flush_err_msg(abort=.true.)

1102  format('Error 1102: Invalid specifications for point source ',I3,&
         '.',/A,' < 0.0. Point sources can only add mass to a system',/&
         'Please correct the input deck.')

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

    if(ps_massflow_g(psv) > zero) then
       write(err_msg,1010) trim(ivar('PS_MASSFLOW_G',psv))
       call flush_err_msg(abort=.true.)
    elseif(ps_u_g(psv) > zero) then
       write(err_msg,1010) trim(ivar('PS_U_g',psv))
       call flush_err_msg(abort=.true.)
    elseif(ps_v_g(psv) > zero) then
       write(err_msg,1010) trim(ivar('PS_V_g',psv))
       call flush_err_msg(abort=.true.)
    elseif(ps_w_g(psv) > zero) then
       write(err_msg,1010) trim(ivar('PS_W_g',psv))
       call flush_err_msg(abort=.true.)
    endif

    call finl_err_msg

1010 format('Error 1010: ',A,' specified in an undefined PS region.')

  end subroutine check_ps_overflow

end module check_point_sources_module
