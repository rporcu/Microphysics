module check_collision_model

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding ,    only: c_int

  use error_manager,  only: init_err_msg, flush_err_msg, finl_err_msg,  &
                            err_msg, ival, ivar

  use param, only: zero, half, one, undefined
  use param, only: is_undefined, is_defined

  implicit none

  private
  public :: check_collision_model_lsd, check_collision_model_hertz

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_COLLISION_MODEL_LSD                               !
!                                                                      !
!  Purpose: Check user input data for DES collision calculations.      !
!                                                                      !
!  References:                                                         !
!   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
!   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
!      65-149 (pages 94-95)                                            !
!   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_collision_model_lsd

    use constant,       only: MMAX
    use discretelement, only: kn, kn_w, kt_fac, kt_w_fac, &
         & des_en_input, des_en_wall_input, &
         & des_etat_fac, des_etat_w_fac

    integer :: m, l, lc

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DEM_COLL_LSD")

    ! Check for particle-particle normal spring constants.
    if(is_undefined(kn)) then
       write(err_msg, 1000) 'KN'
       call flush_err_msg(abort=.true.)
    endif

    ! Check for particle-wall normal spring constants.
    if(is_undefined(kn_w)) then
       write(err_msg, 1000) 'KN_W'
       call flush_err_msg(abort=.true.)
    endif

    ! Check for particle-particle tangential spring constant factors.
    if(kt_fac > one .or. kt_fac < zero) then
       write(err_msg,1001) 'KT_FAC', trim(ival(kt_fac))
       call flush_err_msg(abort=.true.)
    endif

    ! Check for particle-wall tangential spring constant factors.
    if(kt_w_fac > one .or. kt_w_fac < zero) then
       write(err_msg,1001) 'KT_W_FAC', trim(ival(kt_w_fac))
       call flush_err_msg(abort=.true.)
    endif

    ! Check for particle-particle tangential damping coefficients
    if(des_etat_fac > one .or. des_etat_fac < zero) then
       write(err_msg,1001) 'DES_ETAT_FAC', ival(des_etat_fac)
       call flush_err_msg(abort=.true.)
    endif

    ! Check for particle-wall tangential damping coefficients
    if(des_etat_w_fac > one .or. des_etat_w_fac < zero) then
       write(err_msg,1001) 'DES_ETAT_W_FAC', ival(des_etat_w_fac)
       call flush_err_msg(abort=.true.)
    endif

    ! Check particle-particle normal restitution coefficient
    lc = 0
    do m = 1, mmax
       do l = m, mmax
          lc = lc+1
          if(is_undefined(des_en_input(lc))) then
             write(err_msg,1000) trim(ivar('DES_EN_INPUT',lc))
             call flush_err_msg(abort=.true.)
          elseif(des_en_input(lc) > one .or.                         &
            des_en_input(lc) < zero) then
             write(err_msg,1001) trim(ivar('DES_EN_INPUT',lc)),      &
               trim(ival(des_en_input(lc)))
             call flush_err_msg(abort=.true.)
          endif
       enddo

       if(is_undefined(des_en_wall_input(m))) then
          write(err_msg,1000) trim(ivar('DES_EN_WALL_INPUT',m))
          call flush_err_msg(abort=.true.)
       elseif(des_en_wall_input(m) > one .or.                        &
            des_en_wall_input(m) < zero) then
          write(err_msg,1001) trim(ivar('DES_EN_WALL_INPUT',m)),     &
               trim(ival(des_en_wall_input(m)))
          call flush_err_msg(abort=.true.)
       endif
    enddo

    call finl_err_msg

1000 format('Error 1000: Required input not specified: ',A,/&
       'Please correct the input deck.')

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
       'Please correct the input deck.')

  end subroutine check_collision_model_lsd

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_COLLISION_MODEL_HERTZ                             !
!                                                                      !
!  Purpose: Check user input data for Hertzian collisions.             !
!                                                                      !
!  References:                                                         !
!   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
!   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
!      65-149 (pages 94-95)                                            !
!   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine check_collision_model_hertz

  use constant,       only: mmax
  use discretelement, only: des_en_input, des_en_wall_input,  &
    des_et_input, des_et_wall_input, e_young, ew_young, v_poisson, vw_poisson

    integer           :: m, l, lc
    character(len=64) :: msg

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DEM_COLL_HERTZ")

    lc=0
    do m=1,mmax
       if(is_undefined(e_young(m))) then
          MSG=''; write(MSG,"('Phase ',I2,' Youngs modulus')") M
          write(err_msg,1002) 'E_YOUNG', msg
          call flush_err_msg(abort=.true.)
       endif
       if(is_undefined(v_poisson(m))) then
          MSG=''; write(MSG,"('Phase ',I2,' Poissons ratio')") M
          write(err_msg,1002) 'V_POISSON', msg
          call flush_err_msg(abort=.true.)
       elseif(v_poisson(m) > 0.5d0 .or.                              &
            v_poisson(m) <= -one) then
          write(err_msg,1001) trim(ivar('V_POISSON',m)),             &
               ival(v_poisson(m))
          call flush_err_msg(abort=.true.)
       endif
    enddo

    ! check young's modulus and poisson ratio
    if(is_undefined(ew_young)) then
       MSG='Wall value for Youngs modulus'
       write(err_msg,1002) 'EW_YOUNG', msg
       call flush_err_msg(abort=.true.)
    endif

    if(is_undefined(vw_poisson)) then
       MSG='Wall value for Poissons ratio'
       write(err_msg,1002) 'VW_POISSON', msg
       call flush_err_msg(abort=.true.)
    elseif (vw_poisson > 0.5d0 .or. vw_poisson <= -one) then
       write(err_msg,1001) 'VW_POISSON',ival(vw_poisson)
       call flush_err_msg(abort=.true.)
    endif

    do m=1,mmax
       do l=m,mmax
          lc = lc+1

          ! Check particle-particle normal restitution coefficient
          if(is_undefined(des_en_input(lc))) then
             write(err_msg,1000) trim(ivar('DES_EN_INPUT',lc))
             call flush_err_msg(abort=.true.)

          elseif(des_en_input(lc) > one .or.                         &
               des_en_input(lc) < zero) then
             write(err_msg,1001) trim(ivar('DES_EN_INPUT',lc)),      &
                  trim(ival(des_en_input(lc)))
             call flush_err_msg(abort=.true.)
          endif

          ! Check particle-particle tangential restitution coefficient
          if(is_undefined(des_et_input(m))) then
             write(err_msg,1000) trim(ivar('DES_ET_INPUT',m))
             call flush_err_msg(abort=.true.)
          elseif(des_et_input(m) > one .or.                          &
               des_et_input(m) < zero) then
             write(err_msg,1001) trim(ivar('DES_ET_INPUT',m)),       &
                  ival(des_et_input(m))
             call flush_err_msg(abort=.true.)
          endif
       enddo

       ! Check particle-wall normal restitution coefficient.
       if(is_undefined(des_en_wall_input(m))) then
          write(err_msg,1000) trim(ivar('DES_EN_WALL_INPUT',m))
          call flush_err_msg(abort=.true.)
       elseif(des_en_wall_input(m) > one .or.                        &
            des_en_wall_input(m) < zero) then
          write(err_msg,1001) trim(ivar('DES_EN_WALL_INPUT',m)),     &
               trim(ival(des_en_wall_input(m)))
          call flush_err_msg(abort=.true.)
       endif

       ! Check particle-wall tangential restitution coefficient
       if(is_undefined(des_et_wall_input(m))) then
          write(err_msg,1000) trim(ivar('DES_ET_WALL_INPUT',m))
          call flush_err_msg(abort=.true.)
       elseif(des_et_wall_input(m) > one .or.                        &
            des_et_wall_input(m) < zero) then
          write(err_msg,1001) trim(ivar('DES_ET_WALL_INPUT',m)),     &
               trim(ival(des_et_wall_input(m)))
          call flush_err_msg(abort=.true.)
       endif
    enddo

    call finl_err_msg

1000 format('Error 1000: Required input not specified: ',A,/&
       'Please correct the input deck.')

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
          'Please correct the input deck.')

1002 format('Error 1002: Required input not specified: ',A,/          &
         'Description:',A,/'Please correct the input deck.')

  end subroutine check_collision_model_hertz

end module check_collision_model
