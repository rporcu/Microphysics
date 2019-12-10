module check_collision_model

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding ,    only: c_int

  use error_manager,  only: init_err_msg, flush_err_msg, finl_err_msg,  &
                            err_msg, ival, ivar

  use param, only: zero, half, one, undefined
  use param, only: is_undefined, is_defined

  implicit none

  private
  public :: check_collision_model_lsd

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
    use discretelement, only:  kt_fac, kt_w_fac, &
         & des_en_input, des_en_wall_input, &
         & des_etat_fac, des_etat_w_fac

    integer :: m, l, lc

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DEM_COLL_LSD")



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


end module check_collision_model
