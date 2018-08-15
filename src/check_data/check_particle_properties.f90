module check_particle_prop_module

  use error_manager, only: init_err_msg, flush_err_msg, finl_err_msg, &
                           err_msg, ivar, ival

  private
  public :: check_particle_properties

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  SUBROUTINE: CHECK_SOLIDS_PHASES                                     !
!                                                                      !
!  Purpose: Driver routine for calls to solids phase checks.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_particle_properties

! Global Variables:
!---------------------------------------------------------------------//
! Runtime flag specifying DEM solids
    use discretelement, only: mew, mew_w
    use drag, only: drag_type, drag_type_enum, invalid_drag
    use discretelement, only: particle_types
    use discretelement, only: des_coll_model_enum, invalid_coll
    use discretelement, only: des_coll_model, lsd, hertzian

    use param, only: dim_m
    use param, only: zero, one

! Global Module procedures:
!---------------------------------------------------------------------//
    use check_collision_model, only: check_collision_model_lsd
    use check_collision_model, only: check_collision_model_hertz

    use param, only: is_undefined, is_defined

    implicit none

!......................................................................!

! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_PHASES")

    if(particle_types > 0) then

       if(drag_type_enum == invalid_drag) then
          write(err_msg, 1001) 'DRAG_TYPE', trim(adjustl(drag_type))
          call flush_err_msg(abort=.true.)
       endif

       if(is_undefined(mew)) then
          write(err_msg,1000) 'MEW'
          call flush_err_msg(abort=.true.)
       elseif (mew < zero .or. mew_w > one) then
          write(err_msg,1001) 'MEW', trim(ival(mew))
          call flush_err_msg(abort=.true.)
       endif

       if(is_undefined(mew_w)) then
          write(err_msg,1000) 'MEW_W'
          call flush_err_msg(abort=.true.)
       elseif(mew_w < zero .or. mew_w > one) then
          write(err_msg,1001) 'MEW_W', trim(ival(mew_w))
          call flush_err_msg(abort=.true.)
       endif

       if(des_coll_model_enum == invalid_coll) then
          write(err_msg,1001) 'DES_COLL_MODEL', &
            trim(adjustl(des_coll_model))
          call flush_err_msg(abort=.true.)
       else if(des_coll_model_enum == lsd)then
          call check_collision_model_lsd
       else if(des_coll_model_enum == hertzian) then
          call check_collision_model_hertz
       endif
   endif

   call finl_err_msg

   return

1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the input deck.')

1001 format('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the input deck.')

    end subroutine check_particle_properties

  end module check_particle_prop_module
