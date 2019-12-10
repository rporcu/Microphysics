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

    use param, only: dim_m
    use param, only: zero, one

! Global Module procedures:
!---------------------------------------------------------------------//
    use check_collision_model, only: check_collision_model_lsd

    use param, only: is_undefined, is_defined

    implicit none

!......................................................................!

! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_PHASES")


       ! call check_collision_model_lsd

   call finl_err_msg

   return

1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the input deck.')

1001 format('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the input deck.')

    end subroutine check_particle_properties

  end module check_particle_prop_module
