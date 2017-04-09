module check_geometry_prereqs_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use run,            only: IFILE_NAME
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg,    &
                          &  ivar,  ival, err_msg

  implicit none
  private 

  public check_geometry_prereqs

contains
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_GEOMETRY_PREREQS                                  !
  !  Purpose: Check the distributed parallel namelist variables.         !
  !                                                                      !
  !  Author: P. Nicoletti                               Date: 14-DEC-99  !
  !  Reviewer: J.Musser                                 Date: 16-Jan-14  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_geometry_prereqs

    use geometry, only: coordinates

    ! Initialize the error manager.
    call init_err_msg("CHECK_GEOMETRY_PREREQS")
    
    select case(trim(coordinates))
       
    case ('CARTESIAN')
       
    case default
       write(err_msg, 1103) trim(IFILE_NAME)
       call flush_err_msg(abort=.true.)
       
1103   format('Error 1103: Unknown COORDINATES specified. Please ',     &
            'correct the ',/A,' file.')
       
    end select
    
    call finl_err_msg

  end subroutine check_geometry_prereqs

end module check_geometry_prereqs_module
