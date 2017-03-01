module check_solids_common_discrete_module

  use bl_fort_module, only: c_real
  use iso_c_binding , only: c_int
  use run,            only: IFILE_NAME
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg,    &
       &  ivar,  ival, err_msg

  implicit none
  private 

  public check_solids_common_discrete

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_SOLIDS_COMMON_DISCRETE                            !
  !  Author: J.Musser                                   Date: 02-FEB-14  !
  !                                                                      !
  !  Purpose:                                                            !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_solids_common_discrete

    use constant,       only: mmax,  d_p0
    use param1,         only: undefined
    use discretelement, only: des_intg_method, intg_adams_bashforth, &
                            & intg_euler, max_radius, min_radius,    &
                            & do_old

    integer :: lM  ! Solids phase Index

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_COMMON_DISCRETE")

    max_radius = -undefined
    min_radius =  undefined

    ! determine the maximum particle size in the system (max_radius), which
    ! in turn is used for various tasks
    do lm=1, mmax
       max_radius = max(max_radius, 0.5d0*d_p0(lm))
       min_radius = min(min_radius, 0.5d0*d_p0(lm))
    enddo

    ! Check for valid integration method
    select case(trim(DES_INTG_METHOD))
    case ('EULER')
       intg_euler = .true.
       intg_adams_bashforth = .false.
       !DES_INTG_METHOD_ENUM = 1
    case ('ADAMS_BASHFORTH')
       intg_euler = .false.
       intg_adams_bashforth = .true.
       !DES_INTG_METHOD_ENUM = 2
    case DEFAULT
       write(err_msg,2020) trim(des_intg_method), trim(IFILE_NAME)
       call flush_err_msg(abort=.true.)
       
2020   format('Error 2020:Invalid DES_INGT_METHOD: ',A,/'Please ',      &
            'correct the ',A,' file.')
       
    end select
    
    do_old = intg_adams_bashforth
    
    ! Check geometry constrains.
    call check_solids_common_discrete_geometry
    
    call finl_err_msg
    
  end subroutine check_solids_common_discrete
  
  
  
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY                   !
  !  Author: J.Musser                                   Date: 11-DEC-13  !
  !                                                                      !
  !  Purpose: Check user input data                                      !
  !                                                                      !
  !  Comments: Geometry checks were moved here from CHECK_DES_DATA.      !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_solids_common_discrete_geometry

    use geometry,       only: zlength
    use discretelement, only: des_continuum_coupled, max_radius

    real(c_real) :: min_depth
    
    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY")
    
    
    if(des_continuum_coupled)then
       ! check that the depth of the simulation exceeds the largest particle
       ! to ensure correct calculation of volume fraction. this is important
       ! for coupled simulations.
       min_depth = 2.0d0*max_radius
       if(zlength < min_depth)then
          write(err_msg, 1300) trim(IFILE_NAME)
          call flush_err_msg(abort=.false.)
       endif
    endif
    
1300 format('Error 1300: The maximum particle diameter exceeds the ', &
          'simulation',/'depth (ZLENGTH). Please correct the ',A,' ',&
          'file.')
    
    call finl_err_msg
    
  end subroutine check_solids_common_discrete_geometry
  
end module check_solids_common_discrete_module
