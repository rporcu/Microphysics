module get_data_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  subroutine: get_data                                                !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine get_data(mfix_dat, dt)

    use init_namelist_module, only: init_namelist
    use read_namelist_module, only: read_namelist

    use run, only: dem_solids

    use constant, only: mmax
    use discretelement, only: particle_types

    use discretelement, only: des_continuum_coupled, des_oneway_coupled
    use fld_const, only: ro_g0

    use drag, only: drag_type, drag_type_enum, user_drag, wen_yu, &
                    gidaspow, bvk2, wen_yu_pcf, gidaspow_pcf  

    use param, only: is_undefined, is_defined

    implicit none

    character(len=*) :: mfix_dat
    real(rt), intent(  out) :: dt

    ! This module call routines to initialize the namelist variables.
    call init_namelist

    ! Read in the namelist variables from the ascii input file.
    call read_namelist(mfix_dat, dt)

    ! Set flag for coupled simulations
    des_continuum_coupled = (particle_types>0) .and. (abs(ro_g0) > 0.0d0)
    dem_solids = (particle_types > 0)

    mmax = particle_types
    ! Overwrite user settings if no Lagrangian solids
    if(particle_types==0) then
       des_continuum_coupled = .false.
       des_oneway_coupled = .false.
    endif

    select case(trim(adjustl(drag_type)))
    case ('USER_DRAG','USR_DRAG'); drag_type_enum = user_drag
    case ('WEN_YU'); drag_type_enum = wen_yu
    case ('GIDASPOW'); drag_type_enum = gidaspow
    case ('BVK2'); drag_type_enum = bvk2
    case ('WEN_YU_PCF'); drag_type_enum = wen_yu_pcf
    case ('GIDASPOW_PCF'); drag_type_enum = gidaspow_pcf
    end select

  end subroutine get_data

end module get_data_module
