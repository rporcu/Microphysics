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

    use discretelement, only: des_coll_model, des_coll_model_enum, &
                              lsd, hertzian

    use drag, only: drag_type, drag_type_enum, syam_obrien, gidaspow, &
                    gidaspow_pcf, gidaspow_blend, gidaspow_blend_pcf, &
                    wen_yu, wen_yu_pcf, koch_hill, koch_hill_pcf, bvk,&
                    user_drag

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
    case ('SYAM_OBRIEN'); drag_type_enum = syam_obrien
    case ('GIDASPOW'); drag_type_enum = gidaspow
    case ('GIDASPOW_BLEND'); drag_type_enum = gidaspow_blend
    case ('WEN_YU'); drag_type_enum = wen_yu
    case ('KOCH_HILL'); drag_type_enum = koch_hill
    case ('BVK'); drag_type_enum = bvk
    case ('GIDASPOW_PCF'); drag_type_enum = gidaspow_pcf
    case ('GIDASPOW_BLEND_PCF'); drag_type_enum = gidaspow_blend_pcf
    case ('WEN_YU_PCF'); drag_type_enum = wen_yu_pcf
    case ('KOCH_HILL_PCF'); drag_type_enum = koch_hill_pcf
    case ('USER_DRAG','USR_DRAG'); drag_type_enum = user_drag
    end select

    ! Check collision model specific parameters.
    select case (trim(adjustl(des_coll_model)))
    case('LSD'); des_coll_model_enum = lsd
    case('HERTZIAN'); des_coll_model_enum = hertzian
    end select

  end subroutine get_data

end module get_data_module
