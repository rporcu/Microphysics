module check_solids_model_prereqs_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg, &
                          & ivar, ival, err_msg

  implicit none
  private

  public check_solids_model_prereqs

CONTAINS
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  SUBROUTINE: CHECK_SOLIDS_MODEL_PREREQS                              !
  !  Purpose: Check the distributed parallel namelist variables.         !
  !                                                                      !
  !  Author: P. Nicoletti                               Date: 14-DEC-99  !
  !  Reviewer: J.Musser                                 Date: 16-Jan-14  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  SUBROUTINE CHECK_SOLIDS_MODEL_PREREQS

    use constant,       only: mmax
    use fld_const,      only: ro_g0
    use param,          only: dim_m
    use run,            only: tfm_solids, tfm_count, dem_solids, &
                            & dem_count, pic_solids, pic_count,  &
                            & solids_model
    use discretelement, only: des_continuum_coupled, des_oneway_coupled

    integer      :: M ! Phase index

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_MODEL_PREREQS")

    ! Loop over the phases to see what was specified.
    do m=1, dim_m
       solids_model(m) = trim(adjustl(solids_model(m)))
       select case(solids_model(m))
       case ('TFM'); tfm_count = tfm_count + 1
       case ('DEM'); dem_count = dem_count + 1
       case ('PIC'); pic_count = pic_count + 1
       case ('---')
       case default
          write(err_msg,1001) ivar('SOLIDS_MODEL',m), solids_model(m)
          call flush_err_msg(ABORT=.true.)

1001      format('Error 1001: Unknown solids model: ',A,' = ',A)

       end select
    enddo

    ! Clear out the unused phases.
    mmax =  tfm_count + dem_count + pic_count

    ! Set the runtime flags:
    tfm_solids = (tfm_count > 0)
    dem_solids = (dem_count > 0)
    pic_solids = (pic_count > 0)

    if(tfm_solids)then
       write(err_msg, 1002)
       call flush_err_msg(abort=.true.)
    elseif(pic_solids)then
       write(err_msg, 1003)
       call flush_err_msg(abort=.true.)
    endif

1002 format('Error 1002: TFM solids are not supported in this&
          & version of MFIX.',/'Please correct the input deck.')

1003 format('Error 1003: PIC solids are not supported in this&
          & version of MFIX.',/'Please correct the input deck.')

    ! Set flag for coupled simulations
    des_continuum_coupled = dem_solids .and. (abs(ro_g0) > 0.0d0)

    ! Overwrite user settings if no Lagrangian solids
    if(.not.dem_solids) then
       des_continuum_coupled = .false.
       des_oneway_coupled = .false.
    endif

    call finl_err_msg

  end subroutine check_solids_model_prereqs

end module check_solids_model_prereqs_module
