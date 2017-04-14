module get_data_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  subroutine: get_data                                                !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine get_data(dt)

    use check_gas_phase_module, only: check_gas_phase
    use check_run_control_module, only: check_run_control
    use check_solids_phases_module, only: check_solids_phases
    use error_manager  , only: init_error_manager
    use init_namelist_module, only: init_namelist
    use read_namelist_module, only: read_namelist
    use param, only: is_undefined

    ! Cyclic domain flags.
    use bc, only: cyclic_x, cyclic_x_pd, cyclic_x_mf
    use bc, only: cyclic_y, cyclic_y_pd, cyclic_y_mf
    use bc, only: cyclic_z, cyclic_z_pd, cyclic_z_mf

    ! Flag for specificed constant mass flux.
    use bc, only: flux_g

    use run, only: detect_stall, dem_solids
    use discretelement, only: particle_types
    use discretelement, only: des_continuum_coupled, des_oneway_coupled
    use fld_const, only: ro_g0

    use param, only: is_defined
    use constant, only: mmax

    implicit none

    real(c_real), intent(  out) :: dt

    ! This module call routines to initialize the namelist variables.
    call init_namelist

    ! Read in the namelist variables from the ascii input file.
    call read_namelist(dt)

    ! Disable detect stall for steady state
    if(is_undefined(dt)) detect_stall = .false.

    ! Determine the cyclic direction with a specified mass flux
    cyclic_x_mf = (is_defined(flux_g) .and. cyclic_x_pd)
    cyclic_y_mf = (is_defined(flux_g) .and. cyclic_y_pd)
    cyclic_z_mf = (is_defined(flux_g) .and. cyclic_z_pd)

    ! Force the cyclic flag if cyclic with pressure drop.
    if (cyclic_x_pd) cyclic_x = .true.
    if (cyclic_y_pd) cyclic_y = .true.
    if (cyclic_z_pd) cyclic_z = .true.

    ! Set flag for coupled simulations
    des_continuum_coupled = (particle_types>0) .and. (abs(ro_g0) > 0.0d0)
    dem_solids = (particle_types > 0)

    mmax = particle_types
    write(*,*) 'particle types',particle_types
    ! Overwrite user settings if no Lagrangian solids
    if(particle_types==0) then
       des_continuum_coupled = .false.
       des_oneway_coupled = .false.
    endif

    ! Initialize the error manager. This call occurs after the mfix.dat
    ! is read so that message verbosity can be set and the .LOG file
    ! can be opened.
    call init_error_manager

    call check_run_control(dt)

    call check_gas_phase
    call check_solids_phases

  end subroutine get_data

end module get_data_module
