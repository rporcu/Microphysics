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
      use check_solids_model_prereqs_module, only: check_solids_model_prereqs
      use check_solids_phases_module, only: check_solids_phases
      use error_manager  , only: init_error_manager
      use init_namelist_module, only: init_namelist
      use read_namelist_module, only: read_namelist
      use set_parameters_module, only: set_parameters
      use param1, only: is_undefined

      use run, only: detect_stall

      implicit none

      real(c_real), intent(  out) :: dt

      ! This module call routines to initialize the namelist variables.
      call init_namelist

      ! Read in the namelist variables from the ascii input file.
      call read_namelist(dt)

      ! Initialize the error manager. This call occurs after the mfix.dat
      ! is read so that message verbosity can be set and the .LOG file
      ! can be opened.
      call init_error_manager

      ! Check the minimum solids phase requirements.
      call check_solids_model_prereqs

      call check_run_control(dt)

      call check_gas_phase
      call check_solids_phases
      call set_parameters


      if(is_undefined(dt)) detect_stall = .false.

      end subroutine get_data

end module get_data_module
