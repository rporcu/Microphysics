module get_data_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  subroutine: get_data                                                C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine get_data(time, dt)

      use check_gas_phase_module, only: check_gas_phase
      use check_numerics_module, only: check_numerics
      use check_output_control_module, only: check_output_control
      use check_run_control_module, only: check_run_control
      use check_solids_model_prereqs_module, only: check_solids_model_prereqs
      use check_solids_phases_module, only: check_solids_phases
      use error_manager  , only: init_error_manager
      use init_namelist_module, only: init_namelist
      use read_namelist_module, only: read_namelist
      use set_parameters_module, only: set_parameters
      use write_header_module, only: write_header

      implicit none

      real(c_real), intent(  out) :: time, dt

      ! This module call routines to initialize the namelist variables.
      call init_namelist

      ! Read in the namelist variables from the ascii input file.
      call read_namelist(time, dt)

      ! Initialize the error manager. This call occurs after the mfix.dat
      ! is read so that message verbosity can be set and the .LOG file
      ! can be opened.
      call init_error_manager

      ! Write header in the .LOG file and to screen.
      ! Not sure if the values are correct or useful
      call write_header

      ! Check the minimum solids phase requirements.
      call check_solids_model_prereqs

      call check_run_control(time, dt)
      call check_numerics
      call check_output_control

      call check_gas_phase
      call check_solids_phases
      call set_parameters

      end subroutine get_data

end module get_data_module
