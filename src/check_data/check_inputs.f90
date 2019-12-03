!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: check_inputs                                            !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine check_inputs() bind(C, name="check_inputs")

  use amrex_fort_module, only : rt => amrex_real

  use error_manager, only: init_error_manager

  use check_gas_prop_module, only: check_gas_properties
  use check_particle_prop_module, only: check_particle_properties

  ! Initialize the error manager. This call occurs after the mfix.dat
  ! is read so that message verbosity can be set and the .LOG file
  ! can be opened.
  call init_error_manager



  call check_gas_properties
  call check_particle_properties

end subroutine check_inputs
