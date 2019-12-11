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
  subroutine get_data (mfix_dat)

    use init_namelist_module, only: init_namelist
    use read_namelist_module, only: read_namelist

!    use fld_const, only: ro_g0 ! UNUSED_VARIABLE

    use param, only: is_undefined, is_defined

    implicit none

    character(len=*) :: mfix_dat

    ! This module call routines to initialize the namelist variables.
    call init_namelist

    ! Read in the namelist variables from the ascii input file.
    call read_namelist(mfix_dat)

  end subroutine get_data

end module get_data_module
