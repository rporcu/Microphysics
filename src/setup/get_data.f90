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
  subroutine get_data(mfix_dat)

    use init_namelist_module, only: init_namelist
    use read_namelist_module, only: read_namelist

    use run, only: dem_solids

    use constant, only: mmax
    use discretelement, only: particle_types

    use fld_const, only: ro_g0

    use param, only: is_undefined, is_defined

    implicit none

    character(len=*) :: mfix_dat

    ! This module call routines to initialize the namelist variables.
    call init_namelist

    ! Read in the namelist variables from the ascii input file.
    call read_namelist(mfix_dat)

    dem_solids = (particle_types > 0)

    mmax = particle_types

  end subroutine get_data

end module get_data_module
