!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: geometry                                               C
!  Purpose: Common block containing geometry and discretization data   C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      module geometry

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use param, only: DIM_I, DIM_J, DIM_K

! Coordinates: CARTESIAN
      CHARACTER(LEN=16)     COORDINATES

! Reactor length in the x or r direction
      real(c_real) :: XLENGTH
! Reactor length in the y direction
      real(c_real) :: YLENGTH
! Reactor length in the z or theta direction
      real(c_real) :: ZLENGTH

! Number of cells
      integer :: imax, jmax, kmax

      integer :: domlo(3),domhi(3)

      !  one or more periodic boundary condition is used
      logical :: cyclic

      ! Variable to flag periodic boundary condition in X
      logical :: cyclic_x

      ! Variable to flag periodic boundary condition in Y
      logical :: cyclic_y

      ! Variable to flag periodic boundary condition in Z
      logical :: cyclic_z

      ! Variable to flag periodic bc with pressure drop in X
      logical :: cyclic_x_pd

      ! Variable to flag periodic bc with pressure drop in Y
      logical :: cyclic_y_pd

      ! Variable to flag periodic bc with pressure drop in Z
      logical :: cyclic_z_pd

      ! Variable to flag periodic bc with mass flux in X
      logical :: cyclic_x_mf

      ! Variable to flag periodic bc with mass flux in Y
      logical :: cyclic_y_mf

      ! Variable to flag periodic bc with mass flux in Z
      logical :: cyclic_z_mf

      end module geometry
