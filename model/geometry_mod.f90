!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: geometry                                               C
!  Purpose: Common block containing geometry and discretization data   C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE geometry

      use bl_fort_module, only : c_real
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

! Number of cells +/-
      integer :: imax1, jmax1, kmax1
      integer :: imin1, jmin1, kmin1
      integer :: imax2, jmax2, kmax2
      integer :: imin2, jmin2, kmin2

      integer :: domlo(3),domhi(3)

!  one or more periodic boundary condition is used
      LOGICAL :: CYCLIC
! Variable to flag periodic boundary condition in X
      LOGICAL :: CYCLIC_X
! Variable to flag periodic boundary condition in Y
      LOGICAL :: CYCLIC_Y
! Variable to flag periodic boundary condition in Z
      LOGICAL :: CYCLIC_Z
! Variable to flag periodic bc with pressure drop in X
      LOGICAL :: CYCLIC_X_PD
! Variable to flag periodic bc with pressure drop in Y
      LOGICAL :: CYCLIC_Y_PD
! Variable to flag periodic bc with pressure drop in Z
      LOGICAL :: CYCLIC_Z_PD
! Variable to flag periodic bc with mass flux in X
      LOGICAL :: CYCLIC_X_MF
! Variable to flag periodic bc with mass flux in Y
      LOGICAL :: CYCLIC_Y_MF
! Variable to flag periodic bc with mass flux in Z
      LOGICAL :: CYCLIC_Z_MF

! Cell flags.
      integer, DIMENSION(:,:,:,:), ALLOCATABLE :: flag_mod

      END MODULE geometry
