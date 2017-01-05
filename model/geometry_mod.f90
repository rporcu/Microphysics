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

! Starting index in the x or r direction
      INTEGER :: IMIN1
! Starting index in the y direction
      INTEGER :: JMIN1
! Starting index in the z or theta direction
      INTEGER :: KMIN1

! Number of cells in the x or r direction
      INTEGER :: IMAX
! Number of cells in the y direction
      INTEGER :: JMAX
! Number of cells in the z or theta direction
      INTEGER :: KMAX

! Number of cells in the x or r direction + 1
      INTEGER :: IMAX1
! Number of cells in the y direction + 1
      INTEGER :: JMAX1
! Number of cells in the z or theta direction + 1
      INTEGER :: KMAX1

! Number of cells in the x or r direction + 2
      INTEGER :: IMAX2
! Number of cells in the y direction + 2
      INTEGER :: JMAX2
! Number of cells in the z or theta direction + 2
      INTEGER :: KMAX2

! IMAX2 * JMAX2
      INTEGER :: IJMAX2
! IMAX2 * JMAX2 * KMAX2
      INTEGER :: IJKMAX2
! IMAX2 * JMAX2 * KMAX2
      INTEGER :: IJKMAX3
! IJMAX2 + 1
      INTEGER :: IJKMIN1
! IJKMAX2 - IJMAX2
      INTEGER :: IJKMAX1

! For discretization in parallel
      INTEGER :: IMIN2, JMIN2, KMIN2
      INTEGER :: IMIN3, JMIN3, KMIN3
      INTEGER :: IMAX3, JMAX3, KMAX3

! For 4th order discretization in parallel
      INTEGER :: IMIN4, JMIN4, KMIN4
      INTEGER :: IMAX4, JMAX4, KMAX4
      INTEGER :: IJKMAX4, IJKMIN4

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
