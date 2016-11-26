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

      Use param, only: DIM_I, DIM_J, DIM_K


! Coordinates: CARTESIAN
      CHARACTER(LEN=16)     COORDINATES

! Indicates whether x or r direction is not considered
      LOGICAL :: NO_I
! Indicates whether x or r direction is considered
      LOGICAL :: DO_I
! Indicates whether y direction is not considered
      LOGICAL :: NO_J
! Indicates whether y direction is considered
      LOGICAL :: DO_J
! Indicates whether z or theta direction is not considered
      LOGICAL :: NO_K
! Indicates whether z or theta direction is considered
      LOGICAL :: DO_K

! Reactor length in the x or r direction
      DOUBLE PRECISION :: XLENGTH
! Reactor length in the y direction
      DOUBLE PRECISION :: YLENGTH
! Reactor length in the z or theta direction
      DOUBLE PRECISION :: ZLENGTH

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

! Cell sizes in the x or r direction
      DOUBLE PRECISION :: DX, DY, DZ
      DOUBLE PRECISION :: oDX, oDY, oDZ

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

! Cell flags.
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: FLAG
! Flag for the East surface
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: FLAG_E
! Flag for North surface
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: FLAG_N
! Flag for Top surface
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: FLAG_T
! Cell flags (bc/ic conditions)
! Allocatable type causes PG internal error, Ed's soln: pointers
!      CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: ICBC_FLAG
      character(LEN=3),  dimension(:), pointer :: icbc_flag


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


! East face area - scalar cell
      DOUBLE PRECISION :: AYZ
! North face area - scalar cell
      DOUBLE PRECISION :: AXZ
! Top face area - scalar cell
      DOUBLE PRECISION :: AXY
! Cell volume - scalar cell
      DOUBLE PRECISION :: VOL

! Total volume of cell's DES stencil neighbors
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: vol_surr

      END MODULE geometry
