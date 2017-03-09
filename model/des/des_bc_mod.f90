!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_INLET                                              !
!                                                                      !
!  Purpose: Common elements needed for the des mass inflow boundary    !
!  condition.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE DES_BC

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      USE param, only: dimension_bc

      INTEGER :: DEM_BCMI
      INTEGER :: DEM_BCMO

      LOGICAL DEM_MIO  ! either inlet or outlet exists

! Map between DEM MI/MO IDs and the user input BC index.
      INTEGER :: DEM_BCMI_MAP(DIMENSION_BC)
      INTEGER :: DEM_BCMO_MAP(DIMENSION_BC)

! This array contains integers representing the mass/solid phase indices
! present at a specific boundary condtion in proportion to their
! respective number fraction at the inlet (i.e., it represents the
! particle number distribution of incoming solids at the inlet).  The
! array is scaled in size according to the parameter NUMFRAC_LIMIT.
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: DEM_BC_POLY_LAYOUT

! Particle injection time scale; used when pi_factor > 1 to keep track
! of time needed for next injection
      real(c_real), DIMENSION(:), ALLOCATABLE :: DEM_MI_TIME

! Particle injection factor; how many solid time steps (dtsolid) pass
! before the next injection of a particle. if pi_count is greater than
! 1, then pi_factor is set to 1 (i.e. multiple particles enter every
! solids time step).
      INTEGER, DIMENSION(:), ALLOCATABLE :: PI_FACTOR   !(DES_BCMI)

! Particle injection count (injection number); how many particles are
! injected in one solids time step. pi_count is set to one if
! less than 1 particle enters per solids time step.
      INTEGER, DIMENSION(:), ALLOCATABLE :: PI_COUNT   !(DES_BCMI)


! Limit on the total number of divisions (fineness) used to represent
! the particle number distribution at an inlet.
      INTEGER, PARAMETER :: NUMFRAC_LIMIT = 10000


! the dimension of this variable is equal to the number of grid
! cells in the inlet edge/face
      TYPE DEM_MI_
! Array position of next seed location.
         INTEGER :: VACANCY
! Number of positions in the layout grid.
         INTEGER :: OCCUPANTS
! Flag for polydisperse inlets.
         LOGICAL :: POLYDISPERSE
! Uniform grid dimension (width and height).
         real(c_real) :: WINDOW
! Offset for placing particles in ghost cell.
         real(c_real) :: OFFSET
! Fluid cell index associated with each grid. (I/J/K like)
         INTEGER :: L
         INTEGER, ALLOCATABLE :: W(:)
         INTEGER, ALLOCATABLE :: H(:)
! Spatial location of each grid cell's lower, bottom corder.
         real(c_real), ALLOCATABLE :: P(:)
         real(c_real), ALLOCATABLE :: Q(:)
! The rank of the owning process owning the indexed grid cell.
         INTEGER, ALLOCATABLE :: OWNER(:)
      END TYPE DEM_MI_

! Construct an array of integers in values from 1 to a calculated factor
! in a random order, which is used when placing new particles.
!      TYPE(DEM_MI_DATA), DIMENSION(:), ALLOCATABLE :: MI_ORDER

! Array linking all of the reaction data.
      TYPE(DEM_MI_), DIMENSION(:), TARGET, ALLOCATABLE :: DEM_MI

      INTEGER, ALLOCATABLE :: DEM_BCMO_IJKSTART(:)
      INTEGER, ALLOCATABLE :: DEM_BCMO_IJKEND(:)

      INTEGER, ALLOCATABLE :: DEM_BCMO_IJK(:)


      INTEGER, ALLOCATABLE :: DEM_BCMI_IJKSTART(:)
      INTEGER, ALLOCATABLE :: DEM_BCMI_IJKEND(:)

      INTEGER, ALLOCATABLE :: DEM_BCMI_IJK(:)

      CONTAINS

      END MODULE DES_BC
