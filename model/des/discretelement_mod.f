!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!   Module name: DISCRETELEMENT                                        C
!   Purpose: DES mod file                                              C
!            Common Block containing DEM conditions                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE DISCRETELEMENT

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param, only: dim_m
      USE param1, only: undefined, one, zero, half
      IMPLICIT NONE
!-----------------------------------------------

! Total number of particles in simulation: read from input or generated
      INTEGER :: PARTICLES

! Start particle tracking quantities
!----------------------------------------------------------------->>>
! Generally for inlet/outlet related routines but also employed in
! tracking for parallelization

! Dynamic particle states:

! NONEXISTENT: This state is used with the inlet/outlet to skip
! indices that do not represent particles in the system or indices
! that represent particles that have exited the system.

! ENTERING: This state identifies a particle as 'new' if true.
! Particles with a classification of 'new' do not react when in contact
! with a wall or another particle, however existing particles do collide
! and interact with 'new' particles. The classification allows new
! particles to push particles already in the system out of the way when
! entering to prevent overlap.  This flag is also used when the center
! of a particle crosses a dem outlet (i.e. an exiting particle; see
! EXITING) so that the particle will maintain its present trajectory
! until it has fully exited the system

! EXITING: This state identifies a particle as 'exiting' if true.
! If a particle initiates contact with a wall surface designated as a
! des outlet, this flag is set to true. With this classification the
! location of the particle is checked to assess if the particle has
! fully exited the system.  At this point, the particle is removed
! from the list.

! GHOST, ENTERING_GHOST, EXITING_GHOST: for ghost particles

      INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE :: PARTICLE_STATE ! (PARTICLES)

      INTEGER, PARAMETER :: nonexistent=0
      INTEGER, PARAMETER :: normal_particle=1
      INTEGER, PARAMETER :: entering_particle=2
      INTEGER, PARAMETER :: exiting_particle=3
      INTEGER, PARAMETER :: normal_ghost=4
      INTEGER, PARAMETER :: entering_ghost=5
      INTEGER, PARAMETER :: exiting_ghost=6

! PARALLEL PROCESSING: explanation of variables in parallel architecture
! pip - particles in each processor (includes the ghost particles)
! max_pip - maximum allocated particles in processor

! Number of particles in the system (current)
      INTEGER :: PIP
! Global sum of particles (excluding ghost) in the system
      INTEGER :: TOT_PAR
! Maximum particles permitted in the system at once
      INTEGER :: MAX_PIP

! End particle tracking quantities
!-----------------------------------------------------------------<<<

! For parallel processing: global id of particles
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGLOBAL_ID
! Ghost count of particles on each processor
      INTEGER :: IGHOST_CNT
! Maximum global id, new particles global id will be assigned based
! on this value
      Integer :: imax_global_id

! If gener_part_config is true, then the particle_input.dat file does
! not need to be supplied nor does the total number of particles as
! these are determined based on the specified volume fraction (vol_frac)
! in the specified domain
      LOGICAL :: GENER_PART_CONFIG

! Output/debug controls
!----------------------------------------------------------------->>>
! Logic that controls whether to print data dem simulations (granular or
! coupled)
      LOGICAL :: PRINT_DES_DATA
      CHARACTER(LEN=255) :: VTP_DIR

! Output file count for .vtp type files and for tecplot files;
! for vtp output used to label .vtp file names in sequential order
! and is saved so restarts begin at the correct count
      INTEGER VTP_FINDEX, TECPLOT_FINDEX
! End Output/debug controls
!-----------------------------------------------------------------<<<

! DES - Continuum
      LOGICAL DISCRETE_ELEMENT
      LOGICAL DES_CONTINUUM_COUPLED

! DES -
! With this logic the particles see the fluid but the fluid does
! not see the particles.
      LOGICAL DES_ONEWAY_COUPLED

! Collision model, options are as follows
!   linear spring dashpot model (default/undefined)
!   'hertzian' model
      CHARACTER(LEN=64) :: DES_COLL_MODEL
      INTEGER DES_COLL_MODEL_ENUM
      INTEGER,PARAMETER ::  HERTZIAN=0
      INTEGER,PARAMETER ::  LSD=1

! Integration method, options are as follows
!   'euler' first-order scheme (default)
!   'adams_bashforth' second-order scheme (by T.Li)
      CHARACTER(LEN=64) :: DES_INTG_METHOD
      LOGICAL :: INTG_ADAMS_BASHFORTH
      LOGICAL :: INTG_EULER

! Value of solids time step based on particle properties
      DOUBLE PRECISION :: DTSOLID
! Run time value of simulation time used in dem simulation
      DOUBLE PRECISION :: S_TIME


! Neighbor search related quantities
!----------------------------------------------------------------->>>
! Quantities used to determine whether neighbor search should be called
      INTEGER :: NEIGHBOR_SEARCH_N
      DOUBLE PRECISION :: NEIGHBOR_SEARCH_RAD_RATIO
      LOGICAL :: DO_NSEARCH

! Flag on whether to have DES_*_OLD arrays, if either Adams Bashforth or PIC is used
      LOGICAL :: DO_OLD

! Factor muliplied by sum of radii in grid based neighbor search and
! nsquare search method.  increases the effective radius of a particle
! for detecting particle contacts
      DOUBLE PRECISION :: FACTOR_RLM

! Stores number of neighbors based on neighbor search
      INTEGER, DIMENSION(:), ALLOCATABLE :: NEIGHBOR_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NEIGHBOR_INDEX_OLD
      INTEGER, DIMENSION(:), ALLOCATABLE :: NEIGHBORS
      INTEGER, DIMENSION(:), ALLOCATABLE :: NEIGHBORS_OLD
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PFT_NEIGHBOR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PFT_NEIGHBOR_OLD

      INTEGER :: NEIGH_NUM,NEIGH_MAX

! Quantities used for reporting: max no. neighbors and max overlap
! that exists during last solid time step of dem simulation
      DOUBLE PRECISION :: OVERLAP_MAX

! The number of i, j, k divisions in the grid used to perform the
! cell based neighbor search
      INTEGER :: DESGRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, &
                 DESGRIDSEARCH_KMAX

! End neighbor search related quantities
!-----------------------------------------------------------------<<<


! User specified dimension of the system (by default 2D, but if 3D system is
! desired then it must be explicitly specified)
      INTEGER, PARAMETER :: DIMN = 3

! X, Y, Z position of cell faces of computational fluid grid
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XE  !(0:DIMENSION_I)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: YN  !(0:DIMENSION_J)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZT  !(0:DIMENSION_K)


! Gravity vector and magnitude
      DOUBLE PRECISION :: GRAV(3)
      DOUBLE PRECISION :: GRAV_MAG


! Periodic wall BC
      LOGICAL DES_PERIODIC_WALLS
      LOGICAL DES_PERIODIC_WALLS_X
      LOGICAL DES_PERIODIC_WALLS_Y
      LOGICAL DES_PERIODIC_WALLS_Z

! Particle-particle and Particle-wall collision model parameters
!----------------------------------------------------------------->>>
! Spring contants
      DOUBLE PRECISION :: KN, KN_W  !Normal
      DOUBLE PRECISION :: KT, KT_W, KT_FAC, KT_W_FAC

! Damping coefficients
      DOUBLE PRECISION :: ETA_DES_N, ETA_N_W  !Normal
      DOUBLE PRECISION :: ETA_DES_T, ETA_T_W  !Tangential

! Tangential damping factors, eta_t = eta_t_factor * eta_N
      DOUBLE PRECISION :: DES_ETAT_FAC, DES_ETAT_W_FAC

! Damping coeffients in array form
      DOUBLE PRECISION :: DES_ETAN(DIM_M, DIM_M), DES_ETAN_WALL(DIM_M)
      DOUBLE PRECISION :: DES_ETAT(DIM_M, DIM_M), DES_ETAT_WALL(DIM_M)

! Friction coeficients
      DOUBLE PRECISION MEW, MEW_W

! coeff of restituion input in one D array, solid solid
! Tangential rest. coef. are used for hertzian collision model but not linear
      DOUBLE PRECISION :: DES_EN_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)
      DOUBLE PRECISION :: DES_ET_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)

! coeff of restitution input in one D array, solid wall
      DOUBLE PRECISION :: DES_EN_WALL_INPUT(DIM_M)
      DOUBLE PRECISION :: DES_ET_WALL_INPUT(DIM_M)

! Hertzian collision model:
      DOUBLE PRECISION :: E_YOUNG(DIM_M), Ew_YOUNG
      DOUBLE PRECISION :: V_POISSON(DIM_M), Vw_POISSON
      DOUBLE PRECISION :: HERT_KN(DIM_M, DIM_M), HERT_KWN(DIM_M)
      DOUBLE PRECISION :: HERT_KT(DIM_M, DIM_M), HERT_KWT(DIM_M)

! End particle-particle and particle-wall collision model parameters
!-----------------------------------------------------------------<<<


! Particle attributes: radius, density, mass, moment of inertia
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_RADIUS !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RO_Sol     !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PVOL       !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PMASS      !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OMOI       !(PARTICLES)

! Additional quantities
      DOUBLE PRECISION :: MIN_RADIUS, MAX_RADIUS

! Old and new particle positions, velocities (translational and
! rotational)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_NEW  !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_NEW  !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_NEW    !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PPOS         !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_ACC_OLD  !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ROT_ACC_OLD  !(PARTICLES,3)

! Defining user defined allocatable array
      INTEGER :: DES_USR_VAR_SIZE
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_USR_VAR  !(PARTICLES,3)

! Total force and torque on each particle
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FC    !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TOW   !(PARTICLES,3)

! Dynamic information related to computational (eulerian) fluid grid
!----------------------------------------------------------------->>>
! Dynamic variable. for each ijk computational fluid cell store the
! total number of particles and the id's of the particles in that cell
      TYPE iap1
         INTEGER, DIMENSION(:), POINTER:: p
      END TYPE iap1

! particle can collide with at most COLLISION_ARRAY_MAX facets simultaneously
      INTEGER :: COLLISION_ARRAY_MAX = 8

! -1 value indicates no collision
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: wall_collision_facet_id       ! (COLLISION_ARRAY_MAX,PARTICLES)
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: wall_collision_PFT ! (DIMN,COLLISION_ARRAY_MAX,PARTICLES)

! in order to facilitate the parallel processing the PIC is defined
! as single array IJK
      TYPE(iap1), DIMENSION(:,:,:), ALLOCATABLE:: pic

! Store the number of particles in a computational fluid cell
      INTEGER, DIMENSION(:), ALLOCATABLE :: PINC  ! (DIMENSION_3)

! For each particle track its i, j, k & ijk location on the fluid grid
! and solids phase no.:
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PIJK ! (PARTICLES,5)=>I,J,K,IJK,M
!-----------------------------------------------------------------<<<

! drag coefficient between gas phase and discrete particle 'phases'
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: F_GDS

! the following should probably be local to the subroutine
! solve_vel_star they are only needed when invoking the non-interpolated
! version of drag wherein they are used to determine part of the source
! in the gas-phase momentum balances, in the routine to determine
! coefficients for the pressure correction equation and in the partial
! elimination algorithm
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VXF_GDS

! the coefficient add to gas momentum A matrix
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DRAG_AM
! the coefficient add to gas momentum B matrix
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DRAG_BM

! Explicitly calculated fluid-particle drag force.
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DRAG_FC !(PARTICLES,3)

! An intermediate array used in calculation of mean solids velocity
! by backward interpolation, i.e., when INTERP_DES_MEAN_FIELDS is true.
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::DES_VEL_NODE

! An intermediate array used in calculation of solids volume fraction
! by backward interpolation, i.e., when INTERP_DES_MEAN_FIELDS is true.
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  DES_ROPS_NODE

! the gas-particle drag coefficient
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: f_gp
                        !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: wtderivp

! quantities are set in subroutine set_interpolation_scheme
! order = order of the interpolation method, ob2l = (order+1)/2,
! ob2r = order/2
      CHARACTER(LEN=7):: scheme, interp_scheme
      INTEGER:: order, ob2l, ob2r

! END of interpolation related data
!-----------------------------------------------------------------<<<

! Volume of each node. Used to obtain Eulerian fields
      double precision, allocatable, dimension(:) :: des_vol_node

! Variable to track pressure force in computational fluid cell (ijk)
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: P_FORCE

! Bulk density of particles in fluid cell
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_ROP_S

! Granular temperature in a fluid cell
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_THETA

! Flag to turn on/off optimizing the list of facets at each des grid cell
      LOGICAL :: MINIMIZE_DES_FACET_LIST

      LOGICAL :: DES_EXPLICITLY_COUPLED

! particle in cell related variable
      type iap2
         integer :: isize
         integer, dimension(:), pointer:: p
      end type iap2

      type(iap2), dimension(:),allocatable:: dg_pic
      integer, dimension(:),allocatable :: dg_pijk,dg_pijkprv

! variable to clean the ghost cells
      logical,dimension(:),allocatable :: ighost_updated
      integer :: max_isize

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!  Subroutine: DES_CROSSPRDCT                                         !
!  Purpose: Calculate the cross product of two 3D vectors.            !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      FUNCTION DES_CROSSPRDCT(XX,YY)

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Input vectors
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: XX, YY
! Result: cross product of vectors
      DOUBLE PRECISION, DIMENSION(3) :: DES_CROSSPRDCT
!......................................................................!

      DES_CROSSPRDCT(1) = XX(2)*YY(3) - XX(3)*YY(2)
      DES_CROSSPRDCT(2) = XX(3)*YY(1) - XX(1)*YY(3)
      DES_CROSSPRDCT(3) = XX(1)*YY(2) - XX(2)*YY(1)

      END FUNCTION DES_CROSSPRDCT

      END MODULE DISCRETELEMENT
