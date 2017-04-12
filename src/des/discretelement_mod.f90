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

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use param, only: dim_m
      use param1, only: undefined, one, zero, half
      implicit none
!-----------------------------------------------

! Total number of particles in simulation: read from input or generated
      integer :: PARTICLES

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
      integer, parameter :: nonexistent=0
      integer, parameter :: normal_particle=1
      integer, parameter :: entering_particle=2
      integer, parameter :: exiting_particle=3
      integer, parameter :: normal_ghost=4
      integer, parameter :: entering_ghost=5
      integer, parameter :: exiting_ghost=6


! Output/debug controls
!----------------------------------------------------------------->>>
! Logic that controls whether to print data dem simulations (granular or
! coupled)
      logical :: PRINT_DES_DATA
      CHARACTER(LEN=255) :: VTP_DIR

! Output file count for .vtp type files and for tecplot files;
! for vtp output used to label .vtp file names in sequential order
! and is saved so restarts begin at the correct count
      integer :: VTP_FINDEX
! End Output/debug controls
!-----------------------------------------------------------------<<<

! DES - Continuum
      logical :: DES_CONTINUUM_COUPLED
      logical :: DES_EXPLICITLY_COUPLED

! With this logic the particles see the fluid but the fluid does
! not see the particles.
      logical :: DES_ONEWAY_COUPLED

! Collision model, options are as follows
!   linear spring dashpot model (default/undefined)
!   'hertzian' model
      CHARACTER(LEN=64) :: DES_COLL_MODEL
      integer :: DES_COLL_MODEL_ENUM
      integer,PARAMETER :: HERTZIAN=0
      integer,PARAMETER :: LSD=1

! Integration method, options are as follows
!   'euler' first-order scheme (default)
!   'adams_bashforth' second-order scheme (by T.Li)
      CHARACTER(LEN=64) :: DES_INTG_METHOD
      logical :: INTG_ADAMS_BASHFORTH
      logical :: INTG_EULER

! Value of solids time step based on particle properties
      real(c_real) :: DTSOLID
! Run time value of simulation time used in dem simulation
      real(c_real) :: S_TIME


! Neighbor search related quantities
!----------------------------------------------------------------->>>
! Quantities used to determine whether neighbor search should be called

      logical :: DO_NSEARCH

! Flag on whether to have DES_*_OLD arrays, if either Adams Bashforth or PIC is used
      logical :: DO_OLD

! Quantities used for reporting: max no. neighbors and max overlap
! that exists during last solid time step of dem simulation
      real(c_real) :: OVERLAP_MAX

! End neighbor search related quantities
!-----------------------------------------------------------------<<<


! Particle-particle and Particle-wall collision model parameters
!----------------------------------------------------------------->>>
! Spring contants
      real(c_real) :: KN, KN_W  !Normal
      real(c_real) :: KT, KT_W, KT_FAC, KT_W_FAC

! Damping coefficients
      real(c_real) :: ETA_DES_N, ETA_N_W  !Normal
      real(c_real) :: ETA_DES_T, ETA_T_W  !Tangential

! Tangential damping factors, eta_t = eta_t_factor * eta_N
      real(c_real) :: DES_ETAT_FAC, DES_ETAT_W_FAC

! Damping coeffients in array form
      real(c_real) :: DES_ETAN(DIM_M, DIM_M), DES_ETAN_WALL(DIM_M)
      real(c_real) :: DES_ETAT(DIM_M, DIM_M), DES_ETAT_WALL(DIM_M)

! Friction coefficients
      real(c_real) :: MEW, MEW_W

! coeff of restituion input in one D array, solid solid
! Tangential rest. coef. are used for hertzian collision model but not linear
      real(c_real) :: DES_EN_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)
      real(c_real) :: DES_ET_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)

! coeff of restitution input in one D array, solid wall
      real(c_real) :: DES_EN_WALL_INPUT(DIM_M)
      real(c_real) :: DES_ET_WALL_INPUT(DIM_M)

! Hertzian collision model:
      real(c_real) :: E_YOUNG(DIM_M), Ew_YOUNG
      real(c_real) :: V_POISSON(DIM_M), Vw_POISSON
      real(c_real) :: HERT_KN(DIM_M, DIM_M), HERT_KWN(DIM_M)
      real(c_real) :: HERT_KT(DIM_M, DIM_M), HERT_KWT(DIM_M)

! End particle-particle and particle-wall collision model parameters
!-----------------------------------------------------------------<<<

! Additional quantities
      real(c_real) :: MIN_RADIUS, MAX_RADIUS

!-----------------------------------------------------------------<<<





      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!  Subroutine: DES_CROSSPRDCT                                         !
!  Purpose: Calculate the cross product of two 3D vectors.            !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      FUNCTION DES_CROSSPRDCT(XX,YY)

      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! Input vectors
      real(c_real), DIMENSION(3), INTENT(IN) :: XX, YY
! Result: cross product of vectors
      real(c_real), DIMENSION(3) :: DES_CROSSPRDCT
!......................................................................!

      DES_CROSSPRDCT(1) = XX(2)*YY(3) - XX(3)*YY(2)
      DES_CROSSPRDCT(2) = XX(3)*YY(1) - XX(1)*YY(3)
      DES_CROSSPRDCT(3) = XX(1)*YY(2) - XX(2)*YY(1)

      END FUNCTION DES_CROSSPRDCT

      END MODULE DISCRETELEMENT
