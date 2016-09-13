!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          C
!                                                                         C
!     Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!     Comments: Added some interpolation based inputs                     C
!                                                                         C
!  Keyword Documentation Format:                                          C
!<keyword category="category name" required="true/false"                  C
!                                    legacy="true/false">                 C
!  <description></description>                                            C
!  <arg index="" id="" max="" min=""/>                                    C
!  <dependent keyword="" value="DEFINED"/>                                C
!  <conflict keyword="" value="DEFINED"/>                                 C
!  <valid value="" note="" alias=""/>                                     C
!  <range min="" max="" />                                                C
!  MFIX_KEYWORD=INIT_VALUE                                                C
!</keyword>                                                               C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_INIT_NAMELIST

      USE param1
      USE discretelement
      USE des_bc
      USE particle_filter

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

      INCLUDE 'desnamelist.inc'



!#####################################################################!
!                             Run Control                             !
!#####################################################################!




!#####################################################################!
!                           Physical Parameters                       !
!#####################################################################!



!#####################################################################!
!                          Numerical Parameters                       !
!#####################################################################!



!#####################################################################!
!                          Output Control                             !
!#####################################################################!

!<keyword category="Output Control" required="false"
!  dem="true" pic="true">
!  <description>
!    Reports mass based on Lagrangian particles and continuum
!    representation. Useful to ensure mass conservation between
!    Lagrangian and continuum representations. Recommended use for
!    debugging purposes.
!  </description>
!  <dependent keyword="DES_INTERP_MEAN_FIELDS" value=".TRUE."/>
      DES_REPORT_MASS_INTERP = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false"
!  dem="true" pic="true">
!  <description>
!    Allows writing of discrete particle data to output files. Relevant
!    to both granular and coupled simulations.
!  </description>
      PRINT_DES_DATA = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false"
!   dem="true" pic="true">
!  <description> The output file format for DES data.</description>
!  <valid value="PARAVIEW" note="ParaView formatted files (.vtp)"/>
!  <valid value="TECPLOT" note="Tecplot formatted files (.dat)"/>
      DES_OUTPUT_TYPE = "PARAVIEW"
!</keyword>



!#####################################################################!
! DEM/PIC COMMON:      Discrete Element Simulation                    !
!#####################################################################!


!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    Number of particles to be read in from the particle_input.dat file.
!    This value is overwritten when using automatic particle generation.
!    A simulation with a mass inflow BC can start without solids by
!    setting PARTICLES = 0.
!  </description>
!  <range min="0" max="+Inf" />
      PARTICLES = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    Automatically generate the initial particle position and velocity
!    data based on the parameters specified for each initial condition
!    (IC) region.
!  </description>
!  <valid value=".TRUE." note="Generate particle configuration based
!    on the initial condition parameters. Data provided in the
!    particle_input.dat file, if present, is ignored. "/>
!  <valid value=".FALSE." note="Particle position and velocity data are
!    provided in the particle_input.dat file. A runtime error occurs if
!    this file is not provided."/>
      GENER_PART_CONFIG = .FALSE.
!</keyword>

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    To switch between pure granular or coupled simulations of carried
!    and dispersed phase flows.
!  </description>
!  <valid value=".true." note="Performs coupled simulations. "/>
      DES_CONTINUUM_COUPLED = .FALSE.
!</keyword>

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>Run one-way coupled simulations. The fluid does not
! see the particles in terms of drag force. The effect of particle volume
! is still felt by the fluid through non-unity voidage values.
! </description>
      DES_ONEWAY_COUPLED = .FALSE.
!</keyword>

!<keyword category="Discrete Element Simulation" required="false" dem="true">
!  <description>
!    Time stepping scheme.
!  </description>
!  <valid value="EULER"
!    note="First-Order Euler Scheme."/>
!  <valid value="ADAMS BASHFORTH"
!    note="Second order ADAMS BASHFORTH scheme (DEM only)"/>
      DES_INTG_METHOD = 'EULER'
!</keyword>

!<keyword category="Discrete Element Simulation" required="false" dem="true">
!  <description>
!    Defines the size of the particle-based user variable:
!    DES_USR_VAR(SIZE, PARTICLES). Information in this array follows
!    the particle throughout a simulation.
!  </description>
      DES_USR_VAR_SIZE = 0
!</keyword>

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    Number of des grid cells in the I-direction. If left undefined,
!    then it is set by MFIX such that its size equals three times the
!    maximum particle diameter with a minimum of 1 cell.
!  </description>
      DESGRIDSEARCH_IMAX = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    Number of des grid cells in the J-direction. If left undefined,
!    then it is set by MFIX such that its size equals three times
!    the maximum particle diameter with a minimum of 1 cell.
!  </description>
      DESGRIDSEARCH_JMAX = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    Number of des grid cells in the K-direction. If left undefined,
!    then it is set by MFIX such that its size equals three times
!    the maximum particle diameter with a minimum of 1 cell.
!  </description>
      DESGRIDSEARCH_KMAX = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    Specify the scheme used to map data to/from a particle's position
!    and the Eulerian grid. This keyword is required when
!    DES_INTERP_MEAN_FIELDS and/or DES_INTERP_ON are specified. A
!    graphical representation of the schemes is shown below.
!  </description>
!  <valid value="NONE" note="Do not use interpolation."/>
!  <valid value="GARG_2012" note="Interpolate to/from a particle's
!    position using the corners (nodes) of the fluid cells. This was
!    the default behavior prior to the 2015-1 Release.
!    See Garg et al. (2012) Documentation of the open-souce MFIX-DEM
!    software for gas-solids flows."/>
!  <valid value="SQUARE_DPVM" note="Divided Particle Volume Method:
!    Information is interpolated to/from a particles position using
!    a square filter of size DES_INTERP_WIDTH. This scheme is not
!    available to MFIX-PIC simulations."/>
      DES_INTERP_SCHEME = 'NONE'
!</keyword>

!<keyword category="Discrete Element Simulation" required="false" dem="true">
!  <description>
!    The length used in interpolating data to/from a particle's position
!    and the Eulerian grid. The interpolation width is only applicable
!    to the DPVM_SQUARE and DPVM_GAUSS interpolation schemes as the
!    GARG_2012 scheme's interpolation width is determined by the
!    Eulerian grid dimensions.
!    o The interpolation half-width cannot exceed the minimum cell
!      dimension because interpolation is restricted to the 27-cell
!      neighborhood surrounding a particle (9-cell neighborhood in 2D).
!    o It is recommend that the DES_INTERP_WIDTH be set equal to the
!      maximum particle diameter when using STL defined boundaries.
!  </description>
      DES_INTERP_WIDTH = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    Enables/Disables interpolation of field quantities to a particle's
!    position. This is used in calculating gas-particle interactions,
!    such as the drag force.
!  </description>
!  <valid value=".FALSE." note="Use fluid values from the cell containing
!    the particle's center."/>
!  <valid value=".TRUE." note="Interpolate fluid values from the 27-cell
!    neighborhood to a particle's position."/>
      DES_INTERP_ON = .FALSE.
!</keyword>

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    Enables/Disables interpolation of particle data (e.g., solids
!    volume and drag force) from a particle's position to the
!    Eulerian grid.
!  </description>
!  <valid value=".FALSE." note="Assign particle data to the fluid
!    grid cell containing the particle's center."/>
!  <valid value=".TRUE." note="Interpolate particle data from the
!    particle's position to the 27-cell neighborhood surrounding
!    the particle."/>
      DES_INTERP_MEAN_FIELDS = .FALSE.
!</keyword>


!<keyword category="Discrete Element Simulation" required="false" dem="true">
!  <description>
!    Enable/Disable explicit coupling of DEM solids and the fluid. This
!    algorithm is presently limited to hydrodynamic simulations.
!  </description>
!  <valid value=".FALSE." note="The fluid and particles calculate
!    interphase forces at their respective time scales. The fluid phase
!    calculates the interphase coupling forces once per fluid time step.
!    Similarly, DEM particles calculate the interface coupling forces at
!    each solids time-step. The DEM must also bin particles to the fluid
!    grid and recalculate the fluid volume fraction every time-step."/>
!  <valid value=".TRUE." note="Interphase forces are calculated during
!    the fluid time step and stored for each particle. The interphase
!    forces are then distributed among the solids time-steps. This
!    approach can substantially reduce the computational overhead for
!    coupled simulations."/>
      DES_EXPLICITLY_COUPLED = .FALSE.
!</keyword>


!#####################################################################!
! DEM ONLY:            Discrete Element Model                         !
!#####################################################################!

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Maximum number of steps through a DEM loop before a neighbor
!    search will be performed. The search may be called earlier
!    based on other logic.
!  </description>
!  <range min="0.0" max="+Inf" />
      NEIGHBOR_SEARCH_N = 25
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Flag to set the neighbor search algorithm.
!  </description>
!  <valid value="1" note="N-Square search algorithm (most expensive)"/>
!  <valid value="4" note="Grid-Based Neighbor Search (Recommended)"/>
      DES_NEIGHBOR_SEARCH = 4
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Ratio of the distance (imaginary sphere radius) to particle radius
!    that is allowed before a neighbor search is performed. This works
!    in conjunction with the logic imposed by NEIGHBOR_SEARCH_N in
!    deciding calls to the neighbor search algorithm.
!  </description>
      NEIGHBOR_SEARCH_RAD_RATIO = 1.0D0
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Effectively increase the radius of a particle (multiple of the sum
!    of particle radii) during the building of particle neighbor list.
!  </description>
      FACTOR_RLM = 1.2
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Collision model for the soft-sphere approach used in DEM model.
!    All models require specifying the following parameters: DES_EN_INPUT,
!    DES_EN_WALL_INPUT, MEW, and MEW_W.
!  </description>
!  <valid value="LSD" note="The linear spring-dashpot model.
!    Requires: KN, KN_W, KT_FAC, KT_W_FAC, DES_ETAT_FAC, DES_ETAT_W_FAC."/>
!  <valid value="HERTZIAN" note="The Hertzian model.
!    Requires: DES_ET_INPUT, DES_ET_WALL_INPUT, E_YOUNG, EW_YOUNG
!    V_POISSON, VW_POISSON."/>
      DES_COLL_MODEL = 'LSD'
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Normal spring constant [dyne/cm in CGS] for inter-particle collisions.
!    Required when using the linear spring-dashpot collision model.
!  </description>
      KN = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Ratio of the tangential spring constant to normal spring constant
!    for inter-particle collisions. Use it to specify the tangential
!    spring constant for particle-particle collisions as KT_FAC*KN.
!    Required when using the linear spring-dashpot collision model.
!  </description>
!  <dependent keyword="DES_COLL_MODEL" value="LSD"/>
!  <range min="0.0" max="1.0" />
      KT_FAC = 2.d0/7.d0
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem=.true.>
!  <description>
!    Normal spring constant [dyne/cm in CGS] for particle-wall collisions.
!    Required when using the linear spring-dashpot collision model.
!  </description>
      KN_W = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Ratio of the tangential spring constant to normal spring constant
!    for particle-wall collisions. Use it to specify the tangential
!    spring constant for particle-wall collisions as KT_W_FAC*KN_W.
!    Required when using the linear spring-dashpot collision model.
!  </description>
!  <dependent keyword="DES_COLL_MODEL" value="LSD"/>
!  <range min="0.0" max="1.0" />
      KT_W_FAC = 2.d0/7.d0
!</keyword>

!<keyword category="Discrete Element Model" required="false" dem="true"
!  <description>
!    Inter-particle Coulomb friction coefficient.
!  </description>
! <range min="0.0" max="1.0" />
      MEW = UNDEFINED
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Particle-wall Coulomb friction coefficient.
!  </description>
! <range min="0.0" max="1.0" />
      MEW_W = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    The normal restitution coefficient for inter-particle collisions
!    used to determine the inter-particle normal damping factor.
!
!    Values should be defined for a single dimensional array. For
!    example, a simulation with three solids phases (MMAX=3) needs
!    six values: en11, en12, en13; en22 en 23; en33.
!  </description>
!  <range min="0.0" max="1.0" />
      DES_EN_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    The normal restitution coefficient for particle-wall collisions
!    used to determine the particle-wall normal damping factor.
!
!    Values should be defined in a single dimensional array. For
!    example, a simulation with three solids phases (MMAX=3) needs
!    three values: enw1, enw2, enw3.
!  </description>
!  <range min="0.0" max="1.0" />
      DES_EN_WALL_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Tangential restitution coefficient for inter-particle collisions.
!    Values are defined in a one dimensional array. This is required
!    input when using the Hertzian collision model.
! </description>
! <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
! <range min="0.0" max="1.0" />
      DES_ET_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Tangential restitution coefficient for particle wall collisions.
!    Values are defined in a one dimensional array. This is required
!    input when using the Hertzian collision model.
!  </description>
! <range min="0.0" max="1.0" />
! <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      DES_ET_WALL_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Ratio of the tangential damping factor to the normal damping factor
!    for inter-particle collisions.  Required for the linear spring-
!    dashpot model collision model
!  </description>
!  <dependent keyword="DES_COLL_MODEL" value="LSD"/>
!  <range min="0.0" max="1.0" />
!  <valid value="UNDEFINED" note="For LSD model, if left undefined, MFIX
!   reverts to default value of 0.5" />
      DES_ETAT_FAC = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
! <description>
!    Ratio of the tangential damping factor to the normal damping
!    factor for particle-wall collisions. Required for the linear
!    spring-dashpot model for soft-spring collision modelling under
!    DEM. For the Hertzian model, the tangential damping coefficients
!    have to be explicitly specified and specification of this
!    variable is not required.
! </description>
! <dependent keyword="DES_COLL_MODEL" value="LSD"/>
! <range min="0.0" max="1.0" />
! <valid value="UNDEFINED" note="For LSD model, if left undefined, MFIX
! will revert to default value of 0.5" />
      DES_ETAT_W_FAC = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Youngs modulus for the wall [barye in CGS]. Required when using the
!    Hertzian spring-dashpot model.
!  </description>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      EW_YOUNG = UNDEFINED
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Poisson ratio for the wall. Required when using the Hertzian
!    spring-dashpot model.
!  </description>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      VW_POISSON = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Youngs modulus for the particle [barye in CGS]. Required when using
!    the Hertzian spring-dashpot model.
!  </description>
!  <arg index="1" id="Phase" min="1" max="MMAX"/>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      E_YOUNG(:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Poissons ratio for the particle. Required when using the Hertzian
!    spring-dashpot model.
!  </description>
!  <arg index="1" id="Phase" min="1" max="MMAX"/>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      V_POISSON(:DIM_M) = UNDEFINED
!</keyword>



!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Flag to turn on/off optimizing the list of facets at each des grid cell
!  </description>
!  <dependent keyword="USE_STL" value=".TRUE."/>
      MINIMIZE_DES_FACET_LIST =.TRUE.
!</keyword>


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                            UNSUPPORTED KEYWORDS                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

! Logical to force the inlet to operate with an ordered boundary
! condition. This may be useful during long simulations or if the
! inlet appears to be taking a long time to randomly place particles.
      FORCE_ORD_BC = .FALSE.


! These need to be inialized to 0, but they are not part of the namelist
      VTP_FINDEX = 0
      TECPLOT_FINDEX = 0

      RETURN
      END SUBROUTINE DES_INIT_NAMELIST
