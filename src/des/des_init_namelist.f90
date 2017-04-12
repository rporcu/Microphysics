MODULE DES_INIT_NAMELIST_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS
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

      use discretelement, only: des_intg_method
      use discretelement, only: des_coll_model

      use discretelement, only: des_continuum_coupled
      use discretelement, only: des_explicitly_coupled
      use discretelement, only: des_oneway_coupled

      use discretelement, only: kn, kn_w
      use discretelement, only: kt_fac, kt_w_fac

      use discretelement, only: mew, mew_w

      use discretelement, only: des_en_input, des_en_wall_input
      use discretelement, only: des_et_input, des_et_wall_input

      use discretelement, only: e_young, ew_young
      use discretelement, only: v_poisson, vw_poisson

      use discretelement, only: des_etat_fac, des_etat_w_fac

      use discretelement, only: particles

      use discretelement, only: dim_m
      use param1, only: undefined_i, undefined


      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

      include 'desnamelist.inc'

!#####################################################################!
! DEM/PIC COMMON:      Discrete Element Simulation                    !
!#####################################################################!

!<keyword category="Discrete Element Simulation" required="false"
!  dem="true" pic="true">
!  <description>
!    To switch between pure granular or coupled simulations of carried
!    and dispersed phase flows.
!  </description>
!  <valid value=".true." note="Performs coupled simulations. "/>
      des_continuum_coupled = .FALSE.
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




      RETURN
      END SUBROUTINE DES_INIT_NAMELIST
END MODULE DES_INIT_NAMELIST_MODULE
