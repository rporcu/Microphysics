!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: INIT_NAMELIST                                           !
!  Purpose: initialize the NAMELIST variables                          !
!                                                                      !
!  Author: P. Nicoletti                               Date: 26-NOV-91  !
!                                                                      !
!  Keyword Documentation Format:                                       !
!                                                                      !
!<keyword category="category name" required="true"/FALSE               !
!                                    legacy=TRUE/FALSE>                !
!  <description></description>                                         !
!  <arg index="" id="" max="" min=""/>                                 !
!  <dependent keyword="" value="DEFINED"/>                             !
!  <conflict keyword="" value="DEFINED"/>                              !
!  <valid value="" note="" alias=""/>                                  !
!  <range min="" max="" />                                             !
!  MFIX_KEYWORD=INIT_VALUE                                             !
!</keyword>                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE INIT_NAMELIST

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run
      USE output
      USE physprop
      USE geometry
      USE ic
      USE bc
      USE ps
      USE fldvar
      USE constant
      USE indices
      USE toleranc
      USE scales
      USE ur_facs
      USE leqsol
      USE residual
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counters
      INTEGER :: LC

!#####################################################################!
!                             Run Control                             !
!#####################################################################!

!<keyword category="Run Control" required="true">
!  <description> Name used to create output files. The name should
!    generate legal file names after appending extensions.
!    Ex: Given the input, RUN_NAME = "bub01", MFIX will generate
!    the output files: BUB01.LOG, BUB01.OUT, BUB01.RES, etc.
!  </description>
      RUN_NAME = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Problem description. Limited to 60 characters.</description>
      DESCRIPTION = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="true">
!  <description> Simulation input/output units.</description>
!  <valid value="cgs" note="All input and output in CGS units (g, cm, s, cal)."/>
!  <valid value="si" note="All input and output in SI units (kg, m, s, J)."/>
      UNITS = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="true">
!  <description>Type of run.</description>
!  <valid value="new" note="A new run. There should be no .RES, .SPx,
!    .OUT, or .LOG files in the run directory."/>
!  <valid value="RESTART_1" note="Traditional restart. The run continues
!    from the last time the .RES file was updated and new data is added
!    to the SPx files."/>
!  <valid value="RESTART_2"
!    note="Start a new run with initial conditions from a .RES file
!      created from another run. No other data files (SPx) should be
!      in the run directory."/>
      RUN_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Simulation start time. This is typically zero.
!  </description>
!  <range min="0.0" max="+Inf" />
      TIME = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Simulation stop time.
!  </description>
!  <range min="0.0" max="+Inf" />
      TSTOP = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Initial time step size. If left undefined, a steady-state
!    calculation is performed.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Maximum time step size.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MAX = ONE
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Minimum time step size.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MIN = 1.0D-6
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Factor for adjusting time step.
!    * The value must be less than or equal to 1.0.
!    * A value of 1.0 keeps the time step constant which may help overcome
!      initial non-convergence.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="1" />
      DT_FAC = 0.9D0
!</keyword>


!<keyword category="Run Control" required="false">
!  <description>
!    Flag to restart the code when DT < DT_MIN.
!  </description>
      AUTO_RESTART = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to enable/disable solving the X-momentum equations.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve X-momentum equations."/>
!  <valid value=".FALSE." note="The X velocity initial conditions
!   persist throughout the simulation."/>
      MOMENTUM_X_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to enable/disable solving the Y-momentum equations.
! </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve Y-momentum equations."/>
!  <valid value=".FALSE." note="The Y velocity initial conditions
!   persist throughout the simulation."/>
      MOMENTUM_Y_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to enable/disable solving the Z-momentum equations.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve Z-momentum equations."/>
!  <valid value=".FALSE." note="The Z velocity initial conditions
!   persist throughout the simulation."/>
      MOMENTUM_Z_EQ(:DIM_M) = .TRUE.
!</keyword>



!<keyword category="Run Control" required="false">
!  <description>
!     Available gas-solids drag models.
!     Note: The extension _PCF following the specified drag model
!     indicates that the polydisperse correction factor is available.
!     For PCF details see:
!     o Van der Hoef MA, Beetstra R, Kuipers JAM. (2005)
!       Journal of Fluid Mechanics.528:233-254.
!     o Beetstra, R., van der Hoef, M. A., Kuipers, J.A.M. (2007).
!       AIChE Journal, 53:489-501.
!     o Erratum (2007), AIChE Journal, Volume 53:3020
!  </description>
!
!  <valid value="SYAM_OBRIEN" note="Syamlal M, OBrien TJ (1988).
!   International Journal of Multiphase Flow 14:473-481.
!   Two additional parameters may be specified: DRAG_C1, DRAG_D1"/>
!
!  <valid value="GIDASPOW" note="Ding J, Gidaspow D (1990).
!   AIChE Journal 36:523-538"/>
!
!  <valid value="GIDASPOW_BLEND" note="Lathouwers D, Bellan J (2000).
!    Proceedings of the 2000 U.S. DOE
!        Hydrogen Program Review NREL/CP-570-28890."/>
!
!  <valid value="WEN_YU" note="Wen CY, Yu YH (1966).
!   Chemical Engineering Progress Symposium Series 62:100-111."/>
!
!  <valid value="KOCH_HILL" note="Hill RJ, Koch DL, Ladd JC (2001).
!   Journal of Fluid Mechanics, 448: 213-241. and 448:243-278."/>
!
!  <valid value="BVK" note="Beetstra, van der Hoef, Kuipers (2007).
!   Chemical Engineering Science 62:246-255"/>
!
!  <valid value="USER_DRAG" note="Invoke user-defined drag law. (usr_drag.f)"/>
!
!  <valid value="GIDASPOW_PCF" note="see GIDASPOW"/>
!  <valid value="GIDASPOW_BLEND_PCF" note="see GIDASPOW_BLEND"/>
!  <valid value="WEN_YU_PCF" note="see WEN_YU"/>
!  <valid value="KOCH_HILL_PCF" note="see KOCH_HILL"/>
!
      DRAG_TYPE = 'SYAM_OBRIEN'
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Quantity for calibrating Syamlal-O'Brien drag correlation using Umf
!    data.  This is determined using the Umf spreadsheet.
!  </description>
      DRAG_C1 = 0.8d0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Quantity for calibrating Syamlal-O'Brien drag correlation using Umf
!    data.  This is determined using the Umf spreadsheet.
!  </description>
      DRAG_D1 = 2.65d0
!</keyword>


!<keyword category="Run Control" required="false" tfm="true">
!  <description>
!    Subgrid models.
!  </description>
!
!  <valid value="Igci" note="
!   Igci, Y., Pannala, S., Benyahia, S., and Sundaresan S. (2012).
!   Industrial & Engineering Chemistry Research, 2012, 51(4):2094-2103"/>
!
!  <valid value="Milioli" note="
!   Milioli, C.C., Milioli, F. E., Holloway, W., Agrawal, K. and
!   Sundaresan, S. (2013). AIChE Journal, 59:3265-3275."/>
!
      SUBGRID_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false" tfm="true">
!  <description>
!    Ratio of filter size to computational cell size.
!  </description>
      FILTER_SIZE_RATIO = 2.0D0
!</keyword>

!<keyword category="Run Control" required="false" tfm="true">
!  <description>Flag for subgrid wall correction.</description>
!  <valid value=".FALSE." note="Do not include wall correction."/>
!  <valid value=".TRUE." note="Include subgrid wall correction."/>
      SUBGRID_Wall = .FALSE.
!</keyword>


!#####################################################################!
!                           Physical Parameters                       !
!#####################################################################!


!<keyword category="Physical Parameters" required="false">
!  <description>Reference pressure. [0.0]</description>
      P_REF = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>Scale factor for pressure. [1.0]</description>
      P_SCALE = ONE
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>Gravitational acceleration. [980.7 in CGS]</description>
      GRAVITY = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    X-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_X = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Y-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_Y = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Z-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_Z = ZERO
!</keyword>





!#####################################################################!
!                          Numerical Parameters                       !
!#####################################################################!



!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum number of iterations [500].
!  </description>
      MAX_NIT = 500
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Factor to normalize the gas continuity equation residual. The
!    residual from the first iteration is used if NORM_G is left
!    undefined. NORM_G=0 invokes a normalization method based on the
!    dominant term in the continuity equation. This setting may speed up
!    calculations, especially near a steady state and incompressible
!    fluids. But, the number of iterations for the gas phase pressure
!    should be increased, LEQ_IT(1), to ensure mass balance
!  </description>
      NORM_G = 1.0d0
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Maximum residual at convergence (Continuity + Momentum) [1.0d-3].
!  </description>
      TOL_RESID = 1.0D-3
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Minimum residual for declaring divergence [1.0d+4].
!    This parameter is useful for incompressible fluid simulations
!    because velocity residuals can take large values for the second
!    iteration (e.g., 1e+8) before dropping down to smaller values for
!    the third iteration.
!  </description>
      TOL_DIVERGE = 1.0D+4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Reduce the time step if the residuals stop decreasing. Disabling this
!    feature may help overcome initial non-convergence.
!  </description>
!  <valid value=".FALSE." note="Continue iterating if residuals stall."/>
!  <valid value=".TRUE."  note="Reduce time step if residuals stall."/>
      DETECT_STALL = .TRUE.
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>
!    LEQ Solver selection. BiCGSTAB is the default method for all
!    equation types.
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value="1" note="SOR - Successive over-relaxation"/>
!  <valid value="2" note="BiCGSTAB - Biconjugate gradient stabilized."/>
!  <valid value="3" note="GMRES - Generalized minimal residual method"/>
!  <valid value="5" note="CG - Conjugate gradient"/>
      LEQ_METHOD(:) = 2
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Linear Equation tolerance [1.0d-4].
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <dependent keyword="LEQ_METHOD" value="2"/>
!  <dependent keyword="LEQ_METHOD" value="3"/>
      LEQ_TOL(:) = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Number of iterations in the linear equation solver.
!    o 20 iterations for equation types 1-2
!    o  5 iterations for equation types 3-5,10
!    o 15 iterations for equation types 6-9
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
      LEQ_IT(1) =  20
      LEQ_IT(2) =  20
      LEQ_IT(3) =   5
      LEQ_IT(4) =   5
      LEQ_IT(5) =   5
      LEQ_IT(6) =  15
      LEQ_IT(7) =  15
      LEQ_IT(8) =  15
      LEQ_IT(9) =  15
      LEQ_IT(10) =  5
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Linear equation sweep direction. This applies when using GMRES or
!    when using the LINE preconditioner with BiCGSTAB or CG methods.
!    'RSRS' is the default for all equation types.
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value="RSRS" note="(Red/Black Sweep, Send Receive) repeated twice"/>
!  <valid value="ISIS" note="(Sweep in I, Send Receive) repeated twice"/>
!  <valid value="JSJS" note="(Sweep in J, Send Receive) repeated twice"/>
!  <valid value="KSKS" note="(Sweep in K, Send Receive) repeated twice"/>
!  <valid value="ASAS" note="(All Sweep, Send Receive) repeated twice"/>
      LEQ_SWEEP(:) = 'RSRS'
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Linear precondition used by the BiCGSTAB and CG LEQ solvers. 'LINE'
!    is the default for all equation types.
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value="NONE" note="No preconditioner"/>
!  <valid value="LINE" note="Line relaxation"/>
!  <valid value="DIAG" note="Diagonal Scaling"/>
      LEQ_PC(:) = 'LINE'
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Under relaxation factors.
!    o 0.8 for equation types 1,9
!    o 0.5 for equation types 2,3,4,5,8
!    o 1.0 for equation types 6,7,10
!  </description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
      UR_FAC(1)  = 0.8D0     ! pressure
      UR_FAC(2)  = 0.5D0     ! rho, ep
      UR_FAC(3)  = 0.5D0     ! U
      UR_FAC(4)  = 0.5D0     ! V
      UR_FAC(5)  = 0.5D0     ! W
      UR_FAC(6)  = 1.0D0     ! T
      UR_FAC(7)  = 1.0D0     ! X
      UR_FAC(8)  = 0.5D0     ! Th
      UR_FAC(10) = 1.0D0     ! DES Diffusion
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>Discretization scheme of equations.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="DIM_EQS"/>
!  <valid value="0" note="First-order upwinding."/>
!  <valid value="1" note="First-order upwinding (using down-wind factors)."/>
!  <valid value="3" note="Smart."/>
!  <valid value="2" note="Superbee (recommended method)."/>
!  <valid value="5" note="QUICKEST (does not work)."/>
!  <valid value="4" note="ULTRA-QUICK."/>
!  <valid value="7" note="van Leer."/>
!  <valid value="6" note="MUSCL."/>
!  <valid value="8" note="minmod."/>
!  <valid value="9" note="Central (often unstable; useful for testing)."/>
      DISCRETIZE(:) = 0
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>
!    The code declares divergence if the velocity anywhere in the domain
!    exceeds a maximum value.  This maximum value is automatically
!    determined from the boundary values. The user may scale the maximum
!    value by adjusting this scale factor [1.0d0].
!  </description>
      MAX_INLET_VEL_FAC = ONE
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Solve transpose of linear system. (BICGSTAB ONLY).
!  </description>
!  <dependent keyword="LEQ_METHOD" value="2"/>
      DO_TRANSPOSE = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Frequency to check for convergence. (BICGSTAB ONLY)
!  </description>
!  <dependent keyword="LEQ_METHOD" value="2"/>
      icheck_bicgs = 1
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Sets optimal LEQ flags for parallel runs.
!  </description>
      OPT_PARALLEL = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Use do-loop assignment over direct vector assignment.
!  </description>
      USE_DOLOOP = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    Calculate dot-products more efficiently (Serial runs only.)
!  </description>
      IS_SERIAL = .TRUE.
!</keyword>


!#####################################################################!
!                      Geometry and Discretization                    !
!#####################################################################!


!<keyword category="Geometry and Discretization" required="false">
!  <description>Coordinates used in the simulation.</description>
!  <valid value="cartesian" note="Cartesian coordinates."/>
      COORDINATES = UNDEFINED_C
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>(Do not use.)</description>
!  <valid value=".FALSE." note="x (r) direction is considered."/>
!  <valid value=".TRUE." note="x (r) direction is not considered."/>
!     NO_I = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Number of cells in the x (r) direction.</description>
      IMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the x (r) direction. Enter values from DX(0) to
!    DX(IMAX-1).
!    o Use uniform mesh size with higher-order discretization methods.
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_I"/>
      DX(:DIM_I) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description> do not use  </description>
      XMIN = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Reactor length in the x (r) direction.</description>
      XLENGTH = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>(Do not use.)</description>
!  <valid value=".FALSE. note="y-direction is considered."/>
!  <valid value=".TRUE." note="y-direction is not considered."/>
!     NO_J = .FALSE.
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>Number of cells in the y-direction.</description>
      JMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the y-direction. Enter values from DY(0) to
!    DY(IMAX-1). Use uniform mesh size with second-order
!    discretization methods.
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_J"/>
      DY(:DIM_J) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Reactor length in the y-direction.</description>
      YLENGTH = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag to disable the third dimension (i.e., 2D simulation).
!      o Z axis in Cartesian coordinate system
!  </description>
!  <valid value=".FALSE." note="3D simulation."/>
!  <valid value=".TRUE."  note="2D simulation."/>
      NO_K = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Number of cells in the z-direction.</description>
      KMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the z (theta) direction. Enter values from DZ(0) to
!    DZ(IMAX-1). Use uniform mesh size with second-order discretization
!    methods.
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_K"/>
      DZ(:DIM_K) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>Reactor length in the z (theta) direction.</description>
      ZLENGTH = UNDEFINED
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the x-direction cyclic without pressure drop. No other
!    boundary conditions for the x-direction should be specified.
!</description>
!  <valid value=".FALSE." note="No cyclic condition at x-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition at x-boundary."/>
      CYCLIC_X = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the x-direction cyclic with pressure drop. If the
!    keyword FLUX_G is given a value this becomes a cyclic boundary
!    condition with specified mass flux. No other boundary conditions
!    for the x-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at x-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition with pressure drop at x-boundary."/>
      CYCLIC_X_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across XLENGTH when a cyclic boundary condition
!    with pressure drop is imposed in the x-direction.
!  </description>
      DELP_X = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the y-direction cyclic without pressure drop. No
!    other boundary conditions for the y-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at y-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition at x-boundary."/>
      CYCLIC_Y = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the y-direction cyclic with pressure drop. If the
!    keyword FLUX_G is given a value this becomes a cyclic boundary
!    condition with specified mass flux. No other boundary conditions
!    for the y-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at y-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition with pressure drop at y-boundary."/>
      CYCLIC_Y_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across YLENGTH when a cyclic boundary condition
!    with pressure drop is imposed in the y-direction.
!  </description>
      DELP_Y = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the z-direction cyclic without pressure drop. No
!    other boundary conditions for the z-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at z-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition at z-boundary."/>
      CYCLIC_Z = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the z-direction cyclic with pressure drop. If the
!    keyword FLUX_G is given a value this becomes a cyclic boundary
!    condition with specified mass flux. No other boundary conditions
!    for the z-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="No cyclic condition at z-boundary."/>
!  <valid value=".TRUE." note="Cyclic condition with pressure drop at
!    z-boundary."/>
      CYCLIC_Z_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across ZLENGTH when a cyclic boundary condition
!    with pressure drop is imposed in the z-direction.
!  </description>
      DELP_Z = UNDEFINED
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    If a value is specified (in units of g/cm^2.s), the domain-averaged gas
!    flux is held constant at that value in simulations over a periodic
!    domain.  A pair of boundaries specified as periodic with fixed
!    pressure drop is then treated as periodic with fixed mass flux.
!    Even for this case a pressure drop must also be specified, which
!    is used as the initial guess in the simulations.
!  </description>
      Flux_g = UNDEFINED
!</keyword>


!#####################################################################!
!                               Gas Phase                             !
!#####################################################################!

!<keyword category="Gas Phase" required="false">
!  <description>
!    Specified constant gas density [g/cm^3 in CGS]. An equation of
!    state -the ideal gas law by default- is used to calculate the gas
!    density if this parameter is undefined. The value may be set to
!    zero to make the drag zero and to simulate granular flow in a
!    vacuum. For this case, users may turn off solving for gas momentum
!    equations to accelerate convergence.
!  </description>
      RO_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Specified constant gas viscosity [g/(cm.s) in CGS].
!  </description>
      MU_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>
!    Average molecular weight of gas [(g/mol) in CGS]. Used in
!    calculating the gas density for non-reacting flows when the gas
!    composition is not defined.
!  </description>
      MW_AVG = UNDEFINED
!</keyword>



!#####################################################################!
!                            Solids Phase                             !
!#####################################################################!

!<keyword category="Solids Phase" required="false">
!  <description>
!    Defines the model used for the solids phase. For TFM/DEM
!    hybrid simulations, first define all TFM solids, then
!    define the DEM solids phases.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <valid value='TFM' note='Two-fluid Model (continuum)' />
!  <valid value='DEM' note='Discrete Element Model' />
!  <valid value='PIC' note='Multiphase-Particle in Cell' />
      SOLIDS_MODEL(:DIM_M) = 'TFM'
!</keyword>

!<keyword category="Solids Phase" required="false"
!  tfm="true" dem="true" pic="true">
!  <description>Number of solids phases.</description>
      MMAX = 1
!</keyword>

!<keyword category="Solids Phase" required="false"
!  tfm="true" dem="true" pic="true">
!  <description>
!    Initial particle diameters [cm in CGS].
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      D_P0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false"
!  tfm="true" dem="true" pic="true">
!  <description>
!    Specified constant solids density [g/cm^3 in CGS]. Reacting flows
!    may use variable solids density by leaving this parameter
!    undefined and specifying X_S0 and RO_XS0 as well as the index
!    of the inert species.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      RO_S0(:DIM_M) = UNDEFINED
!</keyword>


!#####################################################################!
!                   Initial Conditions Section                        !
!#####################################################################!


      DO LC = 1, DIMENSION_IC

!<keyword category="Initial Condition" required="false">
!  <description>X coordinate of the west face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>X coordinate of the east face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Y coordinate of the south face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Y coordinate of the north face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Z coordinate of the bottom face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Z coordinate of the top face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>I index of the west-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>I index of the east-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>J index of the south-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>J index of the north-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>K index of the bottom-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>K index of the top-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Type of initial condition. Mainly used in restart runs to overwrite
!    values read from the .RES file by specifying it as _PATCH_. The
!    user needs to be careful when using the _PATCH_ option, since the
!    values from the .RES file are overwritten and no error checking is
!    done for the patched values.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial void fraction in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial gas pressure in the IC region. If this quantity is not
!    specified, MFIX will set up a hydrostatic pressure profile,
!    which varies only in the y-direction.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial bulk density (rop_s = ro_s x ep_s) of solids phase-m in the
!    IC region. Users need to specify this IC only for polydisperse flow
!    (MMAX > 1). Users must make sure that summation of ( IC_ROP_s(ic,m)
!    / RO_s(m) ) over all solids phases is equal to ( 1.0 - IC_EP_g(ic)).
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial solids volume fraction of solids phase-m in the IC region.
!    This may be specified in place of IC_ROP_s.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_EP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial gas phase temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids phase-m temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Initial Condition" required="false">
!  <description>Flag for inflating initial lattice distribution
! to the entire IC region. </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
          IC_DES_FIT_TO_REGION(LC) = .FALSE.
!</keyword>

      ENDDO




!#####################################################################!
!                        Boundary Conditions                          !
!#####################################################################!
      DO LC = 1, DIMENSION_BC


!<keyword category="Boundary Condition" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y coordinate of the south face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y coordinate of the north face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>I index of the west-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>I index of the east-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>J index of the south-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>J index of the north-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>K index of the bottom-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>K index of the top-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Type of boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!
!  <valid value='DUMMY'
!    note='The specified boundary condition is ignored. This is
!      useful for turning off some boundary conditions without having
!      to delete them from the file.' />
!
!  <valid value='MASS_INFLOW' alias='MI'
!    note='Mass inflow rates for gas and solids phases are
!      specified at the boundary.'/>
!
!  <valid value='MASS_OUTFLOW' alias='MO'
!    note='The specified values of gas and solids mass outflow
!      rates at the boundary are maintained, approximately. This
!      condition should be used sparingly for minor outflows, when
!      the bulk of the outflow is occurring through other constant
!      pressure outflow boundaries.' />
!
!  <valid value='P_INFLOW' alias='PI'
!    note='Inflow from a boundary at a specified constant
!      pressure. To specify as the west, south, or bottom end of
!      the computational region, add a layer of wall cells to the
!      west, south, or bottom of the PI cells. Users need to specify
!      all scalar quantities and velocity components. The specified
!      values of fluid and solids velocities are only used initially
!      as MFIX computes these values at this inlet boundary.' />
!
!  <valid value='P_OUTFLOW' alias='PO'
!    note='Outflow to a boundary at a specified constant pressure.
!      To specify as the west, south, or bottom end of the computational
!      region, add a layer of wall cells to the west, south, or bottom of
!      the PO cells.' />
!
!  <valid value='FREE_SLIP_WALL' alias='FSW'
!    note='Velocity gradients at the wall vanish./>
!
!  <valid value='NO_SLIP_WALL' alias='NSW'
!    note='All components of the velocity vanish at the wall./>
!
!  <valid value='PAR_SLIP_WALL' alias='PSW'
!    note='Partial slip at the wall implemented as
!      dv/dn + hw (v - vw) = 0, where n is the normal pointing from the
!      fluid into the wall. The coefficients hw and vw should be
!      specified. For free slip set hw = 0. For no slip leave hw
!      undefined (hw=+inf) and set vw = 0. To set hw = +inf, leave it
!      unspecified. />
         BC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_HW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_UW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_WW_G(LC) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_HW_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified gas phase wall temperature, Tw_g, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_TW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas phase heat flux, C, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_C_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_HW_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified wall gas species mass fraction, Xw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Gas phase Xw for mass transfer.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_XW_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas species mass flux, C, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_C_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Solid phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_s)/dn + Hw (X_s - Xw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_HW_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Void fraction at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas pressure at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Bulk density of solids phase at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids volume fraction at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_EP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase-m temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of gas species at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of solids species at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VOLFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VOLFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval at the beginning when the normal
!    velocity at the boundary is equal to BC_Jet_g0. When restarting,
!    run this value and BC_Jet_g0 should be specified such that the
!    transient jet continues correctly. MFIX does not store the jet
!    conditions. For MASS_OUTFLOW boundary conditions, BC_DT_0 is
!    the time period to average and print the outflow rates. The
!    adjustment of velocities to get a specified mass or volumetric
!    flow rate is based on the average outflow rate.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the initial interval BC_DT_0.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_G0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval when normal velocity is equal to BC_Jet_gh.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_H(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the interval BC_DT_h.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_GH(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval when normal velocity is equal to BC_JET_gL.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_L(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the interval BC_DT_L.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_GL(LC) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Magnitude of gas velocity in a specified boundary region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VELMAG_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Magnitude of gas velocity in a specified boundary region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VELMAG_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Flag to specify the constant number
! of computational particles per cell for the PIC solids inflow BC.
!Statistical weight of parcels will be calculated by the code.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <conflict keyword="BC_PIC_CONST_STATWT" value="DEFINED"/>
!  <dependent keyword="SOLIDS_MODEL" value="PIC"/>
          BC_PIC_MI_CONST_NPC(LC, :DIM_M) = 0
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Flag to specify the constant statistical
! weight for inflowing computational particles/parcels. Actual number of
! parcels will be automatically computed. </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <conflict keyword="IC_PIC_CONST_NPC" value="DEFINED"/>
          BC_PIC_MI_CONST_STATWT(LC, :DIM_M) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Flag to make the PO BC invisible to discrete solids.
! Set this flag to.FALSE.to remove this BC for discrete solids. </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
         BC_PO_APPLY_TO_DES(LC) = .TRUE.
!</keyword>


         BC_ROP_G(LC) = UNDEFINED
      ENDDO




!#####################################################################!
!                     Point Source Mass Inlets                        !
!#####################################################################!
      DO LC = 1, DIMENSION_PS

!<keyword category="Point Source" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y coordinate of the south face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y coordinate of the north face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>I index of the west-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>I index of the east-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>J index of the south-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>J index of the north-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>K index of the bottom-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>K index of the top-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_K_T(LC) = UNDEFINED_I
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>X-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas mass flow rate through the point source.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming gas.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas phase incoming species n mass fraction.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         PS_X_G(LC,:DIM_N_g) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>X-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Solids mass flow rate through the point source.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming solids.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Solids phase incoming species n mass fraction.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         PS_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

      ENDDO


!#####################################################################!
!                          Output Control                             !
!#####################################################################!

!<keyword category="Output Control" required="true">
!  <description>
!    Interval at which restart (.res) file is updated.
!  </description>
      RES_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Interval at which a backup copy of the restart file is created.
!  </description>
      RES_BACKUP_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    The number of backup restart files to retain.
!  </description>
      RES_BACKUPS = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Interval at which .SPX files are updated.
!    o SP1: void fraction (EP_G)
!    o SP3: Gas velocity (U_G, V_G, W_G)
!  </description>
!  <arg index="1" id="SP Value" min="1" max="N_SPX"/>
      SPX_DT(:N_SPX) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description> Interval at which standard output (.OUT) file is updated.
!    Only run configuration information is written if left undefined. Otherwise
!    all field variables for the entire domain are written in ASCII
!    format to the .OUT file at OUT_DT intervals.
!  </description>
      OUT_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Number of time steps between .LOG file updates.</description>
      NLOG = 25
!</keyword>

!<keyword category="Output Control" required="false">
!  <description> Display the residuals on the screen and provide
!    messages about convergence on the screen and in the .LOG file.
!  </description>
      FULL_LOG = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Specifies the residuals to display. </description>
!  <arg index="1" id="Residual Index" max="8" min="1"/>
!  <valid value="P0" note="Gas pressure"/>
!  <valid value="PM" note="Solids phase M pressure"/>
!  <valid value="R0" note="Gas density"/>
!  <valid value="RM" note="Solids phase M density"/>
!  <valid value="U0" note="Gas phase U-velocity"/>
!  <valid value="V0" note="Gas phase V-velocity"/>
!  <valid value="W0" note="Gas phase W-velocity"/>
!  <valid value="UM" note="Solids phase M U-velocity"/>
!  <valid value="VM" note="Solids phase M V-velocity"/>
!  <valid value="WM" note="Solids phase M W-velocity"/>
!  <valid value="T0" note="Gas temperature"/>
!  <valid value="TM" note="Solids phase M temperature"/>
!  <valid value="X0NN" note="Gas phase species NN mass fraction"/>
!  <valid value="XMNN" note="Solids phase M species NN mass fraction"/>
!  <valid value="K0" note="K-Epsilon model residuals"/>
      RESID_STRING(:8) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Display residuals by equation.  </description>
      GROUP_RESID = .FALSE.
!</keyword>


!<keyword category="Output Control" required="false">
!  <description>
!    Provide detailed logging of negative density errors.
!  </description>
!  <valid value=".FALSE." note="Do not log negative density errors."/>
!  <valid value=".TRUE." note="Log negative density errors."/>
      REPORT_NEG_DENSITY = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Frequency to perform an overall species mass balance. Leaving
!    undefined suppresses the mass balance calculations which can
!    slightly extend run time.
!  </description>
      REPORT_MASS_BALANCE_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Use distributed IO :: Each MPI process generates RES/SPx files.
!  </description>
      bDist_IO = .FALSE.
!</keyword>



!#####################################################################!
!                           UDF  Control                              !
!#####################################################################!

!<keyword category="UDF Control" required="false">
!  <description>
!    Flag to enable user-defined subroutines: USR0, USR1, USR2, USR3,
!    USR0_DES, USR1_DES, USR2_DES, USR3_DES
!  </description>
!  <valid value=".TRUE." note="Call user-defined subroutines."/>
!  <valid value=".FALSE." note="Do NOT call user-defined subroutines."/>
      CALL_USR = .FALSE.
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>User defined constants.</description>
      C(:DIMENSION_C) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Name of user-defined constant. (20 character max)</description>
      C_NAME(:DIMENSION_C) = '....................'
!</keyword>

      DO LC=1, DIMENSION_USR
!<keyword category="UDF Control" required="false">
!  <description>
!    Intervals at which subroutine write_usr1 is called.
!  </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_DT(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: x coordinate of the west face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: x coordinate of the east face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: y coordinate of the south face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: y coordinate of the north face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: z coordinate of the top face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: i index of the west-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: i index of the east-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: j index of the south-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: j index of the north-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: k index of the bottom-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: k index of the top-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: Type of user-defined output: Binary of ASCII.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook:
!    Variables to be written in the user-defined output files.
!  </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_VAR(LC) = UNDEFINED_C
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook:
!    Format for writing user-defined (ASCII) output file.
!  </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_FORMAT(LC) = UNDEFINED_C
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Udf Hook: File extension for the user-defined output.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_EXT(LC) = UNDEFINED_C
!</keyword>
      ENDDO


!#####################################################################!
!                        Chemical Reactions                           !
!#####################################################################!



!#####################################################################!
!                    Parallelization Control                          !
!#####################################################################!


!<keyword category="Parallelization Control" required="false">
!  <description>Number of grid blocks in x-direction.</description>
      NODESI = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Number of grid blocks in y-direction.</description>
      NODESJ = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Number of grid blocks in z-direction.</description>
      NODESK = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Print out additional statistics for parallel runs</description>
      solver_statistics = .FALSE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Group residuals to reduce global collectives.</description>
      DEBUG_RESID = .TRUE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>All ranks write error messages.</description>
      ENABLE_DMP_LOG = .FALSE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Print the index layout for debugging.</description>
      DBGPRN_LAYOUT = .FALSE.
!</keyword>


!#####################################################################!
!                       Batch Queue Environment                       !
!#####################################################################!


!<keyword category="Batch Queue Environment" required="false">
!  <description>
!    Enables controlled termination feature when running under batch
!    queue system to force MFIX to cleanly terminate before the end
!    of wall clock allocated in the batch session.
!  </description>
      CHK_BATCHQ_END = .FALSE.
!</keyword>

!<keyword category="Batch Queue Environment" required="false">
!  <description>Total wall-clock duration of the job, in seconds.</description>
      BATCH_WALLCLOCK = 9000.0    ! set to 2.5 hrs for jaguarcnl w/ nproc<=512
!</keyword>

!<keyword category="Batch Queue Environment" required="false">
!  <description>
!    Buffer time specified to allow MFIX to write out the files and
!    cleanly terminate before queue wall clock time limit is reached
!    such that (BATCH_WALLCLOCK-TERM_BUFFER) is less than then batch
!    queue wall clock time limit, in seconds.
!  </description>
      TERM_BUFFER = 180.0         ! set to 3 minutes prior to end of job
!</keyword>




! ---------------------------------- questionable namelist entries below








!<keyword category="category name" required="false">
!  <description>Variable which triggers an automatic restart.</description>
      AUTOMATIC_RESTART = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>AUTO_RESTART counter.</description>
      ITER_RESTART = 1
!</keyword>




      U_G0 = UNDEFINED
      V_G0 = UNDEFINED
      W_G0 = UNDEFINED
      U_S0(:DIM_M) = UNDEFINED
      V_S0(:DIM_M) = UNDEFINED
      W_S0(:DIM_M) = UNDEFINED



      CALL DES_INIT_NAMELIST


      CALL USR_INIT_NAMELIST

      CALL CARTESIAN_GRID_INIT_NAMELIST

      RETURN
      END SUBROUTINE INIT_NAMELIST
