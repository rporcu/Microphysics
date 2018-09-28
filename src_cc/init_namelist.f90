MODULE INIT_NAMELIST_MODULE
   CONTAINS
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

      use bc
      use ic
      use drag, only: drag_c1, drag_d1
      use constant, only: gravity
      use deprecated_or_unknown_module, only: deprecated_or_unknown
      use des_init_namelist_module, only: des_init_namelist
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar
      use fld_const, only: mu_g0, mw_avg
      use fld_const, only: ro_g0
      use fld_const, only: ro_g0, mu_g0, mw_avg
      use ic, only: ic_ep_g, ic_ep_s, ic_p_g, ic_t_g, ic_t_s, ic_x_w
      use ic, only: ic_u_g, ic_u_s, ic_v_g, ic_v_s, ic_w_g, ic_w_s
      use ic, only: ic_x_e, ic_y_n, ic_y_s, ic_z_b, ic_z_t
      use run, only: full_log, nlog
      use output, only: usr_dt
      use output, only: usr_x_w, usr_x_e, usr_y_n, usr_y_s, usr_z_b, usr_z_t
      use ps, only: dim_ps
      use ps, only: ps_massflow_g
      use ps, only: ps_t_g, ps_u_g, ps_v_g, ps_w_g
      use ps, only: ps_x_e, ps_x_g, ps_y_n, ps_y_s, ps_z_b, ps_z_t,  ps_x_w
      use run, only: call_usr, description
      use run, only: dt_fac, dt_max, dt_min, run_name
      use drag, only: drag_type
      use scales, only: p_ref, p_scale
      use usr
      use utilities, only: blank_line, line_too_big, seek_comment
      use utilities, only: make_upper_case, replace_tab


      use param, only: zero, one
      use param, only: undefined, undefined_c
      use param, only: dim_usr

      implicit none

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

!<keyword category="Run Control" required="false">
!  <description>Maximum time step size.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      dt_max = ONE
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Minimum time step size.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      dt_min = 1.0D-6
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Factor for adjusting time step.
!    * The value must be less than or equal to 1.0.
!    * A value of 1.0 keeps the time step constant which may help overcome
!      initial non-convergence.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <range min="0.0" max="1" />
      dt_fac = 0.9D0
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
!  <valid value="useR_DRAG" note="Invoke user-defined drag law. (usr_drag.f)"/>
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
!  <description>Gravity vector. [0.0, 9.80, 0.0] m/s^2 </description>
      gravity(1) =  0.00000d0
      gravity(2) = -9.80665d0
      gravity(3) =  0.00000d0
!</keyword>

!#####################################################################!
!                      Geometry and Discretization                    !
!#####################################################################!

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across xlength when a cyclic boundary condition
!    with pressure drop is imposed in the x-direction.
!  </description>
      delp_x = zero
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across ylength when a cyclic boundary condition
!    with pressure drop is imposed in the y-direction.
!  </description>
      delp_y = zero
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across zlength when a cyclic boundary condition
!    with pressure drop is imposed in the z-direction.
!  </description>
      delp_z = zero
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
!                   Initial Conditions Section                        !
!#####################################################################!


!<keyword category="Initial Condition" required="false">
!  <description>X coordinate of the west face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_X_W(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>X coordinate of the east face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_X_E(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Y coordinate of the south face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_Y_S(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Y coordinate of the north face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_Y_N(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Z coordinate of the bottom face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_Z_B(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Z coordinate of the top face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_Z_T(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial void fraction in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_EP_G(:) = ONE
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial gas pressure in the IC region. If this quantity is not
!    specified, MFIX will set up a hydrostatic pressure profile,
!    which varies only in the y-direction.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_P_G(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial solids volume fraction of solids phase-m in the IC region.
!    This may be specified in place of IC_ROP_s.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      IC_EP_S(:,:) = ZERO
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial gas phase temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_T_G(:) = 293.15d0
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids phase-m temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      IC_T_S(:,:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_U_G(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      IC_U_S(:,:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_V_G(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      IC_V_S(:,:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_W_G(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      IC_W_S(:,:) = UNDEFINED
!</keyword>

      ic_dp_dist(:,:) = 'CONSTANT'
      ic_dp_mean(:,:) = undefined
      ic_dp_std(:,:) = undefined
      ic_dp_min(:,:) = undefined
      ic_dp_max(:,:) = undefined

      ! Particle density properties
      ic_ro_s_dist(:,:) = 'CONSTANT'
      ic_ro_s_mean(:,:) = undefined
      ic_ro_s_std(:,:) = undefined
      ic_ro_s_min(:,:) = undefined
      ic_ro_s_max(:,:) = undefined

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
!  <valid value='HCP'    note='Hexagonal Close Pack'/>
!  <valid value='RANDOM' note='Random Fill -- requires EPs < 30%'/>
      ic_pack_type = 'HCP'
!</keyword>

!#####################################################################!
!                        Boundary Conditions                          !
!#####################################################################!

!<keyword category="Boundary Condition" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_X_W(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_X_E(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y coordinate of the south face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_Y_S(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y coordinate of the north face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_Y_N(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_Z_B(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_Z_T(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_NORMAL(:,:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_CENTER(:,:) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Type of boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
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
      BC_TYPE(:) = UNDEFINED_C
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_HW_G(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_UW_G(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_VW_G(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_WW_G(:) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_HW_T_G(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified gas phase wall temperature, Tw_g, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_TW_G(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas phase heat flux, C, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_C_T_G(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
      BC_HW_X_G(:,:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified wall gas species mass fraction, Xw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Gas phase Xw for mass transfer.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
      BC_XW_G(:,:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas species mass flux, C, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
      BC_C_X_G(:,:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Void fraction at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_EP_G(:) = ONE
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas pressure at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_P_G(:) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Solids volume fraction at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      BC_EP_S(:,:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_T_G(:) = 293.15d0
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase-m temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      BC_T_S(:,:) = 293.15d0
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of gas species at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
      BC_X_G(:,:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of solids species at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
      BC_X_S(:,:,:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_U_G(:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      BC_U_S(:,:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_V_G(:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      BC_V_S(:,:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_W_G(:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      BC_W_S(:,:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_VOLFLOW_G(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      BC_VOLFLOW_S(:,:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_MASSFLOW_G(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      BC_MASSFLOW_S(:,:) = UNDEFINED
!</keyword>



!#####################################################################!
!                     Point Source Mass Inlets                        !
!#####################################################################!

!<keyword category="Point Source" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_X_W(:) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_X_E(:) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y coordinate of the south face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_Y_S(:) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y coordinate of the north face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_Y_N(:) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_Z_B(:) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_Z_T(:) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>X-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_U_G(:) = ZERO
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Y-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_V_G(:) = ZERO
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Z-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_W_G(:) = ZERO
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas mass flow rate through the point source.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_MASSFLOW_G(:) = ZERO
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming gas.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
      PS_T_G(:) = 293.15d0
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas phase incoming species n mass fraction.</description>
!  <arg index="1" id="PS" min="1" max="DIM_PS"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
      PS_X_G(:,:) = UNDEFINED
!</keyword>



!#####################################################################!
!                          Output Control                             !
!#####################################################################!

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
!  <description>
!    Intervals at which subroutine write_usr1 is called.
!  </description>
!  <arg index="1" id="USR" max="DIM_USR" min="1"/>
      usr_dt(1:DIM_USR) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="USR" min="1" max="DIM_USR"/>
      USR_X_W(1:DIM_USR) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="USR" min="1" max="DIM_USR"/>
      USR_X_E(1:DIM_USR) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Y coordinate of the south face or edge.</description>
!  <arg index="1" id="USR" min="1" max="DIM_USR"/>
      USR_Y_S(1:DIM_USR) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Y coordinate of the north face or edge.</description>
!  <arg index="1" id="USR" min="1" max="DIM_USR"/>
      USR_Y_N(1:DIM_USR) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="USR" min="1" max="DIM_USR"/>
      USR_Z_B(1:DIM_USR) = UNDEFINED
!</keyword>

!<keyword category="UDF Control" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="USR" min="1" max="DIM_USR"/>
      USR_Z_T(1:DIM_USR) = UNDEFINED
!</keyword>



      call des_init_namelist

      end subroutine init_namelist
end module init_namelist_module
