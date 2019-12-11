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

      use ic
      use deprecated_or_unknown_module, only: deprecated_or_unknown
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar
      use ic, only: ic_ep_g, ic_ep_s, ic_p_g, ic_t_g, ic_t_s, ic_x_w
      use ic, only: ic_u_g, ic_u_s, ic_v_g, ic_v_s, ic_w_g, ic_w_s
      use ic, only: ic_x_e, ic_y_n, ic_y_s, ic_z_b, ic_z_t

      use usr
      use utilities, only: blank_line, line_too_big, seek_comment
      use utilities, only: make_upper_case, replace_tab


      use param, only: zero, one
      use param, only: undefined, undefined_c
      use param, only: dim_usr

      implicit none




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




      end subroutine init_namelist
end module init_namelist_module
