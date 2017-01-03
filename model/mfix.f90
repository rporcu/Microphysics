!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine mfix1(time, dt, nstep, u_g, v_g, w_g, u_go, v_go, w_go, &
   p_g, p_go, pp_g, ep_g, ep_go, &
   ro_g, ro_go, rop_g, rop_go, &
   rop_ge, rop_gn, rop_gt, &
   d_e, d_n, d_t, &
   flux_ge, flux_gn, flux_gt, &
   trD_g, lambda_g, mu_g, &
   f_gds, drag_bm,  flag) &
   bind(C, name="mfix_main1")

!-----------------------------------------------
! Modules
!-----------------------------------------------

      use check_data_20_module, only: check_data_20
      use compar, only: myPE, istart3, iend3, jstart3, jend3, kstart3, kend3
      use corner_module, only: get_corner_cells
      use discretelement, only: max_pip
      use error_manager, only: err_msg, finl_err_msg, flush_err_msg, init_err_msg
      use exit_mod, only: mfix_exit
      use fld_const, only: ro_g0, mu_g0
      use funits , only: dmp_log, unit_log
      use geometry, only: dx, dy, dz, ayz, axy, axz, vol, flag_mod
      use set_domain_module, only: set_domain
      use machine, only: wall_time
      use param1 , only: zero, is_defined, is_undefined, undefined
      use run, only: call_usr, run_type, dem_solids
      use run, only: dt_min, dt_max
      use set_bc0_module, only: set_bc0
      use set_bc1_module, only: set_bc1
      use set_constprop_module, only: set_constprop
      use set_flags_module, only: set_flags1
      use set_fluidbed_p_module, only: set_fluidbed_p
      use set_ic_module, only: set_ic
      use set_ps_module, only: set_ps
      use set_ro_g_module, only: set_ro_g
      use write_out0_module, only: write_out0
      use write_out1_module, only: write_out1
      use write_out3_module, only: write_out3
      use zero_norm_vel_module, only: zero_norm_vel
      use output_manager_module, only: init_output_vars
      use calc_coeff_module, only: calc_coeff

      use discretelement, only: nonexistent, do_old

      IMPLICIT NONE

      integer         , intent(inout) :: nstep
      double precision, intent(inout) :: time, dt

      integer         , intent(inout) :: flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      double precision, intent(inout) :: ep_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ep_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: p_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: p_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ro_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ro_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: u_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: u_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: v_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: v_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: w_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: w_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: pp_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_e&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_t&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_n&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: mu_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: lambda_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: trD_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: flux_gE&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: flux_gN&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: flux_gT&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: rop_gE&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_gN&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_gT&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: f_gds&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: drag_bm&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,3)

!---------------------------------------------------------------------//
      flag_mod = flag

! Write the initial part of the standard output file
      CALL WRITE_OUT0(time, dt)

! Write the initial part of the special output file(s)
      CALL WRITE_USR0

      CALL INIT_ERR_MSG('MFIX')

! Set the flags for wall surfaces impermeable and identify flow
! boundaries using FLAG_E, FLAG_N, and FLAG_T
      CALL SET_FLAGS1(flag)
      flag_mod = flag

! Calculate cell volumes and face areas
      VOL = DX*DY*DZ
      AYZ = DY*DZ
      AXY = DX*DY
      AXZ = DX*DZ

! Find corner cells and set their face areas to zero
      CALL GET_CORNER_CELLS(flag)
      flag_mod = flag

! Set constant physical properties
      CALL set_constprop(ro_g, lambda_g, mu_g, flag(:,:,:,1))

! Set initial conditions
      CALL SET_IC(ep_g, p_g, u_g, v_g, w_g, flag)

! Set point sources.
      CALL SET_PS(flag)

! Set boundary conditions
      CALL ZERO_NORM_VEL(u_g,v_g,w_g, flag)
      CALL SET_BC0(p_g,ep_g,u_g,v_g,w_g,ro_g0,flag)

! Set the pressure field for a fluidized bed
      IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P(p_g, ep_g, flag)

! Initialize densities.
      IF (RUN_TYPE == 'NEW') CALL SET_RO_G(ro_g,rop_g,p_g,ep_g,flag)

! Remove undefined values at wall cells for scalars
      where(rop_g == undefined) rop_g = 0.0

! Initialize time dependent boundary conditions
      CALL SET_BC1(time, dt, p_g, ep_g, ro_g, rop_g, u_g, v_g, w_g, &
         flux_ge, flux_gn, flux_gt, flag)

! Check the field variable data and report errors.
      CALL CHECK_DATA_20(ep_g,p_g,ro_g,rop_g,u_g,v_g,w_g,flag)

      IF(IS_UNDEFINED(mu_g0)) CALL CALC_MU_G(lambda_g,mu_g,mu_g0)

      CALL INIT_OUTPUT_VARS(time, dt)

! Parse residual strings
      CALL PARSE_RESID_STRING ()

      end subroutine mfix1

subroutine mfix2(time, dt, nstep, u_g, v_g, w_g, u_go, v_go, w_go, &
   p_g, p_go, pp_g, ep_g, ep_go, &
   ro_g, ro_go, rop_g, rop_go, &
   rop_ge, rop_gn, rop_gt, &
   d_e, d_n, d_t, &
   flux_ge, flux_gn, flux_gt, &
   trD_g, lambda_g, mu_g, &
   f_gds, drag_bm,  flag, &
   particle_state, particle_phase, des_radius, ro_sol, pvol, pmass, &
   omoi, des_pos_new, des_vel_new, des_usr_var, omega_new, des_acc_old,&
   rot_acc_old, drag_fc, fc, tow)&
   bind(C, name="mfix_main2")

!-----------------------------------------------
! Modules
!-----------------------------------------------

      use check_data_20_module, only: check_data_20
      use compar, only: myPE, istart3, iend3, jstart3, jend3, kstart3, kend3
      use corner_module, only: get_corner_cells
      use discretelement, only: max_pip
      use error_manager, only: err_msg, finl_err_msg, flush_err_msg, init_err_msg
      use exit_mod, only: mfix_exit
      use fld_const, only: ro_g0, mu_g0
      use funits , only: dmp_log, unit_log
      use geometry, only: dx, dy, dz, ayz, axy, axz, vol, flag_mod
      use set_domain_module, only: set_domain
      use machine, only: wall_time
      use make_arrays_des_module, only: make_arrays_des
      use param1 , only: zero, is_defined, is_undefined, undefined
      use run, only: call_usr, run_type, dem_solids
      use run, only: dt_min, dt_max
      use set_bc0_module, only: set_bc0
      use set_bc1_module, only: set_bc1
      use set_constprop_module, only: set_constprop
      use set_flags_module, only: set_flags1
      use set_fluidbed_p_module, only: set_fluidbed_p
      use set_ic_module, only: set_ic
      use set_ps_module, only: set_ps
      use set_ro_g_module, only: set_ro_g
      use write_out0_module, only: write_out0
      use write_out1_module, only: write_out1
      use write_out3_module, only: write_out3
      use zero_norm_vel_module, only: zero_norm_vel
      use output_manager_module, only: init_output_vars
      use calc_coeff_module, only: calc_coeff

      use discretelement, only: nonexistent, do_old

      IMPLICIT NONE

      integer         , intent(inout) :: nstep
      double precision, intent(inout) :: time, dt

      integer         , intent(inout) :: flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      double precision, intent(inout) :: ep_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ep_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: p_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: p_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ro_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ro_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: u_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: u_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: v_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: v_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: w_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: w_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: pp_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_e&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_t&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_n&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: mu_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: lambda_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: trD_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: flux_gE&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: flux_gN&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: flux_gT&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: rop_gE&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_gN&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_gT&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: f_gds&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: drag_bm&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,3)

      integer, intent(inout) :: particle_state(max_pip)
      integer, intent(inout) :: particle_phase(max_pip)

      double precision, intent(inout) :: des_radius(max_pip)
      double precision, intent(inout) :: ro_sol(max_pip)
      double precision, intent(inout) :: pvol(max_pip)
      double precision, intent(inout) :: pmass(max_pip)
      double precision, intent(inout) :: omoi(max_pip)

      double precision, intent(inout) :: des_pos_new(max_pip,3)
      double precision, intent(inout) :: des_vel_new(max_pip,3)
      double precision, intent(inout) :: des_usr_var(max_pip,1)
      double precision, intent(inout) :: omega_new(max_pip,3)

      double precision, intent(inout) :: des_acc_old(max_pip,3)
      double precision, intent(inout) :: rot_acc_old(max_pip,3)
      double precision, intent(inout) :: drag_fc(max_pip,3)
      double precision, intent(inout) :: fc(max_pip,3)
      double precision, intent(inout) :: tow(max_pip,3)

!---------------------------------------------------------------------//

! First instance of particle data ......................................
      IF(DEM_SOLIDS) CALL MAKE_ARRAYS_DES(ep_g, &
         flag, particle_state, particle_phase,  &
         des_radius,  ro_sol, pvol, pmass, omoi, &
         des_pos_new, des_vel_new, des_usr_var, omega_new, fc, tow)


! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0

! Calculate all the coefficients once before entering the time loop
      CALL CALC_COEFF(flag, 2, ro_g, p_g, ep_g, rop_g, u_g, v_g, &
         w_g, mu_g, f_gds, drag_bm, particle_phase, particle_state, &
         pvol, des_pos_new, des_vel_new, des_radius)

      CALL FINL_ERR_MSG

      end subroutine mfix2
