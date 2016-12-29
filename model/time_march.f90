module time_march_module

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: TIME_MARCH                                              !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Controlling module for time marching and finding the       !
!           solution of equations from TIME to TSTOP at intervals of   !
!           DT, updating the b.c.'s, and creating output.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

     subroutine time_march(u_g, v_g, w_g, u_go, v_go, w_go, &
        p_g, p_go, ep_g, ep_go, ro_g, ro_go, rop_g, rop_go)&
        bind(C, name="mfix_time_march")

      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
      use iso_c_binding, only: c_double, c_int

      implicit none

      real(c_double), intent(inout) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: u_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: v_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: w_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: p_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: ep_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: ro_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), intent(inout) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)


! Time march
!      do

!         IF (CALL_USR) CALL USR1
!
! Set wall boundary conditions and transient flow b.c.'s
!         CALL SET_BC1(p_g, ep_g, ro_g, rop_g, u_g, v_g,w_g, &
!            flux_ge, flux_gn, flux_gt, flag)
!
!         CALL OUTPUT_MANAGER(ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g, &
!             particle_state, des_radius, ro_sol, des_pos_new, &
!             des_vel_new, des_usr_var, omega_new, EXIT_SIGNAL, FINISH)


! Update previous-time-step values of field variables
         ep_go = ep_g
         p_go =  p_g
         ro_go = ro_g
         rop_go = rop_g
         U_go =  U_g
         V_go =  V_g
         W_go =  W_g

! ! Calculate coefficients
!          CALL CALC_COEFF_ALL(ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, &
!             mu_g, f_gds, drag_bm,  particle_phase,  &
!             particle_state, pvol, des_pos_new, des_vel_new, des_radius,  &
!             flag)

! ! Calculate the stress tensor trace and cross terms for all phases.
!          CALL CALC_TRD_AND_TAU(tau_u_g,tau_v_g,tau_w_g,trd_g,&
!             ep_g,u_g,v_g,w_g,lambda_g,mu_g, flag)

! ! Check for maximum velocity at inlet to avoid convergence problems
!          MAX_INLET_VEL = 100.0d0*MAX_VEL_INLET(u_g,v_g,w_g)
! ! if no inlet velocity is specified, use an upper limit defined in
! ! toleranc_mod.f
!          IF(MAX_INLET_VEL <= SMALL_NUMBER) THEN
!             MAX_INLET_VEL = MAX_ALLOWED_VEL
!             IF (UNITS == 'SI') MAX_INLET_VEL = 1D-2 * MAX_ALLOWED_VEL
!          ENDIF
! ! Scale the value using a user defined scale factor
!          MAX_INLET_VEL = MAX_INLET_VEL * MAX_INLET_VEL_FAC

! Advance the solution in time by iteratively solving the equations
!          DO
!             prev_dt = dt
!             call iterate(u_g, v_g, w_g, u_go, v_go, w_go,&
!                p_g, pp_g, ep_g, ro_g, rop_g, rop_go,&
!                rop_ge, rop_gn, rop_gt, d_e, d_n, d_t,&
!                flux_ge, flux_gn, flux_gt, mu_g,&
!                f_gds, A_m, b_m, drag_bm, &
!                tau_u_g, tau_v_g, tau_w_g, &
!                particle_phase, particle_state, pvol, &
!                des_radius,  des_pos_new, des_vel_new, flag,  IER, NIT)

!             IF (.NOT.ADJUSTDT(ep_g, ep_go, p_g, p_go, ro_g, ro_go, flag, &
!                rop_g, rop_go, U_g,  U_go, V_g, V_go,  W_g,  W_go, mu_g, &
!                f_gds, drag_bm,  particle_phase,  &
!                particle_state,  pvol, des_radius,  &
!                des_pos_new, des_vel_new, IER, NIT)) EXIT
!          ENDDO

!          IF(DT < DT_MIN) THEN
!             WRITE(ERR_MSG,1100)
!             CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
!             CALL MFIX_EXIT(MyPE)
!          ENDIF
! 1100     FORMAT('DT < DT_MIN.  Recovery not possible!')

! ! Other solids model implementations
!          IF(DEM_SOLIDS) THEN
!             call des_time_march(ep_g, p_g, u_g, v_g, w_g, ro_g, rop_g, mu_g, &
!                particle_state, particle_phase, &
!                des_radius,  ro_sol, pvol, pmass, omoi, des_usr_var, &
!                des_pos_new, des_vel_new, omega_new, des_acc_old, rot_acc_old, &
!                drag_fc, fc, tow, pairs, pair_count, flag)
!             IF(.NOT.DES_CONTINUUM_COUPLED) then
!                finish = 1
!                return
!             endif
!          ENDIF

! !      enddo

      END SUBROUTINE TIME_MARCH

end module time_march_module
