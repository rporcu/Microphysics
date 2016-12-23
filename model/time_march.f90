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
                p_g, p_go, pp_g, ep_g, ep_go, &
                ro_g, ro_go, rop_g, rop_go, &
                rop_ge, rop_gn, rop_gt, &
                d_e, d_n, d_t, &
                tau_u_g ,tau_v_g, tau_w_g,&
                flux_ge, flux_gn, flux_gt, &
                trD_g, lambda_g, mu_g, &
                f_gds, A_m, b_m, &
                drag_am, drag_bm, pinc, &
                flag, vol_surr,  &
                pijk,   iglobal_id, &
                particle_state, particle_phase, des_radius, ro_sol, pvol, pmass, &
                omoi, ppos, des_pos_new, des_vel_new, des_usr_var, omega_new, des_acc_old,&
                rot_acc_old, drag_fc, fc, tow, wall_collision_pft)

      USE check_batch_queue_end_module, only: check_batch_queue_end
      USE compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
      USE compar, only: myPE
      USE discretelement, only: des_continuum_coupled, max_pip
      USE error_manager, only: err_msg, flush_err_msg
      USE param1, only: undefined, small_number, zero
      USE run, only: call_usr, dem_solids
      USE run, only: dem_solids
      USE run, only: time, tstop, nstep, dt, dt_min, dt_prev, use_dt_prev, units
      USE toleranc , only: max_allowed_vel, max_inlet_vel_fac, max_inlet_vel

      ! Use function MAX_VEL_INLET to compute max. velocity at inlet
      USE utilities, ONLY: MAX_VEL_INLET

      use adjust_dt, only: adjustdt

      use exit_mod      , only: mfix_exit
      use iterate_module, only: iterate

      use des_time_march_module, only: des_time_march
      use calc_coeff_module    , only: calc_coeff, calc_coeff_all, calc_trd_and_tau
      use set_bc1_module, only: set_bc1
      use run, only: chk_batchq_end
      USE check_batch_queue_end_module, only: check_batch_queue_end
      USE discretelement, only: des_continuum_coupled
      use des_time_march_module, only: des_time_march
      use output_manager_module, only: output_manager

! HACK: defer to global arrays
      use discretelement, only: NEIGHBOR_INDEX, NEIGHBOR_INDEX_OLD

      implicit none

      DOUBLE PRECISION, INTENT(INOUT) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: u_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: v_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: w_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: p_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ep_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: pp_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: d_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: trd_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: lambda_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: b_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      DOUBLE PRECISION, INTENT(INOUT) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
      DOUBLE PRECISION, INTENT(INOUT) :: vol_surr&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)


      integer         , INTENT(INOUT) :: pinc&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      integer, intent(inout) :: pijk(max_pip,3)
      integer, intent(inout) :: iglobal_id(max_pip)
      integer, intent(inout) :: particle_state(max_pip)
      integer, intent(inout) :: particle_phase(max_pip)

      double precision, intent(inout) :: des_radius(max_pip)
      double precision, intent(inout) :: ro_sol(max_pip)
      double precision, intent(inout) :: pvol(max_pip)
      double precision, intent(inout) :: pmass(max_pip)
      double precision, intent(inout) :: omoi(max_pip)

      double precision, intent(inout) :: ppos(max_pip,3)
      double precision, intent(inout) :: des_pos_new(max_pip,3)
      double precision, intent(inout) :: des_vel_new(max_pip,3)
      double precision, intent(inout) :: des_usr_var(max_pip,1)
      double precision, intent(inout) :: omega_new(max_pip,3)

      double precision, intent(inout) :: des_acc_old(max_pip,3)
      double precision, intent(inout) :: rot_acc_old(max_pip,3)
      double precision, intent(inout) :: drag_fc(max_pip,3)
      double precision, intent(inout) :: fc(max_pip,3)
      double precision, intent(inout) :: tow(max_pip,3)

      double precision, intent(inout) :: wall_collision_pft(3,8,max_pip)

!      INTEGER, DIMENSION(:), INTENT(INOUT) :: NEIGHBOR_INDEX, NEIGHBOR_INDEX_OLD

!-----------------------------------------------
! Local variables
!-----------------------------------------------

! Error index
      INTEGER :: IER
! Number of iterations
      INTEGER :: NIT
! Flag to indicate one pass through iterate for steady
! state conditions.
      LOGICAL :: FINISH
! Flag to save results and cleanly exit.
      LOGICAL :: EXIT_SIGNAL = .FALSE.

!-----------------------------------------------

      FINISH  = .FALSE.
      IER = 0

! Time march
      do
! Terminate MFIX normally before batch queue terminates.
         IF (CHK_BATCHQ_END) CALL CHECK_BATCH_QUEUE_END(EXIT_SIGNAL)

         IF (CALL_USR) CALL USR1

! Set wall boundary conditions and transient flow b.c.'s
         CALL SET_BC1(p_g, ep_g, ro_g, rop_g, u_g, v_g,w_g, &
                       flux_ge, flux_gn, flux_gt, flag)

         CALL OUTPUT_MANAGER(ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g, &
            iglobal_id, particle_state, des_radius, ro_sol, des_pos_new, &
            des_vel_new, des_usr_var, omega_new, EXIT_SIGNAL, FINISH)

         IF (DT == UNDEFINED) THEN
            IF (FINISH) THEN
               return
            ELSE
               FINISH = .TRUE.
            ENDIF

! Mechanism to terminate MFIX normally before batch queue terminates.
         ELSEIF (TIME + 0.1d0*DT >= TSTOP .OR. EXIT_SIGNAL) THEN
            return
         ENDIF

! Update previous-time-step values of field variables
         ep_go = ep_g
         p_go =  p_g
         ro_go = ro_g
         rop_go = rop_g
         U_go =  U_g
         V_go =  V_g
         W_go =  W_g

! Calculate coefficients
         CALL CALC_COEFF_ALL (ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g,&
            f_gds, drag_am, drag_bm, pijk, particle_phase, iglobal_id, &
            particle_state, pmass, pvol, des_pos_new, des_vel_new, des_radius, &
            des_usr_var, flag, vol_surr, pinc)

! Calculate the stress tensor trace and cross terms for all phases.
         CALL CALC_TRD_AND_TAU(tau_u_g,tau_v_g,tau_w_g,trd_g,&
            ep_g,u_g,v_g,w_g,lambda_g,mu_g, flag)

! Check for maximum velocity at inlet to avoid convergence problems
         MAX_INLET_VEL = 100.0d0*MAX_VEL_INLET(u_g,v_g,w_g)
! if no inlet velocity is specified, use an upper limit defined in
! toleranc_mod.f
         IF(MAX_INLET_VEL <= SMALL_NUMBER) THEN
            MAX_INLET_VEL = MAX_ALLOWED_VEL
            IF (UNITS == 'SI') MAX_INLET_VEL = 1D-2 * MAX_ALLOWED_VEL
         ENDIF
! Scale the value using a user defined scale factor
         MAX_INLET_VEL = MAX_INLET_VEL * MAX_INLET_VEL_FAC

! Advance the solution in time by iteratively solving the equations
         call iterate(u_g,v_g,w_g,u_go,v_go,w_go,p_g,pp_g,ep_g,ro_g,rop_g,rop_go,&
            rop_ge,rop_gn,rop_gt,d_e,d_n,d_t,&
            flux_ge,flux_gn,flux_gt,mu_g,f_gds, &
            A_m, b_m, drag_am, drag_bm,&
            tau_u_g,tau_v_g,tau_w_g,&
            pijk, particle_phase, particle_state, pvol, &
            des_radius,  des_pos_new, des_vel_new, flag, pinc, IER, NIT)

         DO WHILE (ADJUSTDT(ep_g, ep_go, p_g, p_go, ro_g, ro_go, flag, vol_surr, &
            rop_g, rop_go, U_g,  U_go, V_g, V_go,  W_g,  W_go, mu_g, f_gds, &
            drag_am, drag_bm, pijk, particle_phase, iglobal_id, &
            particle_state, pinc, pmass, pvol, des_radius,  &
            des_pos_new, des_vel_new, des_usr_var, IER, NIT))

            call iterate(u_g,v_g,w_g,u_go,v_go,w_go,p_g,pp_g,ep_g,ro_g,rop_g,rop_go,&
               rop_ge,rop_gn,rop_gt,d_e,d_n,d_t,&
               flux_ge,flux_gn,flux_gt,mu_g,f_gds, &
               A_m, b_m, drag_am, drag_bm,&
               tau_u_g,tau_v_g,tau_w_g,&
               pijk, particle_phase, particle_state, pvol, &
               des_radius,  des_pos_new, des_vel_new, flag, pinc, IER, NIT)
         ENDDO

         IF(DT < DT_MIN) THEN
            WRITE(ERR_MSG,1100)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            CALL MFIX_EXIT(MyPE)
         ENDIF
1100     FORMAT('DT < DT_MIN.  Recovery not possible!')

! Other solids model implementations
         IF(DEM_SOLIDS) THEN
            call des_time_march(ep_g, p_g, u_g, v_g, w_g, ro_g, rop_g, mu_g, &
               pijk,   iglobal_id, particle_state, particle_phase, &
               neighbor_index, neighbor_index_old, des_radius,  ro_sol, pvol, pmass,&
               omoi, des_usr_var, ppos, des_pos_new, des_vel_new, omega_new, &
               des_acc_old, rot_acc_old, drag_fc, fc, tow, wall_collision_pft, &
               flag, vol_surr, pinc)
            IF(.NOT.DES_CONTINUUM_COUPLED) return
         ENDIF

         IF(DT /= UNDEFINED)THEN
            IF(USE_DT_PREV)THEN
               TIME = TIME + DT_PREV
            ELSE
               TIME = TIME + DT
            ENDIF
            USE_DT_PREV = .FALSE.
            NSTEP = NSTEP + 1
         ENDIF

      enddo


      END SUBROUTINE TIME_MARCH

end module time_march_module
