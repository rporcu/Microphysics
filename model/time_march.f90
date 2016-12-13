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

      subroutine time_march(u_g,v_g,w_g,u_go,v_go,w_go,&
         p_g,p_go,pp_g,ep_g,ep_go,&
         ro_g,ro_go,rop_g,rop_go,&
         rop_ge,rop_gn,rop_gt,d_e,d_n,d_t,&
         tau_u_g,tau_v_g,tau_w_g,&
         flux_ge,flux_gn,flux_gt,trd_g,lambda_g,mu_g,&
         f_gds, drag_am, drag_bm,flag, &
         pijk, dg_pijk, dg_pijkprv, iglobal_id, particle_state, particle_phase, &
         des_radius, ro_sol, pvol, pmass, omoi, neighbor_index, neighbor_index_old, &
         ppos, des_pos_new, des_vel_new, des_usr_var, omega_new, des_acc_old, rot_acc_old, fc, tow, wall_collision_pft)

      USE check_batch_queue_end_module, only: check_batch_queue_end
      USE check_data_30_module, only: check_data_30
      USE compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
      USE compar, only: myPE
      USE discretelement, only: des_continuum_coupled, max_pip
      USE error_manager, only: err_msg, flush_err_msg
      USE fld_const, only: mu_g0
      USE leqsol, only: SOLVER_STATISTICS, REPORT_SOLVER_STATS
      USE param1, only: undefined, small_number, zero
      USE run, only: automatic_restart, auto_restart, chk_batchq_end, call_usr
      USE run, only: automatic_restart, auto_restart, chk_batchq_end, dem_solids
      USE run, only: dem_solids
      USE run, only: time, tstop, nstep, dt, dt_min, dt_prev, use_dt_prev, units
      USE time_cpu, only: cpu_io

      USE toleranc , only: max_allowed_vel, max_inlet_vel_fac, max_inlet_vel

      ! Use function MAX_VEL_INLET to compute max. velocity at inlet
      USE utilities, ONLY: MAX_VEL_INLET

      use output   , only: RES_DT
      use adjust_dt, only: adjustdt

      use exit_mod      , only: mfix_exit
      use iterate_module, only: iterate

      use des_time_march_module, only: des_time_march
      use calc_coeff_module    , only: calc_coeff, calc_coeff_all, calc_trd_and_tau
      use set_bc1_module, only: set_bc1
      use output_manager_module, only: init_output_vars, output_manager

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
      DOUBLE PRECISION, INTENT(INOUT) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) :: wall_collision_pft
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: pvol, pmass, des_radius, ro_sol, omoi
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: des_acc_old, rot_acc_old, fc, tow
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: des_vel_new, des_pos_new, ppos, omega_new, des_usr_var
      INTEGER(KIND=1), DIMENSION(:), INTENT(INOUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(INOUT) :: NEIGHBOR_INDEX, NEIGHBOR_INDEX_OLD
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijk, iglobal_id
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijkprv
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Flag to indicate one pass through iterate for steady
! state conditions.
      LOGICAL :: FINISH

! Error index
      INTEGER :: IER
! Number of iterations
      INTEGER :: NIT, NIT_TOTAL
! used for activating check_data_30
      INTEGER :: NCHECK, DNCHECK

! Flag to save results and cleanly exit.
      LOGICAL :: EXIT_SIGNAL = .FALSE.

!-----------------------------------------------

      FINISH  = .FALSE.
      NCHECK  = NSTEP
      DNCHECK = 1
      CPU_IO  = ZERO
      NIT_TOTAL = 0
      IER = 0

      CALL INIT_OUTPUT_VARS

! Parse residual strings
      CALL PARSE_RESID_STRING ()

! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0

! Calculate all the coefficients once before entering the time loop
      CALL CALC_COEFF(flag, 2, ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g, &
         f_gds, drag_am, drag_bm, pijk(1:MAX_PIP,:), particle_phase(1:MAX_PIP), particle_state(1:MAX_PIP), &
         pvol(1:MAX_PIP), des_pos_new(1:MAX_PIP, :), des_vel_new(1:MAX_PIP, :), des_radius(1:MAX_PIP))
      IF(MU_g0 == UNDEFINED) CALL CALC_MU_G(lambda_g,mu_g,mu_g0)

! Remove undefined values at wall cells for scalars
      where(rop_g == undefined) rop_g = 0.0

! Initialize d's and e's to zero
      D_E = 0.0d0
      D_N = 0.0d0
      D_T = 0.0d0

! The TIME loop begins here.............................................
 100  CONTINUE


! Terminate MFIX normally before batch queue terminates.
      IF (CHK_BATCHQ_END) CALL CHECK_BATCH_QUEUE_END(EXIT_SIGNAL)

      IF (CALL_USR) CALL USR1

! Set wall boundary conditions and transient flow b.c.'s
      CALL SET_BC1(p_g, ep_g, ro_g, rop_g, u_g, v_g,w_g, &
         flux_ge, flux_gn, flux_gt, flag)

      CALL OUTPUT_MANAGER(ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g, &
               iglobal_id, particle_state, des_radius, ro_sol, des_pos_new, des_vel_new, des_usr_var, omega_new, &
         EXIT_SIGNAL, FINISH)

      IF (DT == UNDEFINED) THEN
         IF (FINISH) THEN
            RETURN
         ELSE
            FINISH = .TRUE.
         ENDIF

! Mechanism to terminate MFIX normally before batch queue terminates.
      ELSEIF (TIME + 0.1d0*DT >= TSTOP .OR. EXIT_SIGNAL) THEN
         IF(SOLVER_STATISTICS) &
            CALL REPORT_SOLVER_STATS(NIT_TOTAL, NSTEP)
         RETURN
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
         particle_state, pmass, pvol, des_pos_new, des_vel_new, des_radius, des_usr_var, flag)

! Calculate the stress tensor trace and cross terms for all phases.
      CALL CALC_TRD_AND_TAU(tau_u_g,tau_v_g,tau_w_g,trd_g,&
         ep_g,u_g,v_g,w_g,lambda_g,mu_g, flag)

! Check rates and sums of mass fractions every NLOG time steps
      IF (NSTEP == NCHECK) THEN
         IF (DNCHECK < 256) DNCHECK = DNCHECK*2
         NCHECK = NCHECK + DNCHECK
         CALL CHECK_DATA_30(lambda_g,mu_g,flag)
      ENDIF

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
                   flux_ge,flux_gn,flux_gt,mu_g,f_gds, drag_am, drag_bm,&
                   tau_u_g,tau_v_g,tau_w_g,&
                   pijk, particle_phase, particle_state, pvol, des_radius, des_pos_new, des_vel_new, flag, IER, NIT)

      DO WHILE (ADJUSTDT(ep_g, ep_go, p_g, p_go, ro_g, ro_go, flag, rop_g, &
         rop_go, U_g,  U_go, V_g, V_go,  W_g,  W_go, mu_g, f_gds, &
         drag_am, drag_bm, pijk, particle_phase, iglobal_id, &
         particle_state, pmass, pvol, des_radius, des_pos_new, des_vel_new, des_usr_var, IER, NIT))

         call iterate(u_g,v_g,w_g,u_go,v_go,w_go,p_g,pp_g,ep_g,ro_g,rop_g,rop_go,&
                      rop_ge,rop_gn,rop_gt,d_e,d_n,d_t,&
                      flux_ge,flux_gn,flux_gt,mu_g,f_gds, drag_am, drag_bm,&
                      tau_u_g,tau_v_g,tau_w_g,&
                      pijk, particle_phase, particle_state, pvol, des_radius, des_pos_new, des_vel_new, flag, IER, NIT)
      ENDDO

      IF(DT < DT_MIN) THEN
         AUTO_RESTART = .FALSE.
         IF(AUTO_RESTART) THEN
            IF(TIME <= RES_DT) THEN

 1000 FORMAT('Automatic restart not possible as Total Time < RES_DT')

               WRITE(ERR_MSG,1000)
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
               AUTOMATIC_RESTART = .TRUE.
            ENDIF
            RETURN

         ELSE

 1100 FORMAT('DT < DT_MIN.  Recovery not possible!')

            WRITE(ERR_MSG,1100)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            CALL MFIX_EXIT(MyPE)
         ENDIF
      ENDIF



! Other solids model implementations
      IF(DEM_SOLIDS) THEN
         call des_time_march(ep_g, p_g, u_g, v_g, w_g, ro_g, rop_g, mu_g, &
            pijk, dg_pijk, dg_pijkprv, iglobal_id, particle_state, particle_phase, &
            neighbor_index, neighbor_index_old, des_radius, ro_sol, pvol, pmass, omoi, des_usr_var, &
            ppos, des_pos_new, des_vel_new, omega_new, des_acc_old, rot_acc_old, fc, tow, wall_collision_pft, flag)
         IF(.NOT.DES_CONTINUUM_COUPLED) RETURN
      ENDIF

      IF (DT /= UNDEFINED) THEN
         IF(USE_DT_PREV) THEN
            TIME = TIME + DT_PREV
         ELSE
            TIME = TIME + DT
         ENDIF
         USE_DT_PREV = .FALSE.
         NSTEP = NSTEP + 1
      ENDIF

      NIT_TOTAL = NIT_TOTAL+NIT


      FLUSH (6)
! The TIME loop ends here....................................................
      GOTO 100

      IF(SOLVER_STATISTICS) CALL REPORT_SOLVER_STATS(NIT_TOTAL, NSTEP)

      END SUBROUTINE TIME_MARCH

end module time_march_module
