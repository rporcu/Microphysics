module iterate_module

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ITERATE                                                 C
!  Author: M. Syamlal                                 Date: 12-APR-96  C
!                                                                      C
!  Purpose: This module controls the iterations for solving equations  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ITERATE(u_g, v_g, w_g, u_go, v_go, w_go,&
                         p_g, pp_g, ep_g, ro_g, rop_g, rop_go,&
                         rop_ge, rop_gn, rop_gt, d_e, d_n, d_t,&
                         flux_ge, flux_gn, flux_gt, mu_g,&
                         f_gds, drag_am, drag_bm, &
                         tau_u_g, tau_v_g, tau_w_g, &
                         pijk, particle_phase, particle_state, pvol, des_radius, des_pos_new, des_vel_new, flag, IER, NIT)

      USE calc_mflux_module, only: calc_mflux
      USE check_convergence_module, only: check_convergence
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE compar   , only: myPE, PE_IO
      USE conv_rop_module, only: conv_rop
      USE correct_0_module, only: correct_0
      USE fld_const, only: ro_g0
      USE funits   , only: dmp_log, unit_log
      USE geometry , only: cyclic, cyclic_x, cyclic_y, cyclic_z, vol
      USE leqsol   , only: leq_adjust, max_nit
      USE output   , only: full_log, nlog
      USE param1   , only: small_number, undefined, zero, one
      USE physical_prop_module, only: physical_prop
      USE residual , only: resid, resid_p
      USE run, only: call_usr, dt, dt_prev, nstep, time, tstop
      USE time_cpu, only: cpuos, cpu_nlog, cpu0, time_nlog
      USE toleranc, only: norm_g
      USE vavg_mod, ONLY: vavg_g

      USE solve_pp_module, only: solve_pp_g
      USE solve_vel_star_module, only: solve_vel_star
      USE calc_coeff_module, only: calc_coeff
      USE set_bc1_module, only: set_bc1

      USE error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar
      USE tunit_module, only: get_tunit

      implicit none

      DOUBLE PRECISION, INTENT(INOUT) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: pp_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_go&
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
      DOUBLE PRECISION, INTENT(INOUT) :: flux_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pvol, des_radius
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new
      INTEGER(KIND=1), DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_phase
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Number of iterations
      INTEGER, INTENT(OUT) :: NIT
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! current cpu time used
      DOUBLE PRECISION :: CPU_NOW
! cpu time left
      DOUBLE PRECISION :: TLEFT
! flag indicating convergence status with MUSTIT = 0,1,2 implying
! complete convergence, non-covergence and divergence respectively
      INTEGER :: MUSTIT
! Normalization factor for gas & solids pressure residual
      DOUBLE PRECISION :: NORMg
! Set normalization factor for gas and solids pressure residual
      LOGICAL :: SETg
! gas & solids pressure residual
      DOUBLE PRECISION :: RESg
! average velocity
      DOUBLE PRECISION :: Vavg
      CHARACTER(LEN=4) :: TUNIT
! Number of outer iterations goalSeekMassFlux
      INTEGER :: GSMF
      DOUBLE PRECISION :: delP_MF, lMFlux
! Error Message
      CHARACTER(LEN=32) :: lMsg

! initializations
      DT_prev = DT
      NIT = 0
      MUSTIT = 0
      GSMF = 0
      delP_MF = 0.
      lMFlux = 0.

      RESG = ZERO

      SETG = (NORM_G /= ONE)
      NORMG = NORM_G

      LEQ_ADJUST = .FALSE.

! Initialize residuals
      RESID = ZERO

! CPU time left
      IF (FULL_LOG) THEN
         TLEFT = (TSTOP - TIME)*CPUOS
         CALL GET_TUNIT (TLEFT, TUNIT)

         IF (DT == UNDEFINED) THEN
         ELSE
            IF(myPE.eq.PE_IO) THEN
               WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)') &
                  ' Time = ', TIME, '  Dt = ', DT
            ENDIF
         ENDIF   ! if/else(dt==undefined)
      ENDIF   ! if(full_log)

      ! Calculate the face values of densities and mass fluxes
      CALL CONV_ROP(u_g, v_g, w_g, rop_g, rop_ge, rop_gn, rop_gt, flag)
      CALL CALC_MFLUX (u_g, v_g, w_g, rop_ge, rop_gn, rop_gt, &
         flux_ge, flux_gn, flux_gt,flag)
      CALL SET_BC1(p_g,ep_g,ro_g,rop_g,u_g,v_g,w_g,flux_ge,flux_gn,flux_gt,flag)

! Default/Generic Error message
      lMsg = 'Run diverged/stalled'

! Begin iterations
!-----------------------------------------------------------------
   50 CONTINUE
      MUSTIT = 0
      NIT = NIT + 1
! mechanism to set the normalization factor for the correction
! after the first iteration to the corresponding residual found
! in the first iteration
      IF (.NOT.SETG) THEN
         IF (RESG > SMALL_NUMBER) THEN
            NORMG = RESG
            SETG = .TRUE.
         ENDIF
      ENDIF

! Call user-defined subroutine to set quantities that need to be updated
! every iteration
      IF (CALL_USR) CALL USR2

! Calculate coefficients, excluding density and reactions.
      CALL CALC_COEFF(flag, 1, ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g, &
         f_gds, drag_am, drag_bm, pijk, particle_phase, particle_state, pvol, des_pos_new, des_vel_new, des_radius)
      IF (IER_MANAGER()) goto 1000

! Solve starred velocity components
      call solve_vel_star(u_g,v_g,w_g,u_go,v_go,w_go,&
                          p_g,ro_g,rop_g,rop_go,ep_g,&
                          tau_u_g,tau_v_g,tau_w_g,&
                          d_e,d_n,d_t,flux_ge,flux_gn,flux_gt,mu_g,&
                          f_gds, drag_am, drag_bm, flag, IER)


! Calculate densities.
      CALL PHYSICAL_PROP(IER, 0, ro_g, p_g, ep_g, rop_g, ro_g0, flag)
      IF (IER_MANAGER()) goto 1000

! Calculate the face values of densities.
      CALL CONV_ROP(u_g, v_g, w_g, rop_g, rop_ge, rop_gn, rop_gt, flag)

      IF (RO_G0 /= ZERO) THEN
! Solve fluid pressure correction equation
         CALL solve_pp_g (u_g, v_g, w_g, p_g, ep_g, rop_g, rop_go, ro_g, pp_g, &
            rop_ge, rop_gn, rop_gt, d_e, d_n, d_t, flag, NORMG, RESG, IER)

! Correct pressure, velocities, and density
         CALL CORRECT_0 (p_g,pp_g,u_g,v_g,w_g,d_e,d_n,d_t,flag)

      ENDIF

! Recalculate densities.
      CALL PHYSICAL_PROP(IER, 0, ro_g, p_g, ep_g, rop_g, ro_g0, flag)
      IF (IER_MANAGER()) goto 1000

! Update wall velocities:
! modified by sof to force wall functions so even when NSW or FSW are
! declared, default wall BC will still be treated as NSW and no wall
! functions will be used
      CALL SET_WALL_BC (u_g,v_g,w_g, flag)

! Calculate the face values of mass fluxes
      CALL CALC_MFLUX (u_g, v_g, w_g, rop_ge, rop_gn, rop_gt, &
         flux_ge, flux_gn, flux_gt, flag)
      CALL SET_BC1(p_g,ep_g,ro_g,rop_g,u_g,v_g,w_g,flux_ge,flux_gn,flux_gt,flag)

! User-defined linear equation solver parameters may be adjusted after
! the first iteration
      IF (.NOT.CYCLIC) LEQ_ADJUST = .TRUE.


! Check for convergence
      CALL ACCUM_RESID ! Accumulating residuals from all the processors
      RESG = RESID(RESID_P)
      CALL CHECK_CONVERGENCE (NIT, u_g, v_g, w_g, ep_g, 0.0d+0, MUSTIT,flag)

      IF(CYCLIC .AND. (MUSTIT==0 .OR. NIT >= MAX_NIT)) &
         CALL GoalSeekMassFlux(NIT, MUSTIT, GSMF, delP_MF, lMFlux, flux_ge, flux_gn, flux_gt, flag)


!  If not converged continue iterations; else exit subroutine.
 1000 CONTINUE
!-----------------------------------------------------------------

! Display residuals
      CALL DISPLAY_RESID (NIT)

! Determine course of simulation: converge, non-converge, diverge?
      IF (MUSTIT == 0) THEN
! ---------------------------------------------------------------->>>
         IF (DT==UNDEFINED .AND. NIT==1) GOTO 50   !Iterations converged

         IF(CYCLIC .AND. GSMF > 0) THEN
            Write(ERR_MSG,5500) GSMF, delP_MF, lMFlux
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF


! Perform checks and dump to screen every NLOG time steps
         IF (MOD(NSTEP,NLOG) == 0) THEN
            CALL CPU_TIME (CPU_NOW)
            CPUOS = (CPU_NOW - CPU_NLOG)/(TIME - TIME_NLOG)
            CPU_NLOG = CPU_NOW
            TIME_NLOG = TIME
            CPU_NOW = CPU_NOW - CPU0

            WRITE(ERR_MSG,5001) TIME, DT, NIT, 0.0, CPU_NOW
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 5001 FORMAT(1X,'t=',F11.4,' Dt=',G11.4,' NIT=',I3,' Sm=',G12.5, &
         T84,'CPU=',F8.0,' s')

            IF (.NOT.FULL_LOG) THEN
               TLEFT = (TSTOP - TIME)*CPUOS
               CALL GET_TUNIT (TLEFT, TUNIT)
               IF(DMP_LOG)WRITE (UNIT_LOG, '(46X,A,F9.3,1X,A)')
            ENDIF


5500  Format('Mass Flux Iterations:', I0,'   DelP=', &
      G12.5, ' Gas Flux=', G12.5)

            IF (CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z) THEN
               Vavg = VAVG_G(U_G, EP_G, VOL, flag)
               IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'U_g = ', Vavg
               Vavg = VAVG_G(V_G, EP_G, VOL, flag)
               IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'V_g = ',  Vavg
               Vavg = VAVG_G(W_G, EP_G, VOL, flag)
               IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'W_g = ', Vavg
            ENDIF   ! end if cyclic_x, cyclic_y or cyclic_z

         ENDIF   ! end IF (MOD(NSTEP,NLOG) == 0)

         IER = 0
         RETURN   ! for if mustit =0 (converged)
! end converged: go back to time_march
! ----------------------------------------------------------------<<<

! diverged
      ELSEIF (MUSTIT==2 .AND. DT/=UNDEFINED) THEN
! ---------------------------------------------------------------->>>
         IF (FULL_LOG) THEN
            WRITE(ERR_MSG,5200) TIME, DT, NIT, 0.0, trim(adjustl(lMsg))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IER = 1
         RETURN  ! for if mustit =2 (diverged)
      ENDIF
! end diverged: go back to time_march, decrease time step, try again
! ----------------------------------------------------------------<<<

! not converged (mustit = 1, !=0,2 )
! ---------------------------------------------------------------->>>
      IF(NIT < MAX_NIT) THEN
         MUSTIT = 0
         GOTO 50
      ENDIF ! continue iterate
! ----------------------------------------------------------------<<<

      WRITE(ERR_MSG, 5100) TIME, DT, NIT, 0.0
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

! SOF: MFIX will not go the next time step if MAX_NIT is reached,
! instead it will decrease the time step. (IER changed from 0 to 1)
      IER = 1
      RETURN

 5050 FORMAT(5X,'Average ',A,G12.5)
 5100 FORMAT(1X,'t=',F11.4,' Dt=',G11.4,&
         ' NIT>',I3,' Sm= ',G12.5, 'MbErr%=', G11.4)
 5200 FORMAT(1X,'t=',F11.4,' Dt=',G11.4,' NIT=',&
      I3,'MbErr%=', G11.4, ': ',A,' :-(')

      contains

!----------------------------------------------------------------------!
! Function: IER_Manager                                                !
!                                                                      !
! Purpose: Identify and account for errors from called subroutines.    !
!          Returns .TRUE. for lErr >= 100, otherwise .FALSE.           !
!                                                                      !
! Reserved Error Blocks:                                               !
!                                                                      !
! [ 100,  109]: PHYSICAL_PROP                                          !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IER_MANAGER()

! Default case: do nothing.
      IF(IER < 100) THEN
         IER_MANAGER = .FALSE.
         return
      ENDIF

! Errors with an index greater than 100 will force an exit from iterate
! and in turn, reduce the step-size, and restart the time-step.
      IER_MANAGER = .TRUE.
      MUSTIT = 2

! Errors reported from PHYSICAL_PROP
!```````````````````````````````````````````````````````````````````````
      IF(IER <  110) THEN
         IF(IER ==  100) THEN
            lMsg = 'Negative gas density detected'
         ELSEIF(IER ==  101) THEN
            lMsg = 'Negative solids density detected'
         ELSE
            lMsg = 'UCE in PHYSICAL_PROP'
         ENDIF


! Unclassified Errors
!```````````````````````````````````````````````````````````````````````
      ELSE
         lMsg = 'Run diverged/stalled with UCE'
      ENDIF


      IF(DT == UNDEFINED) IER_MANAGER = .FALSE.

      return
      END FUNCTION IER_MANAGER

      END SUBROUTINE ITERATE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!  Purpose:  In the following subroutine the mass flux across a periodic
!            domain with pressure drop is held constant at a
!            user-specified value.  This module is activated only if
!            the user specifies a value for the keyword flux_g in the
!            mfix.dat file.
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      subroutine GoalSeekMassFlux(NIT, MUSTIT, OUTIT, delp_n, mdot_n, flux_ge, flux_gn, flux_gt, flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc, only: delp_x, delp_y, delp_z, flux_g
      USE param1, only: one
      USE compar   ,only: istart3, iend3, jstart3, jend3, kstart3, kend3, myPE, PE_IO
      USE exit_mod, only: mfix_exit
      USE geometry, only: axy, ayz, axz, cyclic_x_mf, cyclic_y_mf, cyclic_z_mf
      USE run     , only: automatic_restart
      USE utilities, ONLY: mfix_isnan
      USE vavg_mod, ONLY: vavg_flux_g

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: flux_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: flux_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer, intent(in   ) :: flag(istart3:iend3,jstart3:jend3,kstart3:kend3,4)
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: NIT, MUSTIT
      INTEGER, INTENT(INOUT) :: OUTIT
      DOUBLE PRECISION, INTENT(INOUT) :: delp_n, mdot_n
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER, PARAMETER :: MAXOUTIT = 500
      DOUBLE PRECISION, PARAMETER :: omega = 0.9
      DOUBLE PRECISION, PARAMETER :: TOL = 1E-03

      DOUBLE PRECISION :: mdot_nm1, delp_nm1
      DOUBLE PRECISION :: delp_xyz

! Store previous values (only used for OUTIT>1)
      mdot_nm1 = mdot_n
      delp_nm1 = delp_n

! Calculate the average gas mass flux and error
      IF(CYCLIC_X_MF)THEN
         delp_n = delp_x
         mdot_n = VAVG_Flux_G(flux_ge, ayz, flag)
      ELSEIF(CYCLIC_Y_MF)THEN
         delp_n = delp_y
         mdot_n = VAVG_Flux_G(flux_gn, axz, flag)
      ELSEIF(CYCLIC_Z_MF)THEN
         delp_n = delp_z
         mdot_n = VAVG_Flux_G(flux_gt, axy, flag)
      ELSE
         RETURN
      ENDIF

      OUTIT = OUTIT + 1

      IF(OUTIT > MAXOUTIT) THEN
         IF (myPE.EQ.PE_IO) write(*,5400) MAXOUTIT
         CALL mfix_exit(0)
      ENDIF

      IF (mfix_isnan(mdot_n) .OR. mfix_isnan(delp_n)) THEN
         IF (myPE.eq.PE_IO) write(*,*) mdot_n, delp_n, &
            ' NaN being caught in GoalSeekMassFlux '
         AUTOMATIC_RESTART = .TRUE.
         RETURN
      ENDIF


! correct delP
      if(OUTIT > 1)then
! Fail-Safe Newton's method
         delp_xyz = delp_n - omega * (delp_n - delp_nm1) * &
            ((mdot_n - FLUX_g)/(mdot_nm1 - FLUX_g)) / &
            ((mdot_n - FLUX_g)/(mdot_nm1 - FLUX_g) - ONE)
      else
         delp_xyz = delp_n*0.99
      endif

! Check for convergence
      IF(abs((mdot_n - FLUX_g)/FLUX_g) < TOL) THEN
         MUSTIT = 0
         delp_n = delp_xyz
      ELSE
        MUSTIT = 1
        NIT = 1
      ENDIF


      IF(CYCLIC_X_MF)THEN
        delp_x = delp_xyz
      ELSEIF(CYCLIC_Y_MF)THEN
        delp_y = delp_xyz
      ELSEIF(CYCLIC_Z_MF)THEN
        delp_z = delp_xyz
      ENDIF

      RETURN

5400 FORMAT(/1X,70('*')//' From: GoalSeekMassFlux',/&
      ' Message: Number of outer iterations exceeded ', I4,/1X,70('*')/)

      END SUBROUTINE GoalSeekMassFlux
end module iterate_module
