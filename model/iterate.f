!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ITERATE                                                 C
!  Author: M. Syamlal                                 Date: 12-APR-96  C
!                                                                      C
!  Purpose: This module controls the iterations for solving equations  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ITERATE(IER, NIT)

      USE compar
      USE cutcell
      USE dashboard
      USE discretelement
      USE fldvar
      USE funits
      USE geometry
      USE indices
      USE leqsol
      USE machine, only: start_log, end_log
      USE mpi_utility
      USE output
      USE param
      USE param1
      USE physprop
      USE residual
      USE run
      USE time_cpu
      USE toleranc
      USE vavg_mod, ONLY: vavg_g
      USE vtk

      use error_manager

      IMPLICIT NONE
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

! Error Message
      CHARACTER(LEN=32) :: lMsg

! initializations
      DT_prev = DT
      NIT = 0
      MUSTIT = 0
      RESG = ZERO

      SETG = (NORM_G /= ONE)
      NORMG = NORM_G

      LEQ_ADJUST = .FALSE.

! Initialize residuals
      RESID = ZERO

! Initialize the routine for holding gas mass flux constant with cyclic bc
      IF(CYCLIC) CALL GoalSeekMassFlux(NIT, MUSTIT, .false.)

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

! Calculate the face values of densities and mass fluxes for the first
! solve_vel_star call.
      CALL CONV_ROP()
      CALL CALC_MFLUX ()
      CALL SET_BC1

! JFD: modification for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_OUTFLOW

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
      CALL CALC_COEFF(IER, 1)
      IF (IER_MANAGER()) goto 1000

! Solve starred velocity components
      CALL SOLVE_VEL_STAR(IER)

! Calculate densities.
      CALL PHYSICAL_PROP(IER, 0)
      IF (IER_MANAGER()) goto 1000

! Calculate the face values of densities.
      CALL CONV_ROP()

      IF (RO_G0 /= ZERO) THEN
! Solve fluid pressure correction equation
         CALL SOLVE_PP_G (NORMG, RESG, IER)
! Correct pressure, velocities, and density
         CALL CORRECT_0 ()
      ENDIF

! Recalculate densities.
      CALL PHYSICAL_PROP(IER, 0)
      IF (IER_MANAGER()) goto 1000

! Update wall velocities:
! modified by sof to force wall functions so even when NSW or FSW are
! declared, default wall BC will still be treated as NSW and no wall
! functions will be used
      CALL SET_WALL_BC ()

! Calculate the face values of mass fluxes
      CALL CALC_MFLUX ()
      CALL SET_BC1

! JFD: modification for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_OUTFLOW


! User-defined linear equation solver parameters may be adjusted after
! the first iteration
      IF (.NOT.CYCLIC) LEQ_ADJUST = .TRUE.


! Check for convergence
      CALL ACCUM_RESID ! Accumulating residuals from all the processors
      RESG = RESID(RESID_P,0)
      CALL CHECK_CONVERGENCE (NIT, 0.0d+0, MUSTIT)

      IF(CYCLIC)THEN
        IF(MUSTIT==0 .OR. NIT >= MAX_NIT) &
           CALL GoalSeekMassFlux(NIT, MUSTIT, .true.)
      ENDIF


!  If not converged continue iterations; else exit subroutine.
 1000 CONTINUE
!-----------------------------------------------------------------

! Display residuals
      CALL DISPLAY_RESID (NIT)

! Determine course of simulation: converge, non-converge, diverge?
      IF (MUSTIT == 0) THEN
! ---------------------------------------------------------------->>>
         IF (DT==UNDEFINED .AND. NIT==1) GOTO 50   !Iterations converged

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

            IF (CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z) THEN
               IF (DO_I) THEN
                 Vavg = VAVG_G(U_G, VOL_U)
                 IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'U_g = ', Vavg
               ENDIF
               IF (DO_J) THEN
                 Vavg = VAVG_G(V_G, VOL_V)
                 IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'V_g = ',  Vavg
               ENDIF
               IF (DO_K) THEN
                 Vavg = VAVG_G(W_G, VOL_W)
                 IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'W_g = ', Vavg
               ENDIF
            ENDIF   ! end if cyclic_x, cyclic_y or cyclic_z

            CALL END_LOG
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



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Purpose:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_TUNIT(TLEFT, TUNIT)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      DOUBLE PRECISION, INTENT(INOUT) :: TLEFT
      CHARACTER(LEN=4) :: TUNIT
!-----------------------------------------------

      IF (TLEFT < 3600.0d0) THEN
         TUNIT = 's'
      ELSE
         TLEFT = TLEFT/3600.0d0
         TUNIT = 'h'
         IF (TLEFT >= 24.) THEN
            TLEFT = TLEFT/24.0d0
            TUNIT = 'days'
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE GET_TUNIT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!  Purpose:  In the following subroutine the mass flux across a periodic
!            domain with pressure drop is held constant at a
!            user-specified value.  This module is activated only if
!            the user specifies a value for the keyword flux_g in the
!            mfix.dat file.
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine GoalSeekMassFlux(NIT, MUSTIT, doit)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE geometry
      USE mflux, ONLY: flux_ge, flux_gn, flux_gt
      USE run
      USE time_cpu
      USE utilities, ONLY: mfix_isnan
      USE vavg_mod, ONLY: vavg_flux_g
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: NIT, MUSTIT
      LOGICAL, INTENT(IN) :: doit
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER, PARAMETER :: MAXOUTIT = 500
      DOUBLE PRECISION, PARAMETER          :: omega = 0.9
      DOUBLE PRECISION, PARAMETER          :: TOL = 1E-03
      INTEGER, SAVE :: OUTIT
      LOGICAL, SAVE :: firstPass = .true.

      DOUBLE PRECISION, SAVE  :: mdot_n, mdot_nm1, delp_n, delp_nm1, err
      DOUBLE PRECISION        :: mdot_0, delp_xyz

      IF(CYCLIC_X_MF)THEN
         delp_n = delp_x
      ELSEIF(CYCLIC_Y_MF)THEN
         delp_n = delp_y
      ELSEIF(CYCLIC_Z_MF)THEN
         delp_n = delp_z
      ELSE
         RETURN
      ENDIF

      IF(.NOT.doit) THEN
         OUTIT = 0
         RETURN
      ENDIF

      OUTIT = OUTIT + 1
      IF(OUTIT > MAXOUTIT) THEN
         IF (myPE.EQ.PE_IO) write(*,5400) MAXOUTIT
         CALL mfix_exit(0)
      ENDIF

      mdot_0 = Flux_g

      ! calculate the average gas mass flux and error
      IF(CYCLIC_X_MF)THEN
        mdot_n = VAVG_Flux_G(flux_ge, ayz)
      ELSEIF(CYCLIC_Y_MF)THEN
        mdot_n = VAVG_Flux_G(flux_gn, axz)
      ELSEIF(CYCLIC_Z_MF)THEN
        mdot_n = VAVG_Flux_G(flux_gt, axy)
      ENDIF

      IF (mfix_isnan(mdot_n) .OR. mfix_isnan(delp_n)) THEN
         IF (myPE.eq.PE_IO) write(*,*) mdot_n, delp_n, &
            ' NaN being caught in GoalSeekMassFlux '
         AUTOMATIC_RESTART = .TRUE.
         RETURN
      ENDIF

      err = abs((mdot_n - mdot_0)/mdot_0)
      IF( err < TOL) THEN
         MUSTIT = 0
      ELSE
        MUSTIT = 1
        NIT = 1
      ENDIF

! correct delp
      if(.not.firstPass)then
!        delp_xyz = delp_n - omega * (delp_n - delp_nm1) * (mdot_n - mdot_0) &
!                          / (mdot_n - mdot_nm1)
! Fail-Safe Newton's method (below) works better than the regular
! Newton method (above)

         delp_xyz = delp_n - omega * (delp_n - delp_nm1) * &
                     ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0)) / &
                     ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0) - ONE)
      else
         firstPass=.false.
         delp_xyz = delp_n*0.99
      endif

      IF(MUSTIT == 0) then
        IF(myPE.eq.PE_IO) Write(*,5500) TIME, OUTIT, delp_xyz, mdot_n
      ENDIF

      mdot_nm1 = mdot_n
      delp_nm1 = delp_n

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
5500  Format('  Time=', G12.5, ' MassFluxIterations=', I4, ' DelP=', &
      G12.5, ' Gas Flux=', G12.5)

      END SUBROUTINE GoalSeekMassFlux
