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

      SUBROUTINE TIME_MARCH

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE cutcell
      USE dashboard
      USE discretelement
      USE drag
      USE fldvar
      USE funits
      USE geometry
      USE indices
      USE leqsol, only: SOLVER_STATISTICS, REPORT_SOLVER_STATS
      USE param
      USE param1
      USE physprop
      USE run
      USE time_cpu
      USE toleranc
! use function MAX_VEL_INLET to compute max. velocity at inlet
      USE utilities, ONLY: MAX_VEL_INLET
      USE vtk
      USE vtp
      use output, only: RES_DT
      use adjust_dt

      IMPLICIT NONE
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

! C Function
      INTERFACE
         SUBROUTINE CHECK_SOCKETS() BIND ( C )
           use, INTRINSIC :: iso_c_binding
         END SUBROUTINE CHECK_SOCKETS
      END INTERFACE

!-----------------------------------------------

      FINISH  = .FALSE.
      NCHECK  = NSTEP
      DNCHECK = 1
      CPU_IO  = ZERO
      NIT_TOTAL = 0

      CALL INIT_OUTPUT_VARS

! Parse residual strings
      CALL PARSE_RESID_STRING ()

! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0

! Calculate all the coefficients once before entering the time loop
      CALL CALC_COEFF(IER, 2)
      IF(MU_g0 == UNDEFINED) CALL CALC_MU_G()

! Remove undefined values at wall cells for scalars
      where(rop_g == undefined) rop_g = 0.0

! Initialize d's and e's to zero
      F_gs = 0.0d0
      D_E = 0.0d0
      D_N = 0.0d0
      D_T = 0.0d0

! The TIME loop begins here.............................................
 100  CONTINUE


! Terminate MFIX normally before batch queue terminates.
      IF (CHK_BATCHQ_END) CALL CHECK_BATCH_QUEUE_END(EXIT_SIGNAL)

      IF (CALL_USR) CALL USR1

! Set wall boundary conditions and transient flow b.c.'s
      CALL SET_BC1

      CALL OUTPUT_MANAGER(EXIT_SIGNAL, FINISH)

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
      CALL UPDATE_OLD

! Calculate coefficients
      CALL CALC_COEFF_ALL (0, IER)

! Calculate the stress tensor trace and cross terms for all phases.
      CALL CALC_TRD_AND_TAU()

! Check rates and sums of mass fractions every NLOG time steps
      IF (NSTEP == NCHECK) THEN
         IF (DNCHECK < 256) DNCHECK = DNCHECK*2
         NCHECK = NCHECK + DNCHECK
         CALL CHECK_DATA_30
      ENDIF

! Check for maximum velocity at inlet to avoid convergence problems
      MAX_INLET_VEL = 100.0d0*MAX_VEL_INLET()
! if no inlet velocity is specified, use an upper limit defined in
! toleranc_mod.f
      IF(MAX_INLET_VEL <= SMALL_NUMBER) THEN
         MAX_INLET_VEL = MAX_ALLOWED_VEL
         IF (UNITS == 'SI') MAX_INLET_VEL = 1D-2 * MAX_ALLOWED_VEL
      ENDIF
! Scale the value using a user defined scale factor
      MAX_INLET_VEL = MAX_INLET_VEL * MAX_INLET_VEL_FAC

! Advance the solution in time by iteratively solving the equations
      CALL ITERATE (IER, NIT)


      DO WHILE (ADJUSTDT(IER,NIT))
         CALL ITERATE (IER, NIT)
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

            IF(WRITE_DASHBOARD) THEN
               RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
               CALL UPDATE_DASHBOARD(NIT,0.0d0,'    ')
            ENDIF
            CALL MFIX_EXIT(MyPE)
         ENDIF
      ENDIF



! Other solids model implementations
      IF(DEM_SOLIDS) THEN
         CALL DES_TIME_MARCH
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
