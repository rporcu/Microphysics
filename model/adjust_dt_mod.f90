      MODULE ADJUST_DT

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ADJUST_DT(IER, NIT)                                    !
!  Author: M. Syamlal                                 Date: FEB-10-97  !
!                                                                      !
!  Purpose: Automatically adjust time step.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         LOGICAL FUNCTION ADJUSTDT (ep_g, ep_go, p_g, p_go, ro_g, ro_go, flag, vol_surr, &
            rop_g, rop_go, U_g,  U_go, V_g, V_go,  W_g,  W_go, mu_g, &
            f_gds, drag_am, drag_bm,  particle_phase, iglobal_id, &
            particle_state, pinc, pmass, pvol, des_radius,  &
            des_pos_new, des_vel_new, des_usr_var, IER, NIT)

! Global Variables:
!---------------------------------------------------------------------//
! User defined aximum number of iterations
      use leqsol, only: MAX_NIT
! User defined: min, max DT and adjustment factor
      use run, only: DT_MIN, DT_MAX, DT_FAC
! Flag: Use stored DT value for advancing TIME
      use run, only: USE_DT_PREV
! Current DT (1/DT) and direction of last change (+/-)
      use run, only: DT, oDT, DT_DIR

      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, ONE, UNDEFINED

! Module proceedures:
!---------------------------------------------------------------------//
! Routine to break successive time step reductions.
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ival

      use calc_coeff_module, only: calc_coeff_all

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
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
      DOUBLE PRECISION, INTENT(IN   ) :: p_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ep_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(  OUT) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(  OUT) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(  OUT) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      integer         , INTENT(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
      DOUBLE PRECISION, INTENT(in   ) :: vol_surr&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer         , INTENT(inout) :: pinc&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pvol, pmass, des_radius
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new, des_usr_var
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: iglobal_id
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase

! Integer flag: 0=Good, 100=initialize, otherwise bad.
      INTEGER, INTENT(INOUT) :: IER
! Number of iterations for current time step
      INTEGER, INTENT(IN) :: NIT


! Local Variables:
!---------------------------------------------------------------------//
! Number of steps in between DT adjustments.
      INTEGER, PARAMETER :: STEPS_MIN = 5
! Number of time steps since last DT adjustment
      INTEGER, SAVE :: STEPS_TOT=0
! number of iterations since last DT adjustment
      INTEGER, SAVE :: NIT_TOT=0
! Iterations per second for last dt
      DOUBLE PRECISION, SAVE :: NIToS=0.0
! Current number of iterations per second
      DOUBLE PRECISION :: NITOS_NEW
!......................................................................!

! Initialize the function result.
      ADJUSTDT = .FALSE.
      USE_DT_PREV = .FALSE.

! Steady-state simulation.
      IF (DT==UNDEFINED .OR. DT<ZERO) RETURN

! Iterate successfully converged.
!---------------------------------------------------------------------//
      IF(IER == 0) THEN

! Calculate a new DT every STEPS_MIN time steps.
         IF(STEPS_TOT >= STEPS_MIN) THEN
            NITOS_NEW = DBLE(NIT_TOT)/(STEPS_TOT*DT)
            IF (NITOS_NEW > NITOS) DT_DIR = DT_DIR*(-1)
            STEPS_TOT = 0
            NITOS = NITOS_NEW
            NIT_TOT = 0
            IF (DT_DIR > 0) THEN
               IF(NIT < MAX_NIT) DT = MIN(DT_MAX,DT/DT_FAC)
            ELSE
               DT = DT*DT_FAC
            ENDIF

! DT was modified. Use the stored DT should be used to update TIME.
            USE_DT_PREV = .TRUE.

! Write the convergence stats to the screen/log file.
            WRITE(ERR_MSG,"('DT=',g11.4,3x,'NIT/s=',A)")  &
               DT, trim(iVal(nint(NITOS)))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., &
               FOOTER=.FALSE., LOG=.FALSE.)

         ELSE
            STEPS_TOT = STEPS_TOT + 1
            NIT_TOT = NIT_TOT + NIT
         ENDIF
! No need to iterate again
         ADJUSTDT = .FALSE.

! Iterate failed to converge.
!---------------------------------------------------------------------//
      ELSE

! Clear the error flag.
         IER = 0

! Reset counters.
         STEPS_TOT = 0
         NITOS = 0.
         NIT_TOT = 0

! Reduce the step size.
         DT = DT*DT_FAC

! The step size has decreased to the minimum.
         IF (DT_FAC >= ONE) THEN

            WRITE(ERR_MSG,"(3X,A)") &
               'DT_FAC >= 1. Recovery not possible!'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE., &
               HEADER=.FALSE., FOOTER=.FALSE.)

         ELSEIF (DT > DT_MIN) THEN

            WRITE(ERR_MSG,"(3X,'Recovered: Dt=',G12.5,' :-)')") DT
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

            ep_g = ep_go
            p_g =  p_go
            ro_g = ro_go
            rop_g = rop_go
            U_g =  U_go
            V_g =  V_go
            W_g =  W_go

            ! Recalculate all coefficients
            CALL CALC_COEFF_ALL (ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g,&
               f_gds, drag_am, drag_bm,  particle_phase, iglobal_id, &
               particle_state, pmass, pvol, des_pos_new, des_vel_new, des_radius,  &
               des_usr_var, flag, vol_surr, pinc)
! Iterate again with new dt
            ADJUSTDT = .TRUE.

! Set the return flag stop iterating.
         ELSE

! Prevent DT from dropping below DT_MIN.
            ADJUSTDT = .FALSE.
         ENDIF

      ENDIF

      ! Update ONE/DT variable.
      ODT = ONE/DT

      RETURN
      END FUNCTION ADJUSTDT

      END MODULE ADJUST_DT
