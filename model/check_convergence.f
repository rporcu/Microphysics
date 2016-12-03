!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_CONVERGENCE                                       C
!  Purpose: Monitor convergence                                        C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 8-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_CONVERGENCE(NIT, u_g, v_g, w_g, ep_g, errorpercent, MUSTIT)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE geometry
      USE param
      USE param1
      USE constant
      USE residual
      USE run
      USE toleranc
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE utilities, ONLY: check_vel_bound

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! Maximum % error allowed in fluid continuity
!      DOUBLE PRECISION, PARAMETER :: MaxErrorPercent = 1.0E-6
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Iteration number
      INTEGER, INTENT(IN) :: NIT
! %error in fluid mass balance
      DOUBLE PRECISION, INTENT(IN) :: errorpercent
! value tells whether to iterate (1) or not (0).
      INTEGER, INTENT(INOUT) :: MUSTIT
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! sum of residuals
      DOUBLE PRECISION :: SUM
! max of residuals
      DOUBLE PRECISION :: maxres
! index
      INTEGER :: L, maxL, maxM, maxN
!-----------------------------------------------

! sum the residuals from correction equation (pressure and/or
! solids), continuity equations (gas and/or solids) and momentum
! equations (gas and solids)

! add pressure correction residual
      SUM = RESID(RESID_P)

! add continuity equation residuals
      SUM = SUM + RESID(RESID_RO)

! add momentum equation residuals
      SUM = SUM + RESID(RESID_U)
      SUM = SUM + RESID(RESID_V)
      SUM = SUM + RESID(RESID_W)


! find the variable with maximum residual
      MAXM = 0
      IF (RESID_INDEX(MAX_RESID_INDEX,1) == UNDEFINED_I) THEN
         MAXRES = ZERO
         DO L = 1, NRESID
            IF (RESID(L) >= MAXRES) THEN
               MAXRES = RESID(L)
               MAXL = L
               IF (L >= RESID_X) THEN
                  MAXN = L - RESID_X + 1
               ELSE
                  MAXN = UNDEFINED_I
               ENDIF
            ENDIF
         ENDDO
         IF (MAXN == UNDEFINED_I) THEN
            WRITE (RESID_STRING(MAX_RESID_INDEX), '(A1,I1)') &
               RESID_PREFIX(MAXL), MAXM
         ELSE
            WRITE (RESID_STRING(MAX_RESID_INDEX), '(A1,I1,I2.0)') &
               'X', MAXM, MAXN
         ENDIF
      ENDIF
      IF (GROUP_RESID) RESID_GRP(HYDRO_GRP) = SUM

! Every 5 iterations detect whether the run is stalled by checking
! that the total residual has decreased.
      IF(DETECT_STALL .AND. MOD(NIT,5) == 0) THEN
         IF(NIT > 10) THEN
            IF(SUM5_RESID <= SUM) THEN
! The run is stalled. Reduce the time step.
               MUSTIT = 2
               RETURN
            ENDIF
         ENDIF
         SUM5_RESID = SUM
      ENDIF

! Require at least two iterations.
      IF(NIT == 1) THEN
         MUSTIT = 1
         RETURN
      ENDIF

! total residual
      IF(SUM<=TOL_RESID) THEN
         MUSTIT = 0                              !converged
      ELSEIF (SUM>=TOL_DIVERGE ) THEN
         IF (NIT /= 1) THEN
            MUSTIT = 2                           !diverged
         ELSE
            MUSTIT = 1                           !not converged
         ENDIF
      ELSE
         MUSTIT = 1                              !not converged
      ENDIF

! Check upper bound (speed of sound) limit for gas velocity components.
      IF(MOMENTUM_X_EQ(0) .OR. MOMENTUM_Y_EQ(0) .OR. &
          MOMENTUM_Z_EQ(0)) THEN
         IF(CHECK_VEL_BOUND(u_g,v_g,w_g,ep_g)) MUSTIT = 2     !divergence
      ENDIF

      RETURN
      END SUBROUTINE CHECK_CONVERGENCE
