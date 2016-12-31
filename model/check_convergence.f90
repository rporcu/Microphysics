MODULE CHECK_CONVERGENCE_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_CONVERGENCE                                       C
!  Author: M. Syamlal                                 Date: 8-JUL-96   C
!                                                                      C
!  Purpose: Monitor convergence                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      integer function check_convergence(nit) &
         bind(C, name="check_convergence")

      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE param1, only: zero, undefined_i, is_undefined
      USE residual, only: max_resid_index, nresid
      use residual, only: resid_p, resid_u, resid_v, resid_w
      USE residual, only: resid, resid_index, resid_string, resid_x
      use residual, only: sum5_resid, group_resid, resid_prefix, resid_grp, hydro_grp
      USE run, only: detect_stall
      USE toleranc, only: tol_resid, tol_diverge
      USE utilities, ONLY: check_vel_bound
      use iso_c_binding, only: c_double, c_int

      IMPLICIT NONE

! Dummy arguments
! Iteration number
      integer(c_int), intent(in) :: nit

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

! add momentum equation residuals
      SUM = SUM + RESID(RESID_U)
      SUM = SUM + RESID(RESID_V)
      SUM = SUM + RESID(RESID_W)


! find the variable with maximum residual
      MAXM = 0
      IF (IS_UNDEFINED(RESID_INDEX(MAX_RESID_INDEX,1))) THEN
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
         IF (IS_UNDEFINED(MAXN)) THEN
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
               check_convergence = 2
               RETURN
            ENDIF
         ENDIF
         SUM5_RESID = SUM
      ENDIF

! Require at least two iterations.
      IF(NIT == 1) THEN
         CHECK_CONVERGENCE = 0
         RETURN
      ENDIF

! total residual
      IF(SUM<=TOL_RESID) THEN
         CHECK_CONVERGENCE = 1          ! converged
      ELSEIF (SUM>=TOL_DIVERGE ) THEN
         IF (NIT /= 1) THEN
            CHECK_CONVERGENCE = 2       ! diverged
         ELSE
            CHECK_CONVERGENCE = 0       ! not converged
         ENDIF
      ELSE
         CHECK_CONVERGENCE = 0          ! not converged
      ENDIF

! Check upper bound (speed of sound) limit for gas velocity components.
!      IF(MOMENTUM_X_EQ(0) .OR. MOMENTUM_Y_EQ(0) .OR. &
!          MOMENTUM_Z_EQ(0)) THEN
!         IF(CHECK_VEL_BOUND(u_g,v_g,w_g,ep_g,flag)) CHECK_CONVERGENCE = 2     !divergence
!      ENDIF

      RETURN
   END function check_convergence
END MODULE CHECK_CONVERGENCE_MODULE
