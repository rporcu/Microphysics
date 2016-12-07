!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_LIN_EQ                                            C
!  Purpose: Interface for linear equation solver                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_LIN_EQ(VNAME, Vno, VAR, A_M, B_M, M, ITMAX,&
                              METHOD, SWEEP, TOL1, PC, IER)

         use compar, only: istart3, iend3
         use compar, only: jstart3, jend3
         use compar, only: kstart3, kend3
         use compar, only: mype
         use exit_mod, only: mfix_exit
         use param, only: dimension_3
         use residual, only: resid
         use toleranc, only: tol_resid

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number
!     Note: not really used beyond this subroutine. here it is
!     used for potentially adjusting the tolerances but it is
!     currently disabled code.
!     1 = pressure correction equation
!     3 = gas u-momentum
!     4 = gas v-momentum
!     5 = gas w-momentum
      INTEGER, INTENT(IN) :: Vno
! variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! phase index
      INTEGER, INTENT(IN) :: M
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! linear equation solver method (generally leq_method)
!     2 = bicgstab (default)
      INTEGER, INTENT(IN) :: METHOD
! sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=4), INTENT(IN) :: SWEEP
! convergence tolerance for leq solver (leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL1
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'
      CHARACTER(LEN=4), INTENT(IN) :: PC
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! Adjust LEQ tolerance flag
      LOGICAL, PARAMETER :: adjust_leq_tol = .FALSE.
      LOGICAL, PARAMETER :: leq_tol_scheme1 = .FALSE.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!
      DOUBLE PRECISION :: max_resid_local, tol_resid_max
! convergence tolerance for leq solver
      DOUBLE PRECISION :: TOL
! for constructing local character strings
      CHARACTER(LEN=80) :: LINE0, LINE1
!-----------------------------------------------

! Adjusting the tolerances
! ---------------------------------------------------------------->>>
      IF(adjust_leq_tol) THEN
         max_resid_local = maxval(resid(:),1)
         tol_resid_max   = TOL_RESID
         IF(leq_tol_scheme1.AND.resid(Vno).LT.1.0D-1) THEN
            if(Vno.le.5) then
               TOL = MAX(TOL1,TOL1*RESID(Vno)/TOL_RESID)
            endif
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, resid(Vno)
         ELSEIF(max_resid_local.LT.1.0D-1) THEN
            TOL = MAX(TOL1,TOL1*max_resid_local/TOL_RESID_max)
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, max_resid_local
         ENDIF
      ELSE
        TOL = TOL1
      ENDIF
! ----------------------------------------------------------------<<<


! Solve the linear system of equations
! ---------------------------------------------------------------->>>
      SELECT CASE (METHOD)

      CASE (2)
            call leq_bicgs(VNAME, VNO, VAR, A_M, B_M,&
                           SWEEP, TOL, PC, ITMAX, IER)

      CASE DEFAULT
         LINE0(1:14) = 'SOLVE_LIN_EQ: '
         LINE0(15:80)= VName
         WRITE(LINE1,'(A, I2, A)') &
             'Error: LEQ_METHOD = ', METHOD, ' is invalid'
         CALL WRITE_ERROR(LINE0, LINE1, 1)
         CALL mfix_exit(myPE)
      END SELECT
! ----------------------------------------------------------------<<<

      RETURN
      END SUBROUTINE SOLVE_LIN_EQ
