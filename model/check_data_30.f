!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: CHECK_DATA_30                                          !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Check whether the sum of reaction rates is zero and the    !
!  sum of mass fractions is 1.0 and EP_g >= EP_Star.                   !
!           and EP_g >= EP_star. Set miscellaneous constants           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DATA_30(lambda_g,mu_g)

! Global variables: (common to sub-functions)
!---------------------------------------------------------------------//
      use compar, only: istart2,iend2, jstart2, jend2, kstart2, kend2
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use param1, only: zero, one

      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: lambda_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) ::     mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Local variables:
!---------------------------------------------------------------------//
! Flag for error message headers.
      INTEGER :: ERR_TAG
!......................................................................!


! Check physical properties in inflow/outflow cells.
      CALL CHECK_FLOW_CELL_PROPS(lambda_g,mu_g)

! Verify physical values for field variables.
      CALL CHECK_PHYSICAL_BOUNDS

      RETURN

      CONTAINS

!----------------------------------------------------------------------!
!  Module name: CHECK_FLOW_PROPS                                       !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Verify that inflow/outflow cells do not contain physical   !
!  properties for specified variables.                                 !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_FLOW_CELL_PROPS(lambda_g,mu_g)

! Global variables:
!---------------------------------------------------------------------//

      use functions, only: FLOW_AT
      use compar   , only: DEAD_CELL_AT
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: lambda_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Local variables:
!---------------------------------------------------------------------//
! Loop counters
      INTEGER :: I, J, K
! Integer error flag.
      INTEGER :: IER
!......................................................................!

      CALL INIT_ERR_MSG("CHECK_FLOW_CELL_PROPS")

      IER = 0
      ERR_TAG=2000

      DO K = KSTART2, KEND2
      DO J = JSTART2, JEND2
      DO I = ISTART2, IEND2

         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

         IF(FLOW_AT(i,j,k)) THEN

! Turbulent viscosity of fluid phase.
            IF(MU_g(I,J,K) /= ZERO) CALL REPORT_ERROR                   &
               (IER, I, J, K, MU_G(I,J,K), '/=', ZERO, 'MU_G')

! Granular second coefficient of viscosity.
            IF(LAMBDA_G(i,j,k) /= ZERO) CALL REPORT_ERROR               &
               (IER, I, J, K, LAMBDA_G(i,j,k), '/=', ZERO, 'LAMBDA_G')
! Gas conductivity.

         ENDIF

      ENDDO
      ENDDO
      ENDDO

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,"('End of Report.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
! Close DMP logs when running in interactive mode.
         CALL CLOSE_PE_LOG
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_FLOW_CELL_PROPS



!----------------------------------------------------------------------!
!                                                                      !
!  Module name: CHECK_PHYSICAL_BOUNDS                                  !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Verify that fluid cells have physical values for the       !
!  specified variables.                                                !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_PHYSICAL_BOUNDS

! Global variables:
!---------------------------------------------------------------------//
      use toleranc, only: TMIN, TMAX, TOL_COM
      use fldvar, only: MU_G

      use functions, only: WALL_AT
      use compar, only: DEAD_CELL_AT

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
      INTEGER :: I, J, K
! Integer error flag.
      INTEGER :: IER

      CALL INIT_ERR_MSG("CHECK_PHYSICAL_BOUNDS")

      IER = 0
      ERR_TAG = 3000

      DO K = KSTART2, KEND2
      DO J = JSTART2, JEND2
      DO I = ISTART2, IEND2
         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
         IF (.NOT.WALL_AT(i,j,k)) THEN

! Gas viscosity must be non-negative.
            IF(MU_G(I,J,K) < ZERO) CALL REPORT_ERROR                     &
               (IER, I, J, K, MU_G(I,J,K), '<', ZERO, 'MU_G')

         ENDIF ! IF(.NOT.WALL_AT(i,j,k))
      ENDDO ! DO I = ISTART2, IEND2
      ENDDO ! DO J = JSTART2, JEND2
      ENDDO ! DO K = KSTART2, KEND2

      ! CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,"('End of Report.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
         CALL CLOSE_PE_LOG
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_PHYSICAL_BOUNDS



!----------------------------------------------------------------------!
!  Subroutine: REPORT_ERROR                                            !
!                                                                      !
!  Purpose: Manage error messages for CHECK_DATA_30.                   !
!----------------------------------------------------------------------!
      SUBROUTINE REPORT_ERROR(pIER, pI, pJ, pK, VAL,  RELATION, BND, &
         VAR, LC1, LC2)

      INTEGER, INTENT(INOUT) :: pIER
      INTEGER, INTENT(IN) :: pI, pJ, pK
      DOUBLE PRECISION, INTENT(IN) :: BND, VAL
      CHARACTER(LEN=*), INTENT(IN) :: RELATION
      CHARACTER(LEN=*), INTENT(IN) :: VAR
      INTEGER, INTENT(IN), OPTIONAL :: LC1, LC2
      CHARACTER(LEN=32) :: VAR_FULL

      INTEGER :: lIER

      VAR_FULL=''
      IF(PRESENT(LC2)) THEN
         VAR_FULL = iVAR(VAR,LC1,LC2)
      ELSEIF(PRESENT(LC1)) THEN
         VAR_FULL = iVAR(VAR,LC1)
      ELSE
         VAR_FULL = VAR
      ENDIF

      IF(pIER == 0) THEN
         CALL OPEN_PE_LOG(lIER)
         SELECT CASE(ERR_TAG)
         CASE(2000); WRITE(ERR_MSG,2000)
         CASE(3000); WRITE(ERR_MSG,3000)
         END SELECT
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
         pIER = 1
      ENDIF

 2000 FORMAT('Error 2000: Physical properties detected in flow cells.',&
         2/3X,'I',6x,'J',6x,'K',5x,'Value',8x,A,2x,'Bound',5x,'Variable')

 3000 FORMAT('Error 3000: Unphysical field variables detected.',&
         2/3X,'I',6x,'J',6x,'K',5x,'Value',8x,A,2x,'Bound',5x,'Variable')

      WRITE(ERR_MSG,9000) pI, pJ, pK, VAL, RELATION, BND, trim(VAR_FULL)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 9000 FORMAT(3(I6,1X),G12.4,1X,A,G12.4,1X,A)

      RETURN
      END SUBROUTINE REPORT_ERROR
      END SUBROUTINE CHECK_DATA_30
