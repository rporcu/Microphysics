MODULE CHECK_DATA_30_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: CHECK_DATA_30                                          !
!  Author: M. Syamlal                                 Date: 27-OCT-92  !
!                                                                      !
!  Purpose: Set miscellaneous constants                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DATA_30(lambda_g,mu_g,flag)

! Global variables: (common to sub-functions)
!---------------------------------------------------------------------//
      use compar, only: istart2,iend2, jstart2, jend2, kstart2, kend2
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use param1, only: zero

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: lambda_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) ::     mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

! Local variables:
!---------------------------------------------------------------------//
! Flag for error message headers.
      INTEGER :: ERR_TAG
!......................................................................!


! Check physical properties in inflow/outflow cells.
      CALL CHECK_FLOW_CELL_PROPS(lambda_g,mu_g, flag)

! Verify physical values for field variables.
      CALL CHECK_PHYSICAL_BOUNDS(mu_g, flag)

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
      SUBROUTINE CHECK_FLOW_CELL_PROPS(lambda_g,mu_g,flag)

! Global variables:
!---------------------------------------------------------------------//

      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE open_files_mod, only: open_pe_log, close_pe_log

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: lambda_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

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

         if(flag(i,j,k,1) >= 10 .and. flag(i,j,k,1) <= 31) then

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
      SUBROUTINE CHECK_PHYSICAL_BOUNDS(mu_g, flag)

! Global variables:
!---------------------------------------------------------------------//
      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use open_files_mod, only: open_pe_log, close_pe_log

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) ::     mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) ::     flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

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
         IF (flag(i,j,k,1)<100) THEN

! Gas viscosity must be non-negative.
            IF(MU_G(I,J,K) < ZERO) CALL REPORT_ERROR                     &
               (IER, I, J, K, MU_G(I,J,K), '<', ZERO, 'MU_G')

         ENDIF ! IF(flag(i,j,k,1)<100)
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

      USE open_files_mod, only: open_pe_log, close_pe_log

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
END MODULE CHECK_DATA_30_MODULE
