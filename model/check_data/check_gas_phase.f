!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_GAS_PHASE                                        !
!  Purpose: Check the gas phase input section                          !
!                                                                      !
!  Author: P.Nicoletti                                Date: 02-DEC-91  !
!          J.Musser                                   Date: 01-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GAS_PHASE


! Global Variables:
!---------------------------------------------------------------------//
! User specified: Constant gas viscosity
      use fldvar, only: MU_G0
! User specified: Constant gas density
      use fldvar, only: RO_G0
! User specified: Constant gas mixture molecular weight
      use fldvar, only: MW_AVG


! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: UNDEFINED, ZERO

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! NONE

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GAS_PHASE")


! CHECK MU_g0
      IF (MU_G0 <= ZERO) THEN
         WRITE(ERR_MSG,1001) 'MU_G0', iVal(MU_G0)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! CHECK MW_AVG
! When the species equations are not solved and the gas phase is
! compressible, verify that the user provided average molecular weight
! has a physical value. (This does not include the case where MW_AVG
! is UNDEFINED.)
      IF (RO_G0 == UNDEFINED) THEN
         IF (MW_AVG == UNDEFINED) THEN
            WRITE(ERR_MSG, 1001) 'MW_AVG', iVal(MW_AVG)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ELSE
! Gas density for incompressible flows must be positive.
         IF (RO_G0 < ZERO) THEN
            WRITE(ERR_MSG, 1001) 'RO_G0', iVal(RO_G0)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
! Incompressible simulations do not need MW_AVG. Notify the user that
! the provided data is ignored.
         IF (MW_AVG /= UNDEFINED)THEN
            WRITE(ERR_MSG, 1100) 'RO_g0 is specified'
            CALL FLUSH_ERR_MSG
         ENDIF

      ENDIF


! Finalize the error manager
      CALL FINL_ERR_MSG


      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')


 1100 FORMAT('Message 2000: MW_AVG is not needed when ',A,'.')

      END SUBROUTINE CHECK_GAS_PHASE
