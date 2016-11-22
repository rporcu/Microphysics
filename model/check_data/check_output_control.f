!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_OUTPUT_CONTROL                                    !
!  Purpose: Check the output control namelist section                  !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_OUTPUT_CONTROL

! Global Variables:
!---------------------------------------------------------------------//
! Time intervalue between updating the RES file.
      use output, only: RES_DT
! Time-step intervalue between updating the .LOG file.
      use output, only: NLOG

! Global Parameters:
!---------------------------------------------------------------------//
! Number aliases
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, LARGE_NUMBER

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_OUTPUT_CONTROL")


! Check the values specified for the RES file.
      IF (RES_DT==UNDEFINED)THEN
         WRITE(ERR_MSG,1000) 'RES_DT'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(RES_DT <= ZERO) THEN
         WRITE(ERR_MSG,1002) 'RES_DT', RES_DT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Verify that the LOG frequency is valid.
      IF(NLOG <= 0) THEN
         WRITE(ERR_MSG,1003) 'NLOG', NLOG
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Finalize the error manager.
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',E14.6,/  &
         'Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: Illegal or unknown input: ',A,' = ',I4,/     &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_OUTPUT_CONTROL
