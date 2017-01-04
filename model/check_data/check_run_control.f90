MODULE CHECK_RUN_CONTROL_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_RUN_CONTROL                                       !
!  Purpose: Check the run control namelist section                     !
!                                                                      !
!  Author: P.Nicoletti                                Date: 27-NOV-91  !
!          J.Musser                                   Date: 31-JAN-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_RUN_CONTROL(time, dt)


! Global Variables:
!---------------------------------------------------------------------//
! New or restart
      USE run, only: RUN_TYPE
! Brief description of simulation.
      USE run, only: DESCRIPTION
! Simulation units: SI, CGS
      USE run, only: UNITS
! Simulation start/stop times.
      USE run, only: TSTOP

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED_C, IS_UNDEFINED
      USE param1, only: ZERO

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg
      use error_manager, only: ivar, ival, err_msg

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
      double precision, intent(inout) :: time, dt

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_RUN_CONTROL")


! Clear out the run description if not specified.
      IF (DESCRIPTION == UNDEFINED_C) DESCRIPTION = ' '

! Verify UNITS input.
      IF(UNITS == UNDEFINED_C) THEN
         WRITE(ERR_MSG,1000) 'UNITS'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF((UNITS /= 'CGS') .AND. (UNITS /= 'SI')) THEN
         WRITE(ERR_MSG,1001) 'UNITS', UNITS
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Verify that DT is valid.
      IF (DT < ZERO) THEN
         WRITE(ERR_MSG,1002) 'DT', DT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

! Steady-state simulation.
      ELSEIF(IS_UNDEFINED(DT) .OR. IS_UNDEFINED(DT)) THEN
         TIME = ZERO

! Transient simulation.
      ELSE
! Verify the remaining time settings.
         IF (IS_UNDEFINED(TIME)) THEN
            WRITE(ERR_MSG,1000) 'TIME'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (IS_UNDEFINED(TSTOP)) THEN
            WRITE(ERR_MSG,1000) 'TSTOP'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (TIME < ZERO) THEN
            WRITE(ERR_MSG,1002)'TIME', TIME
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (TSTOP < ZERO) THEN
            WRITE(ERR_MSG,1002) 'TSTOP', TSTOP
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

! Verify the run type.
      IF(.NOT.(RUN_TYPE=='NEW')) THEN
         WRITE(ERR_MSG,1001) 'RUN_TYPE', RUN_TYPE
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Clear the error manager
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',G14.4,/  &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_RUN_CONTROL
END MODULE CHECK_RUN_CONTROL_MODULE
