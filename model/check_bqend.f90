MODULE CHECK_BATCH_QUEUE_END_MODULE
   CONTAINS
!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_BATCH_QUEUE_END                                   !
!  Author: A.Gel                                      Date:            !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_BATCH_QUEUE_END(pEXIT_SIGNAL)

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use run, only: TERM_BUFFER

      use time_cpu, only: WALL_START

      use machine, only: WALL_TIME

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      use tunit_module, only: get_tunit

      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: pEXIT_SIGNAL

! Logical flags for hault cases.
      LOGICAL :: USER_HAULT
! Elapsed wall time, and fancy formatted buffer/batch queue times.
      real(c_real) :: WALL_STOP, FANCY_BUFF, FANCY_BATCH
! Time units for formatted output.
      CHARACTER(LEN=4) :: WT_UNIT, BF_UNIT, BC_UNIT

! Calculate the current elapsed wall time.
      WALL_STOP = WALL_TIME()
      WALL_STOP = WALL_STOP - WALL_START

! Set flags for wall time exceeded or user specified hault.
      INQUIRE(file="MFIX.STOP", exist=USER_HAULT)

 1100 FORMAT(2/,15('='),' REQUESTED CPU TIME LIMIT REACHED ',('='),/   &
         'Batch Wall Time:',3X,F9.2,1X,A,/'Elapsed Wall Time: ',F9.2,  &
         1X,A,/'Term Buffer:',7X,F9.2,A,/15('='),' REQUESTED CPU ',    &
         'TIME LIMIT REACHED ',('='))

! Report that the hault signal was detected.
      IF(USER_HAULT) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1200 FORMAT(2/,19('='),' MFIX STOP SIGNAL DETECTED ',19('='),/'MFIX.',&
         'STOP file detected in run directory. Terminating MFIX.',/    &
         'Please DO NOT FORGET to erase the MFIX.STOP file before ',   &
         'restarting',/19('='),'MFIX STOP SIGNAL DETECTED ',19('='))

! This routine was restructured so all MPI ranks to the same action. As
! a result, broadcasting the BATCHQ flag may not be needed.
      pEXIT_SIGNAL = (USER_HAULT) .OR. pEXIT_SIGNAL
      ! call bcast (pEXIT_SIGNAL,PE_IO)

      END SUBROUTINE CHECK_BATCH_QUEUE_END
END MODULE CHECK_BATCH_QUEUE_END_MODULE
