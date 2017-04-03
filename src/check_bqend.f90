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

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use time_cpu, only: WALL_START

      use machine, only: WALL_TIME

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      use tunit_module, only: get_tunit

      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: pEXIT_SIGNAL

! Logical flags for hault cases.
      LOGICAL :: useR_HAULT
! Elapsed wall time, and fancy formatted buffer/batch queue times.
      real(c_real) :: WALL_STOP

! Calculate the current elapsed wall time.
      WALL_STOP = WALL_TIME()
      WALL_STOP = WALL_STOP - WALL_START

! Set flags for wall time exceeded or user specified hault.
      INQUIRE(file="MFIX.STOP", exist=useR_HAULT)

! Report that the hault signal was detected.
      IF(useR_HAULT) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1200 FORMAT(2/,19('='),' MFIX STOP SIGNAL DETECTED ',19('='),/'MFIX.',&
         'STOP file detected in run directory. Terminating MFIX.',/    &
         'Please DO NOT FORGET to erase the MFIX.STOP file before ',   &
         'restarting',/19('='),'MFIX STOP SIGNAL DETECTED ',19('='))

! This routine was restructured so all MPI ranks to the same action. As
! a result, broadcasting the BATCHQ flag may not be needed.
      pEXIT_SIGNAL = (useR_HAULT) .OR. pEXIT_SIGNAL
      ! call bcast (pEXIT_SIGNAL,PE_IO)

      END SUBROUTINE CHECK_BATCH_QUEUE_END
END MODULE CHECK_BATCH_QUEUE_END_MODULE
