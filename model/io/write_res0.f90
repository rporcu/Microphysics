!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: WRITE_RES0                                              C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!                                                                      C
!  Purpose: write out the initial restart records (namelist data)      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE WRITE_RES0

      USE funits, only: unit_res
      USE machine, only: id_month, id_day, id_year, id_hour, id_minute, id_second
      USE run, only: run_name

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Pointer to the next record
      INTEGER :: NEXT_RECA
! file version id
      CHARACTER(LEN=512) :: VERSION
!-----------------------------------------------

      NEXT_RECA = 5

      VERSION = 'RES = 01.8'

!------------------------------------------------------------------------
      WRITE (UNIT_RES, REC=1) VERSION
      WRITE (UNIT_RES, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, &
         ID_MINUTE, ID_SECOND
      WRITE (UNIT_RES, REC=3) NEXT_RECA

      FLUSH (UNIT_RES)

      RETURN
      END SUBROUTINE WRITE_RES0
