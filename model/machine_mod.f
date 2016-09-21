      MODULE machine

! Record length used in open statement for unformatted, direct access
! file, with 512 bytes per record.
      INTEGER, PARAMETER :: OPEN_N1 =   512

! Number words in 512 bytes
      INTEGER, PARAMETER :: NWORDS_DP =  64  ! double precision
      INTEGER, PARAMETER :: NWORDS_R  = 128  ! real
      INTEGER, PARAMETER :: NWORDS_I  = 128  ! integer

! computer node name/id
      CHARACTER(LEN=64) :: ID_NODE

! RUN ID info
      INTEGER :: ID_MONTH
      INTEGER :: ID_DAY
      INTEGER :: ID_YEAR
      INTEGER :: ID_HOUR
      INTEGER :: ID_MINUTE
      INTEGER :: ID_SECOND

    CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_RUN_ID                                             C
!  Author: P. Nicoletti                               Date: 16-DEC-91  C
!                                                                      C
!  Purpose: get the run id for this run                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_RUN_ID

      USE param
      IMPLICIT NONE

! temporary array to hold time data
      INTEGER DAT(8)
      CHARACTER(LEN=10) DATE, TIM, ZONE

      CALL DATE_AND_TIME(DATE, TIM, ZONE, DAT)
      ID_YEAR   = DAT(1)
      ID_MONTH  = DAT(2)
      ID_DAY    = DAT(3)
      ID_HOUR   = DAT(5)
      ID_MINUTE = DAT(6)
      ID_SECOND = DAT(7)

      CALL GET_ENVIRONMENT_VARIABLE('HOSTNAME',ID_NODE)

      RETURN
      END SUBROUTINE GET_RUN_ID

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function name: WALL_TIME (CPU)                                      C
!  Author: P. Nicoletti                               Date: 10-JAN-92  C
!                                                                      C
!  Purpose: returns wall time since start of the run                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION WALL_TIME()

      IMPLICIT NONE

      INTEGER, SAVE :: COUNT_OLD=0, WRAP=0

! local variables
!----------------------------------------------------------------------//
! Clock cycle
      INTEGER :: COUNT
! Number of cycles per second
      INTEGER :: COUNT_RATE
! Max number of cycles, after which count is reset to 0
      INTEGER :: COUNT_MAX

      CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
      IF(COUNT_OLD .GT. COUNT) THEN
        WRAP = WRAP + 1
      ENDIF
      COUNT_OLD = COUNT

      WALL_TIME      = DBLE(COUNT)/DBLE(COUNT_RATE) &
                     + DBLE(WRAP) * DBLE(COUNT_MAX)/DBLE(COUNT_RATE)
      END FUNCTION WALL_TIME


      END MODULE machine
