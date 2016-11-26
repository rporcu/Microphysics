!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: WRITE_RES0                                              C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!                                                                      C
!  Purpose: write out the initial restart records (namelist data)      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE WRITE_RES0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE funits
      USE geometry
      USE ic
      USE in_binary_512i
      USE leqsol
      USE machine
      USE output
      USE param
      USE param1
      USE physprop
      USE run
      USE scales
      USE toleranc
      USE ur_facs
      USE fldvar

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
