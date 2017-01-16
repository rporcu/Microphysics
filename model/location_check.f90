MODULE LOCATION_CHECK_MODULE
   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: LOCATION_CHECK                                          !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!                                                                      !
!  Purpose: Check calculated and given cell locations for consistency  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LOCATION_CHECK(CELL_SPECIFIED, CELL_CALCULATED,       &
         COUNTER, MESSAGE)

      ! Module procedure for error message management.
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      ! Cell index specified in the input file.
      INTEGER, INTENT(IN) :: CELL_SPECIFIED

      ! Cell index calculated for location coordinate.
      INTEGER, INTENT(IN) :: CELL_CALCULATED

      ! Index for BC, IC, or IS
      INTEGER, INTENT(IN) :: COUNTER

      ! Error message to print out
      CHARACTER(len=*) :: MESSAGE
!......................................................................!

! Check that the cell_specified in the data input equals to the cell
! calculated.

      IF(CELL_SPECIFIED == CELL_CALCULATED) RETURN

      IF(MESSAGE(6:6)=='b' .OR. MESSAGE(6:6)=='t') RETURN
      IF(MESSAGE(6:6)=='s' .OR. MESSAGE(6:6)=='n') RETURN
      IF(MESSAGE(6:6)=='w' .OR. MESSAGE(6:6)=='e') RETURN

      CALL INIT_ERR_MSG('LOCATION_CHECK')

      WRITE(ERR_MSG, 1000) MESSAGE, COUNTER, CELL_SPECIFIED,           &
         CELL_CALCULATED
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1000 FORMAT('Error 1000: IC, BC, OR IS consistency error for: ',A,/,  &
      'IC/BC/IS No',5X,'= ',I6,/,'Cell specified',2x,'= ',I6,/,        &
      'Cell calculated',1x,'= ',I6)

      CALL FINL_ERR_MSG()

      RETURN
      END SUBROUTINE LOCATION_CHECK
END MODULE LOCATION_CHECK_MODULE
