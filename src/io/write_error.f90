MODULE WRITE_ERROR_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Write_error(Name, Line, L)                             C                     C
!  Purpose: Write an error message                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_ERROR(NAME, LINE, LMAX)

      use funits, only: unit_log, dmp_log

      IMPLICIT NONE

!                      Subroutine name
      CHARACTER(LEN=*)    Name

!                      Message
      CHARACTER(LEN=*)    LINE(*)

!                      Dimension of message array
      INTEGER          LMAX

!                      Index
      INTEGER          L

!-----------------------------------------------

      IF(DMP_LOG)WRITE (UNIT_LOG, 1000) NAME
      DO L = 1, LMAX
         IF(DMP_LOG)WRITE (UNIT_LOG, 1010) LINE(L)
      END DO
      IF(DMP_LOG)WRITE (UNIT_LOG, 1020)
      RETURN
 1000 FORMAT(1X,70('*'),/,/,1X,'From : ',A)
 1010 FORMAT(1X,A)
 1020 FORMAT(/,/,1X,70('*'))
      END SUBROUTINE WRITE_ERROR
END MODULE WRITE_ERROR_MODULE
