!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USR0                                             !
!  Purpose: Write initial part of user-defined output                  !
!                                                                      !
!  Author:                                            Date: dd-mmm-yy  !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine write_usr0() &
        bind(C, name="write_usr0")

      implicit none

      CALL WRITE_DAT_HEADER('POST_TIME.dat','TIME')
      CALL WRITE_DAT_HEADER('POST_TVEL.dat','TVEL')
      CALL WRITE_DAT_HEADER('POST_AVEL.dat','AVEL')

      contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_DAT_HEADER(FNAME, VAR)

      use run, only: DESCRIPTION


      IMPLICIT NONE

      CHARACTER(len=*) :: FNAME
      CHARACTER(len=*) :: VAR

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: fUNIT = 2030

      INQUIRE(FILE=FNAME,EXIST=EXISTS)
      IF (.NOT.EXISTS) THEN
         OPEN(UNIT=fUNIT,FILE=FNAME,STATUS='NEW')
         WRITE(fUNIT, 1000) trim(DESCRIPTION)

         WRITE(fUNIT, 1300) 'MEW', VAR, VAR
      ENDIF


 1000 FORMAT(2/,25x,A)

 1100 FORMAT(2/,7x,'Time Stepping Scheme: ',A)

 1300 FORMAT(2/11X,A,17X,A,12X,A,'_MFIX',9X,'%REL ERR')

      CLOSE(fUNIT)
      RETURN
      END SUBROUTINE WRITE_DAT_HEADER

      END SUBROUTINE WRITE_USR0
