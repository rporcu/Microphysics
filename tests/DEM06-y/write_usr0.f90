!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR0                                             C
!  Purpose: Write initial part of user-defined output                  C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine write_usr0() &
        bind(C, name="write_usr0")

      implicit none

      CALL WRITE_DAT_HEADER('POST_POS.dat','Pos')
      CALL WRITE_DAT_HEADER('POST_VEL.dat','Vel')

      contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_DAT_HEADER(FNAME, VAR)

      use discretelement, only: DES_ONEWAY_COUPLED
      use discretelement, only: DES_EXPLICITLY_COUPLED

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
      ELSE
         OPEN(UNIT=fUNIT,FILE=FNAME,POSITION="APPEND",STATUS='OLD')
      ENDIF

      WRITE(fUNIT, 1100) DES_ONEWAY_COUPLED

      WRITE(fUNIT, 1140, ADVANCE='YES') DES_EXPLICITLY_COUPLED

      WRITE(fUNIT, 1250) VAR, VAR

 1000 FORMAT(2/,25x,A)

 1100 FORMAT(2/,7x,'DES_ONEWAY_COUPLED =',6x,L1)

 1140 FORMAT(7x,'DES_EXPLICITLY_COUPLED =',2x,L1)

 1250 FORMAT(/10X,'Time',17X,A,13X,A,'_MFIX',9X,'%REL DIFF')

      CLOSE(fUNIT)
      RETURN
      END SUBROUTINE WRITE_DAT_HEADER

      END SUBROUTINE WRITE_USR0
