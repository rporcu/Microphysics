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

      call write_dat_header('POST_POS.dat','Pos')
      call write_dat_header('POST_VEL.dat','Vel')

      CONTAINS

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      subroutine write_dat_header(FNAME, VAR)

      use run, only: DESCRIPTION

      use discretelement, only: KN, KN_W
      use discretelement, only: DES_EN_WALL_INPUT

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

      WRITE(fUNIT, 1110) KN, KN_W
      WRITE(fUNIT, 1120) DES_EN_WALL_INPUT(1), DES_EN_WALL_INPUT(1)

      WRITE(fUNIT, 1200) VAR, VAR

 1000 FORMAT(2/,25x,A)

 1100 FORMAT(2/,7x,'Time Stepping Scheme: ',A)

 1110 FORMAT(/7x,'Normal collision spring constant. (N/m)',/&
         10x,'KN = ',T30,G12.4,/&
         10x,'KN_W = ',T30,G12.4)

 1120 FORMAT(/7x,'Restitution coefficient. (1)',/&
         10x,'DES_EN_INPUT = ',T30,G12.4,/&
         10x,'DES_EN_WALL_INPUT = ',T30,G12.4)

 1200 FORMAT(/7X,'Time',7x,'Stage',11X,A,13X,A,'_MFIX')

      CLOSE(fUNIT)

      end subroutine write_dat_header

      end subroutine write_usr0
