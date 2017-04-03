module exit_mod

   contains

      subroutine mfix_exit(myID)

! File unit for .OUT file
      use funits, only : UNIT_OUT
! File unit for .LOG files
      use funits, only : UNIT_LOG

      use compar, only: myPE, pe_io
      use funits, only: dmp_log

      implicit none

! Rank ID
      INTEGER, INTENT(IN) :: myID
! Logical showing that a file unit is open.
      LOGICAL :: isOpen
! Process ID (myPE) converted to a character string.
      CHARACTER(len=64) :: myID_c

! Set the ID of the caller.
      myID_c=''; WRITE(myID_c,*) myID

! Write out that this routine was called.
      IF(myPE == PE_IO) WRITE(*,1000)
      IF(DMP_LOG) THEN
         INQUIRE(UNIT=UNIT_LOG,OPENED=isOpen)
         IF(isOPEN) WRITE(UNIT_LOG,1001) trim(adjustl(myID_c))
      ENDIF

! Terminate MPI.
      ! CALL exitMPI(myID_l)

! Close any open files.
      CALL CLOSE_FILE(UNIT_OUT)
      CALL CLOSE_FILE(UNIT_LOG)

! Last gasp...
      IF(myPE == PE_IO) WRITE(*,1002)

! Hard Stop.
      STOP 1

 1000 FORMAT(2/,1x,70('*'),/' Fatal error reported on one or more',    &
        ' processes. The .LOG file',/' may contain additional',        &
        ' information about the failure.',/1x,70('*'))

 1001 FORMAT(2/,1x,70('*'),/' Fatal error reported on PE ',  &
         A,'. The .LOG file may contain',/' additional ',     &
        'information about the failure.',/1x,70('*'))

 1002 FORMAT(2/,1x,'Program Terminated.',2/)

      end subroutine mfix_exit

!``````````````````````````````````````````````````````````````````````!
! Subroutine: close_file                                               !
!                                                                      !
! Purpose: Close a file if it is open.                                 !
!......................................................................!
      subroutine close_file(unit_l)

! Global Variables.
!---------------------------------------------------------------------//
! NONE

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: UNIT_l

! Local Variables.
!---------------------------------------------------------------------//
! Retruned status of the specifed file unit
      INTEGER :: IOS
! Logical indicating if the file is open
      LOGICAL :: FOPEN

! If the file is open...
      INQUIRE(UNIT=UNIT_l, OPENED=FOPEN, IOSTAT=IOS )
! Close it.
      IF(FOPEN) CLOSE(UNIT_l)

      end subroutine close_file
end module exit_mod
