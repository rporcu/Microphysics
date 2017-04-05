module exit_mod

   contains

      subroutine mfix_exit(myID)

      use compar, only: myPE, pe_io

      implicit none

! Rank ID
      integer, INTENT(IN) :: myID
! Logical showing that a file unit is open.
! Process ID (myPE) converted to a character string.
      CHARACTER(len=64) :: myID_c

! Set the ID of the caller.
      myID_c=''; WRITE(myID_c,*) myID

! Write out that this routine was called.
      IF(myPE == PE_IO) WRITE(*,1000)

! Terminate MPI.

! Close any open files.

! Last gasp...
      IF(myPE == PE_IO) WRITE(*,1002)

! Hard Stop.
      STOP 1

 1000 FORMAT(2/,1x,70('*'),/' Fatal error reported on one or more',    &
        ' processes. The .LOG file',/' may contain additional',        &
        ' information about the failure.',/1x,70('*'))


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
      integer, INTENT(IN) :: UNIT_l

! Local Variables.
!---------------------------------------------------------------------//
! Retruned status of the specifed file unit
      integer :: IOS
! Logical indicating if the file is open
      logical :: FOPEN

! If the file is open...
      INQUIRE(UNIT=UNIT_l, OPENED=FOPEN, IOSTAT=IOS )
! Close it.
      IF(FOPEN) CLOSE(UNIT_l)

      end subroutine close_file
end module exit_mod
