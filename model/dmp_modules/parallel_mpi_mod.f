module parallel_mpi

!       A module to carry out init, finalize and check for any parallel errors

   use geometry
   use compar
   implicit none

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: parallel_init                                          !
!  Author:                                            Date: XX-XXX-XX  !
!                                                                      !
!  Purpose: Wrapper function to initiailize OpenMP and MPI.            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine parallel_init()

      use compar, only: fbname
      USE funits, only: dmp_log, unit_log
      USE machine, only: get_run_id

      implicit none

      integer :: ierr

! Variables for generating file basename with processor id
      integer :: i1, i10, i100, i1000, i10000

      numPEs = 1
      myPE = 0

#ifdef MPI
      call MPI_Init(ierr)
      call MPI_Check( 'parallel_init:MPI_Init ', ierr)

      call MPI_COMM_SIZE( MPI_COMM_WORLD, numPEs, ierr )
      call MPI_Check( 'parallel_init:MPI_Comm_size ', ierr )

      call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
      call MPI_Check( 'parallel_init:MPI_Comm_size ', ierr )
#endif

! Only PE_IO to write out common error messages
      DMP_LOG = (myPE == PE_IO)

! PAR_I/O Generate file basename for LOG files
      i10000 = int(myPE/10000)
      i1000  = int((myPE-i10000*10000)/1000)
      i100   = int((myPE-i10000*10000-i1000*1000)/100)
      i10    = int((myPE-i10000*10000-i1000*1000-i100*100)/10)
      i1     = int((myPE-i10000*10000-i1000*1000-i100*100-i10*10)/1)

      i10000 = i10000 + 48
      i1000  = i1000  + 48
      i100   = i100   + 48
      i10    = i10    + 48
      i1     = i1     + 48

      fbname=char(i10000)//char(i1000)//char(i100)//char(i10)//char(i1)

 2000 format(/1x,'Number of SMP threads: ',I0,2/)

! Get the date and time. They give the unique run_id in binary output
! files
      CALL GET_RUN_ID

      return
      end subroutine parallel_init


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: parallel_fin                                           !
!  Author:                                            Date: XX-XXX-XX  !
!                                                                      !
!  Purpose: Wrapper for MPI_Finalize calls.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine parallel_fin()
      implicit none

#ifdef MPI
      integer :: ierr

      call MPI_Finalize(ierr)
      call MPI_Check( 'parallel_init:MPI_Finalize ', ierr)
#endif

      return
      end subroutine parallel_fin


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MPI_Check                                              !
!  Author:                                            Date: XX-XXX-XX  !
!                                                                      !
!  Purpose: Check error flags returned by MPI calls.                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine MPI_Check( msg, ierr )

#ifdef MPI
      use mpi, only: MPI_SUCCESS
#endif

      implicit none
      character(len=*),intent(in) :: msg
      integer, intent(in) :: ierr

#ifdef MPI
      character(len=512) :: errmsg
      integer :: resultlen, ierror

      if (ierr .ne. MPI_SUCCESS ) then
         call MPI_Error_string( ierr, errmsg, resultlen, ierror )
         print*, 'Error: ', msg
         print*, errmsg(1:resultlen)
         stop '** ERROR ** '
      endif
#endif

      return
      end subroutine MPI_Check
      end module parallel_mpi
