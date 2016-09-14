!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_RES1                                             C
!  Purpose: write out the time-dependent restart records               C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_RES1
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar           !//
      USE fldvar
      USE funits
      USE geometry
      USE machine, only: flush_res
      USE mpi_utility      !//d pnicol : for gather
      USE output
      USE param
      USE param1
      USE physprop
      USE run
      USE sendrecv         !//d pnicol : for gather
!//d pnicol  ... not needed    USE tmp_array
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
!
!
!//d pnicol : allocate arrays for gather/convert_to_io_dp
      double precision, allocatable :: array1(:)
      double precision, allocatable :: array2(:)

!             pointer to first time-dependent record in restart file
      INTEGER :: NEXT_REC
!-----------------------------------------------
!

!//d pnicol
!      if (myPE.eq.PE_IO .and. .not.distio) then
         allocate (array1(ijkmax2))
         allocate (array2(ijkmax3))
!      else
!         allocate (array1(1))
!         allocate (array2(1))
!      end if


      if (myPE.eq.PE_IO) then
         READ (UNIT_RES, REC=3) NEXT_REC
         WRITE (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP
         NEXT_REC = NEXT_REC + 1
      end if
!
!\\SP Local Send Receive - need to be moved to source later!!
      call send_recv(EP_g,2)
      call send_recv(P_g,2)
      call send_recv(RO_g,2)
      call send_recv(ROP_g,2)
      call send_recv(U_g,2)
      call send_recv(V_g,2)
      call send_recv(W_g,2)


      call gatherWriteRes (EP_g,array2, array1, NEXT_REC)  !//d pnicol
      call gatherWriteRes (P_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (RO_g,array2, array1, NEXT_REC)  !//d pnicol
      call gatherWriteRes (ROP_g,array2, array1, NEXT_REC)  !//d pnicol
!
      call gatherWriteRes (U_g,array2, array1, NEXT_REC)  !//d pnicol
      call gatherWriteRes (V_g,array2, array1, NEXT_REC)  !//d pnicol
      call gatherWriteRes (W_g,array2, array1, NEXT_REC)  !//d pnicol
!
!
!     Version 1.6
!---------------------------------------------------------------------

      if(myPE.eq.PE_IO)CALL FLUSH_res (UNIT_RES)

!      call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
      deallocate (array1)  !//d pnicol
      deallocate (array2)  !//d pnicol
!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here

!
      RETURN
      END SUBROUTINE WRITE_RES1

      subroutine gatherWriteRes(VAR, array2, array1, NEXT_REC)

      USE geometry
      USE funits
      USE compar           !//
      USE mpi_utility      !//d pnicol : for gather
      USE sendrecv         !//d pnicol : for gather

      USE cutcell
      USE in_binary_512

      IMPLICIT NONE

      double precision, dimension(ijkmax2) :: array1
      double precision, dimension(ijkmax3) :: array2
      double precision, dimension(DIMENSION_3) :: VAR

      INTEGER :: NEXT_REC

!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here

      CALL gather (VAR,array2,root)

      if (myPE.eq.PE_IO) then
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC)
      end if


      End subroutine gatherWriteRes
