!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_RES1                                             C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!                                                                      C
!  Purpose: write out the time-dependent restart records               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_RES1
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE fldvar
      USE funits
      USE geometry
      USE output
      USE param
      USE param1
      USE physprop
      USE run

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
! Allocatable arrays for gather/convert_to_io_dp
      double precision, allocatable :: array1(:)
      double precision, allocatable :: array2(:)

! Pointer to first time-dependent record in restart file
      INTEGER :: NEXT_REC
!-----------------------------------------------
!
      allocate (array1(ijkmax2))
      allocate (array2(ijkmax3))

      if (myPE.eq.PE_IO) then
         READ (UNIT_RES, REC=3) NEXT_REC
         WRITE (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP
         NEXT_REC = NEXT_REC + 1
      end if

! Local Send Receive - need to be moved to source later!!
      ! call send_recv(EP_g,2)
      ! call send_recv(P_g,2)
      ! call send_recv(RO_g,2)
      ! call send_recv(ROP_g,2)
      ! call send_recv(U_g,2)
      ! call send_recv(V_g,2)
      ! call send_recv(W_g,2)

      call gatherWriteRes (EP_g,array2, array1, NEXT_REC)
      call gatherWriteRes (P_g,array2, array1, NEXT_REC)

      call gatherWriteRes (RO_g,array2, array1, NEXT_REC)
      call gatherWriteRes (ROP_g,array2, array1, NEXT_REC)

      call gatherWriteRes (U_g,array2, array1, NEXT_REC)
      call gatherWriteRes (V_g,array2, array1, NEXT_REC)
      call gatherWriteRes (W_g,array2, array1, NEXT_REC)

!     Version 1.6
!---------------------------------------------------------------------

      if(myPE.eq.PE_IO) FLUSH(UNIT_RES)

      deallocate (array1)
      deallocate (array2)

      RETURN
      END SUBROUTINE WRITE_RES1

      subroutine gatherWriteRes(VAR, array2, array1, NEXT_REC)

      USE geometry
      USE funits
      USE compar
      USE param

      USE in_binary_512

      IMPLICIT NONE

      double precision, dimension(ijkmax2) :: array1
      double precision, dimension(ijkmax3) :: array2
      double precision, dimension(DIMENSION_3) :: VAR

      INTEGER :: NEXT_REC

!     call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here

      array2 = var
      ! CALL gather (VAR,array2,root)

      if (myPE.eq.PE_IO) then
         call convert_to_io_dp(array2,array1,ijkmax2)
         CALL OUT_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC)
      end if


      End subroutine gatherWriteRes
