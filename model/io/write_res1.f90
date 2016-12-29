MODULE WRITE_RES1_MOD

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_RES1                                             C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!                                                                      C
!  Purpose: write out the time-dependent restart records               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_RES1(dt, nstep, time, &
         ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g)

      use compar  , only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      USE compar  , only: mype, pe_io
      USE funits  , only: unit_res
      USE geometry, only: ijkmax2, ijkmax3
      USE out_bin_512_mod, only: out_bin_512
      USE in_binary_512, only: in_bin_512, convert_to_io_dp

      IMPLICIT NONE

      double precision, intent(in   ) :: dt, time
      integer, intent(in   ) :: nstep

      double precision, intent(in   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

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

      call gatherWriteRes (ep_g,array2, array1, NEXT_REC)
      call gatherWriteRes (p_g,array2, array1, NEXT_REC)

      call gatherWriteRes (ro_g,array2, array1, NEXT_REC)
      call gatherWriteRes (rop_g,array2, array1, NEXT_REC)

      call gatherWriteRes (u_g,array2, array1, NEXT_REC)
      call gatherWriteRes (v_g,array2, array1, NEXT_REC)
      call gatherWriteRes (w_g,array2, array1, NEXT_REC)

      if(myPE.eq.PE_IO) FLUSH(UNIT_RES)

      deallocate (array1)
      deallocate (array2)

      end subroutine write_res1

!---------------------------------------------------------------------

      subroutine gatherWriteRes(VAR, array2, array1, NEXT_REC)

      use param   , only: dimension_3
      USE compar  , only: mype, pe_io
      USE funits  , only: unit_res
      USE geometry, only: ijkmax2, ijkmax3

      USE out_bin_512_mod, only: out_bin_512
      USE in_binary_512  , only: convert_to_io_dp

      implicit none

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
END MODULE WRITE_RES1_MOD
