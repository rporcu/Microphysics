!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_RES1                                              C
!  Purpose: read in the time-dependent restart records                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Local variables: TIME_READ, LC, NEXT_REC                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE READ_RES1
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE physprop
      USE run
      USE funits
      USE compar
      USE in_binary_512

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
!             pointer to the next record
      INTEGER NEXT_REC
!
!                file version id
      CHARACTER(LEN=512) :: VERSION
!
!                version number
      REAL       VERSION_NUMBER
!
!                      Dummy array
      DOUBLE PRECISION DT_SAVE

!//PAR_I/O declare global scratch arrays
      double precision, allocatable :: array1(:)
      double precision, allocatable :: array2(:)
!-----------------------------------------------
!
      allocate (array1(ijkmax2))
      allocate (array2(ijkmax3))

      IF (DT_FAC == ONE) DT_SAVE = DT

      if (myPE == PE_IO) then
         READ (UNIT_RES, REC=1) VERSION
         READ (VERSION(6:512), *) VERSION_NUMBER

         READ (UNIT_RES, REC=3) NEXT_REC
         IF (VERSION_NUMBER >= 1.12) THEN
            READ (UNIT_RES, REC=NEXT_REC) TIME, DT, NSTEP
         ELSE
            READ (UNIT_RES, REC=NEXT_REC) TIME, NSTEP
         ENDIF
         NEXT_REC = NEXT_REC + 1
      end if


      !call bcast(VERSION, PE_IO)        !//PAR_I/O BCAST0c
      !call bcast(VERSION_NUMBER, PE_IO) !//PAR_I/O BCAST0r
      !call bcast(TIME, PE_IO)           !//PAR_I/O BCAST0d
      !call bcast(NSTEP, PE_IO)          !//PAR_I/O BCAST0i
      ! if (VERSION_NUMBER >= 1.12) !call bcast(DT, PE_IO)   !//PAR_I/O BCAST0d

!
      call readScatterRes(EP_G,array2, array1, 0, NEXT_REC)
      call readScatterRes(P_G,array2, array1, 0, NEXT_REC)
      call readScatterRes(RO_G,array2, array1, 0, NEXT_REC)
      call readScatterRes(ROP_G,array2, array1, 0, NEXT_REC)
!
      call readScatterRes(U_G, array2, array1, 0, NEXT_REC)
      call readScatterRes(V_G, array2, array1, 0, NEXT_REC)
      call readScatterRes(W_G, array2, array1, 0, NEXT_REC)
!

!------------------------------------------------------------------------
!

!      call MPI_barrier(MPI_COMM_WORLD,mpierr)
      deallocate( array1 )
      deallocate( array2 )
!      call MPI_barrier(MPI_COMM_WORLD,mpierr)

      ! call send_recv(rop_g)
      ! call send_recv(ro_g)


      IF (DT_FAC == ONE) DT = DT_SAVE
!

      RETURN
      END SUBROUTINE READ_RES1



      subroutine readScatterRes(VAR, array2, array1, init, NEXT_REC)

      use param1, only: zero, undefined
      use param, only: dimension_3
      USE geometry
      USE funits
      USE compar
      USE in_binary_512
      IMPLICIT NONE
      double precision, dimension(ijkmax2) :: array1
      double precision, dimension(ijkmax3) :: array2
      double precision, dimension(DIMENSION_3) :: VAR
      INTEGER :: init  ! define VAR initialization, 0: undefin, 1: zero
      INTEGER :: NEXT_REC

      if( init==0 ) then
        array1(:) = Undefined
        array2(:) = Undefined
      else
        array1(:) = zero
        array2(:) = zero
      endif

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, array1, IJKMAX2, NEXT_REC)
         CALL convert_from_io_dp(array1, array2, IJKMAX2)
      end if
      var = array2
      ! call scatter(VAR, array2, PE_IO)

      End subroutine readScatterRes
