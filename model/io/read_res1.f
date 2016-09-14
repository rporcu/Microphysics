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
      USE mpi_utility
      USE sendrecv
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
      if (myPE .eq. PE_IO) then
         allocate (array1(ijkmax2))
         allocate (array2(ijkmax3))
      else
         allocate (array1(1))
         allocate (array2(1))
      end if

!      call MPI_barrier(MPI_COMM_WORLD,mpierr)
!
!     Use DT from data file if DT_FAC is set to 1.0
      IF (DT_FAC == ONE) DT_SAVE = DT
!

!//PAR_I/O only PE_IO reads the restart file
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


      call bcast(VERSION, PE_IO)        !//PAR_I/O BCAST0c
      call bcast(VERSION_NUMBER, PE_IO) !//PAR_I/O BCAST0r
      call bcast(TIME, PE_IO)           !//PAR_I/O BCAST0d
      call bcast(NSTEP, PE_IO)          !//PAR_I/O BCAST0i
      if (VERSION_NUMBER >= 1.12) call bcast(DT, PE_IO)   !//PAR_I/O BCAST0d

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

      call send_recv(rop_g)
      call send_recv(ro_g)


      IF (DT_FAC == ONE) DT = DT_SAVE
!

!     We may no longer need PATCH_AFTER_RESTART
!     CALL PATCH_AFTER_RESTART

      RETURN
      END SUBROUTINE READ_RES1



      subroutine readScatterRes(VAR, array2, array1, init, NEXT_REC)

      use param1, only: zero, undefined
      use param, only: dimension_3
      USE geometry
      USE funits
      USE compar
      USE mpi_utility
      USE sendrecv
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
      call scatter(VAR, array2, PE_IO)

      End subroutine readScatterRes




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PATCH_AFTER_RESTART                                    C
!  Purpose: Patch new fluid cells after a restart                      C
!           This could occur when restarting with a different          C
!           grid partition when the RESTART file was generated         C
!           prior to the Dec. 4th 2014 bug fix.                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 14-APR-15  C
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

      SUBROUTINE PATCH_AFTER_RESTART
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
      USE mpi_utility
      USE sendrecv
      USE cutcell
      use functions

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
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I,J,K, IJK, IJKNB
      INTEGER :: NB
      INTEGER, DIMENSION(6) :: NBCELL
      LOGICAL :: NB_FOUND


!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK).AND.EP_G(IJK)==UNDEFINED) THEN

! Detects new fluid cells that used to be blocked cells with undefined
! values. When a fluid cell has undefined void fraction, this means all
! variables need to be patched. Typically, this will be a fairly small
! cut cell that was flagged as BLOCKED cell with a different partition.
! If a valid fluid cell is found next to this undefined cell, all field
! variables will be copied over. If no valid fluid cell is found, the
! code will continue and will likely stop during the check_data_30
! (zero species mass fractions will yield a zero specific heat).
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            NBCELL(1) = IM_OF(IJK)
            NBCELL(2) = JM_OF(IJK)
            NBCELL(3) = KM_OF(IJK)
            NBCELL(4) = IP_OF(IJK)
            NBCELL(5) = JP_OF(IJK)
            NBCELL(6) = KP_OF(IJK)

            NB_FOUND = .FALSE.

            DO NB = 1,6

               IJKNB = NBCELL(NB)

               IF(FLUID_AT(IJKNB).AND.EP_G(IJKNB)/=UNDEFINED) THEN
                  NB_FOUND = .TRUE.
                  WRITE (*, 1010) MyPE, I,J,K

                  EP_G(IJK)   = EP_G(IJKNB)
                  P_G(IJK)    = P_G(IJKNB)
                  RO_G(IJK)   = RO_G(IJKNB)
                  ROP_G(IJK)  = ROP_G(IJKNB)

                  U_G(IJK) = U_G(IJKNB)
                  V_G(IJK) = V_G(IJKNB)
                  W_G(IJK) = W_G(IJKNB)


                  EXIT ! Exit as soon as first valid neighbor cell is found
               ENDIF  ! NB is a fluid cell

            ENDDO ! NB Loop

            IF(.NOT.NB_FOUND) WRITE (*, 1020) MyPE, I,J,K   ! NO FLUID CELL AMONG NEIGBHORS

         ENDIF ! New fuid cell

      ENDDO ! IJK loop

1010  FORMAT(1X,'PATCHING NEW FLUID CELL UPON RESTART: MyPE,I,J,K =' ,I6,I6,I6,I6)
1020  FORMAT(1X,'UNABLE TO PATCH NEW FLUID CELL UPON RESTART: MyPE,I,J,K =' ,I6,I6,I6,I6)
      END SUBROUTINE PATCH_AFTER_RESTART
