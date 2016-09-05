!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SPX1                                             C
!  Purpose: write out the time-dependent restart records (REAL)        C
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
!  Variables referenced: TIME, NSTEP, EP_g, P_g, U_g                   C
!                        IJKMAX2, MMAX                                 C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables:  LC, N, NEXT_REC, NUM_REC                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_SPX1(L, unit_add)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE cutcell
      USE fldvar
      USE funits
      USE geometry
      USE machine
      USE mpi_utility
      USE output
      USE param
      USE param1
      USE physprop
      USE run
      use discretelement, only: PRINT_DES_DATA
      use discretelement, only: DISCRETE_ELEMENT
      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: PARTICLES

!//       USE tmp_array
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!             flag whether to write a particular SPx file
      INTEGER L

!              offset for use in post_mfix
      INTEGER  unit_add
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
! local variables
!
!//
      double precision, allocatable :: array1(:)     !//
      double precision, allocatable :: array2(:)     !//

!             loop counters
      INTEGER LC, N
!
!             Pointer to the next record
      INTEGER NEXT_REC
!
!              Number of records written each time step
      INTEGER  NUM_REC

      INTEGER  uspx   ! UNIT_SPX + offset from post_mfix
      CHARACTER(LEN=50), DIMENSION(1) :: LINE   !error message
      double precision, dimension(:), allocatable :: TMP_VAR

      allocate(TMP_VAR(DIMENSION_3))

!-----------------------------------------------
      uspx = UNIT_SPX + unit_add

!
      if (myPE .eq.PE_IO) then
         allocate (array1(ijkmax2))   !//
         allocate (array2(ijkmax3))   !//
      else
         allocate (array1(1))   !//
         allocate (array2(1))   !//
      end if
!
! ".SP1" FILE         EP_g    [ ROP_g, RO_g  must be calculated ...
!                                        not written out ]
!
      SELECT CASE (L)
      CASE (1)
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if

! Explicitly coupled simulations do not need to rebin particles to
! the fluid grid every time step. However, this implies that the
! fluid cell information and interpolation weights become stale.
         IF(DISCRETE_ELEMENT .AND. .NOT.DES_CONTINUUM_COUPLED) THEN
! Bin particles to fluid grid.
            CALL PARTICLES_IN_CELL
! Calculate interpolation weights
            CALL CALC_INTERP_WEIGHTS
! Calculate mean fields (EPg).
            CALL COMP_MEAN_FIELDS
         ENDIF


         call gatherWriteSpx (EP_g,array2, array1, uspx+L, NEXT_REC)

         if (myPE .eq. PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if


! remove this redundant call here to write_des_data in case of new
! coupled runs that contain at least a particle and involve some initial
! settling of the system.
! the call made in des_time_march is a better call for capturing the
! initial state of such a des continuum coupled system
         IF(DISCRETE_ELEMENT.AND.PRINT_DES_DATA .AND. &
            .NOT.(TRIM(RUN_TYPE)=='NEW' .AND. PARTICLES /=0 .AND. &
            TIME == ZERO)) THEN
               CALL WRITE_DES_DATA
         ENDIF



!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
! ".SP2" FILE         P_g 
!
      CASE (2)
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         call gatherWriteSpx (P_g,array2, array1, uspx+L, NEXT_REC)

         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if


! ".SP3" FILE         U_g , V_g , W_g
!
      CASE (3)
         if (myPE.eq.PE_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if
         call gatherWriteSpx (U_g,array2, array1, uspx+L, NEXT_REC)   !//
         call gatherWriteSpx (V_g,array2, array1, uspx+L, NEXT_REC)   !//
         call gatherWriteSpx (W_g,array2, array1, uspx+L, NEXT_REC)   !//
         if (myPE.eq.PE_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
!
      CASE (4)
      CASE (5)
      CASE (6)
      CASE (7)
      CASE (8)
      CASE (9)
      CASE (10)  ! Reaction rates
      CASE (11)
!
!
      CASE DEFAULT
            LINE(1) = 'Unknown SPx file index'
            CALL WRITE_ERROR ('WRITE_SPX1', LINE, 1)
            CALL MFIX_EXIT(myPE)
      END SELECT

!//      call unlock_tmp_array
!
      deallocate (array1)
      deallocate (array2)
      deallocate (TMP_VAR)
!
      RETURN
      END SUBROUTINE WRITE_SPX1

      subroutine gatherWriteSpx(VAR, array2, array1, uspxL, NEXT_REC)
        USE geometry
        USE compar           !//
        USE mpi_utility      !//d pnicol : for gatherWriteSpx
        USE sendrecv         !//d pnicol : for gatherWriteSpx
        USE cutcell
        USE in_binary_512
        IMPLICIT NONE
        integer uspxL, NEXT_REC
        double precision, dimension(ijkmax2) :: array1
        double precision, dimension(ijkmax3) :: array2
        double precision, dimension(DIMENSION_3) :: VAR
        double precision, dimension(:), allocatable :: TMP_VAR

        allocate(TMP_VAR(DIMENSION_3))

!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
        IF(RE_INDEXING) THEN
           TMP_VAR = UNDEFINED
           CALL gather (TMP_VAR,array2,root)
        ELSE
           CALL gather (VAR,array2,root)
        ENDIF
!        call gather (VAR,array2,root)  !//d pnicol

!       call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here
        if (myPE.eq.PE_IO) then
           call convert_to_io_dp(array2,array1,ijkmax2)
           CALL OUT_BIN_R (uspxL, array1, IJKMAX2, NEXT_REC)
        end if

        deallocate(TMP_VAR)

      End subroutine gatherWriteSpx






