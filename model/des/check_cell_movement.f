!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT                                     !
!                                                                      !
!  Purpose: Check to see if particles have moved into ghost cells.     !
!                                                                      !
!  Note: This routine needs a global communicator to identify errors.  !
!  The collection could get expensive so the call frequency of this    !
!  routine should probably be reduced.                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT

! Global Variables:
!---------------------------------------------------------------------//
! Max number of particles in process
      use discretelement, only: MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! Run time flag indicating DEM or PIC solids.
      use run, only: DEM_SOLIDS

      use functions, only: IS_NORMAL
      use mpi_utility
      use error_manager

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! Loop indicies:
      INTEGER :: L, I, J, K
! Integer error flag.
      INTEGER :: IER

! Initialize local variables.
      IER = 0

! Set an error flag if any errors are found. Preform a global collection
! to sync error flags. If needed, reort errors.
!.......................................................................
!!$omp parallel default(shared) private(L, I, J, K, IJK)
!!$omp do reduction(+:IER) schedule (guided,50)
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.IS_NORMAL(L)) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF(I > IEND1 .OR. I < ISTART1) IER = 1
         IF(J > JEND1 .OR. J < JSTART1) IER = 1
         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) IER = 1
      ENDDO
!!$omp end parallel

      CALL GLOBAL_ALL_SUM(IER)
      IF(IER == 0) RETURN

      IF(DEM_SOLIDS) CALL CHECK_CELL_MOVEMENT_DEM

      RETURN
      END SUBROUTINE CHECK_CELL_MOVEMENT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT_DEM                                 !
!                                                                      !
!  Purpose: Report which DEM particles have moved into ghost cells.    !
!  This is a dead-end routine. Once called, the simulation will exit.  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT_DEM

! Global Variables:
!---------------------------------------------------------------------//
! The global ID of a particle.
      use discretelement, only: iGlobal_ID
! Max number of particles in process
      use discretelement, only: MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! Particle positions and velocities
      use discretelement, only: DES_POS_NEW, DES_VEL_NEW

      use functions, only: IS_NORMAL
      use mpi_utility
      USE error_manager

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! Loop indicies:.
      INTEGER :: L, I, J, K
! Integer error flag
      INTEGER :: IER


      CALL INIT_ERR_MSG("CHECK_CELL_MOVEMENT_DEM")
      CALL OPEN_PE_LOG(IER)

      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1100 FORMAT('Error 1100: Particles detected in a ghost cell:',/' ')

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.IS_NORMAL(L)) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_ID(L))),'I',        &
               trim(iVal(I)),'X',DES_POS_NEW(1,L),'X',DES_VEL_NEW(1,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_id(L))),'J',        &
               trim(iVal(J)),'Y',DES_POS_NEW(2,L),'Y',DES_VEL_NEW(2,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF (DO_K .AND. (K.GT.KEND1 .OR. K.LT.KSTART1)) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_ID(L))),'K',        &
               trim(iVal(K)),'Z',DES_POS_NEW(3,L),'Z',DES_VEL_NEW(3,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

 1101 FORMAT('Particle ',A,' moved into cell with ',A,' index ',A,/ &
         3x,A,'-Position: ',g11.4,6x,A,'-Velocity:',g11.4,/' ')

      WRITE(ERR_MSG, 1102)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
 1102 FORMAT('This is a fatal error. A particle output file (vtp) ',   &
         'will be written',/'to aid debugging.')


      CALL WRITE_DES_DATA
      CALL MFIX_EXIT(myPE)

      END SUBROUTINE CHECK_CELL_MOVEMENT_DEM
