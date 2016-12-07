MODULE CHECK_CELL_MOVEMENT_MODULE
CONTAINS

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
      SUBROUTINE CHECK_CELL_MOVEMENT(pijk, iglobal_id, des_vel_new, des_pos_new)

! Global Variables:
!---------------------------------------------------------------------//
! Max number of particles in process
      use discretelement, only: MAX_PIP
! Run time flag indicating DEM or PIC solids.
      use run, only: DEM_SOLIDS

      use compar, only: mype
      use compar, only: istart1, jstart1, kstart1
      use compar, only: iend1, jend1, kend1
      use exit_mod, only: mfix_exit

      use discretelement, only: normal_particle, particle_state
      use error_manager, only: err_msg, flush_err_msg, ival, init_err_msg

      IMPLICIT NONE

      INTEGER, DIMENSION(:), INTENT(OUT) :: iglobal_id
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new


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
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(L)) CYCLE

         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(L)) CYCLE
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF(I > IEND1 .OR. I < ISTART1) IER = 1
         IF(J > JEND1 .OR. J < JSTART1) IER = 1
         IF(K > KEND1 .OR. K < KSTART1) IER = 1
      ENDDO

      ! CALL GLOBAL_ALL_SUM(IER)
      IF(IER == 0) RETURN

      IF(DEM_SOLIDS) CALL CHECK_CELL_MOVEMENT_DEM

      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT_DEM                                 !
!                                                                      !
!  Purpose: Report which DEM particles have moved into ghost cells.    !
!  This is a dead-end routine. Once called, the simulation will exit.  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT_DEM

! Max number of particles in process
      use discretelement, only: MAX_PIP

      USE open_files_mod, only: open_pe_log

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
         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(L)) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(L)) CYCLE
         K = PIJK(L,3)

         IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_ID(L))),'I',        &
               trim(iVal(I)),'X',DES_POS_NEW(L,1),'X',DES_VEL_NEW(L,1)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_id(L))),'J',        &
               trim(iVal(J)),'Y',DES_POS_NEW(L,2),'Y',DES_VEL_NEW(L,2)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF (K.GT.KEND1 .OR. K.LT.KSTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_ID(L))),'K',        &
               trim(iVal(K)),'Z',DES_POS_NEW(L,3),'Z',DES_VEL_NEW(L,3)
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

      END SUBROUTINE CHECK_CELL_MOVEMENT
END MODULE CHECK_CELL_MOVEMENT_MODULE
