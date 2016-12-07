MODULE SET_PHASE_INDEX_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_PHASE_INDEX                                         !
!                                                                      !
!  Purpose: Set the index of all particles based on their diameter and !
!  density.                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_PHASE_INDEX(PIJK)

      USE discretelement, only: DES_RADIUS, RO_SOL
      USE discretelement, only: MAX_PIP
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar
      USE discretelement, only: nonexistent, normal_ghost, entering_ghost, exiting_ghost, particle_state
      USE mpi_funs_des, only: des_par_exchange
      USE open_files_mod, only: open_pe_log
      USE param1, only: small_number
      USE constant, only: MMAX, D_P0, RO_s0

      IMPLICIT NONE

      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER :: L
! solids phase no.
      INTEGER :: M
! IER for error reporting
      INTEGER :: IER
! Difference between a particles diameter (density) and the diameter
! (density) of a phase specified in the data file.
      DOUBLE PRECISION dDp, dRho

! The restart file contains the phase index for reacting cases as the
! diameter and/or density of the particle may have changed.

! Initialize the error flag.
      IER = 0

! solids phase index of particle.
! ---------------------------------------------------------------->>>
      DO L = 1, MAX_PIP
         IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE
         IF(NORMAL_GHOST==PARTICLE_STATE(L) .OR. ENTERING_GHOST==PARTICLE_STATE(L) .OR. EXITING_GHOST==PARTICLE_STATE(L)) CYCLE

! Determining the solids phase of each particle by matching the diameter
! and density to those specified in the data file.
         M_LP: DO M = 1, MMAX
            dDp  = ABS(2.0d0*DES_RADIUS(L)-D_P0(M))
            dRho = ABS( RO_Sol(L)-RO_S0(M))
            IF( dDp < SMALL_NUMBER .AND. dRho < SMALL_NUMBER) THEN
               PIJK(L,5) = M
               EXIT M_LP
            ENDIF
         ENDDO M_LP
! Flag error if no match is found.
         IF(PIJK(L,5).EQ.0) IER = 1
      ENDDO

! Sync up the error flag across all processes.
      ! CALL GLOBAL_ALL_SUM(IER)
      IF(IER == 0) RETURN

! Point of no return: Report errors and abort
!----------------------------------------------------------------------
      CALL INIT_ERR_MSG("SET_PHASE_INDEX")

      CALL OPEN_PE_LOG(IER)

      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1100 FORMAT('Error 1100: Unable to determine the phase of one or ',&
         'more particles.',/8x,'ID',4X,'Diameter',6x,'Density',/)

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE
         IF(NORMAL_GHOST==PARTICLE_STATE(L) .OR. ENTERING_GHOST==PARTICLE_STATE(L) .OR. EXITING_GHOST==PARTICLE_STATE(L)) CYCLE

! Flag as an error if no match is found.
         IF(PIJK(L,5).EQ.0) THEN
            WRITE(ERR_MSG,9000) L,  2.0*DES_RADIUS(L), Ro_Sol(L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

      WRITE(ERR_MSG, 1101)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1101 FORMAT(' ',/'Defined phase parameters from mfix.dat:',/3x,'ID',&
         5X,'Diameter',5x,'Density')

      DO M = 1, MMAX
         WRITE(ERR_MSG, 9000) M, D_P0(M), RO_S0(M)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO

      WRITE(ERR_MSG, 1102)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

 1102 FORMAT('Please correct the mfix.dat or particle_input.dat files.')

 9000 FORMAT(I10,2(2x,g12.5))

      END SUBROUTINE SET_PHASE_INDEX
END MODULE SET_PHASE_INDEX_MODULE
