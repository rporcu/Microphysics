MODULE PARTICLES_IN_CELL_MODULE

      use discretelement, only: PINC
      USE discretelement, only: MAX_PIP
      USE discretelement, only: XE, YN, ZT
      use discretelement, only: entering_ghost, exiting_ghost, nonexistent, normal_ghost
      use mpi_funs_des, only: des_par_exchange

      use check_cell_movement_module, only: check_cell_movement

! Number of particles in the I/J/K direction
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K

      USE error_manager, only: init_err_msg, finl_err_msg
      USE desgrid, only: desgrid_pic
      USE geometry, only: imin2, jmin2, kmin2
      USE geometry, only: imax2, jmax2, kmax2

! The number of particles on the current process.
      use discretelement, only: PIP, MAX_PIP
! The number and list of particles in each fluid cell IJK.
      use discretelement, only: PINC, PIC
! The East/North/Top face location of a given I/J/K index.
      use discretelement, only: XE, YN, ZT
! The Upper and Loper indices covered by the current process.
      use compar, only: ISTART3, IEND3
      use compar, only: JSTART3, JEND3
      use compar, only: KSTART3, KEND3
! Fluid grid cell dimensions and mesh size
      USE geometry, only: IMIN2, IMAX2
      USE geometry, only: JMIN2, JMAX2
      USE geometry, only: KMIN2, KMAX2
! Fixed array sizes in the I/J/K direction
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: PARTICLES_IN_CELL                                       !
!                                                                      !
!  Purpose:                                                            !
!     - For each particle find the computational fluid cell            !
!       containing the particle center.                                !
!     - Calculate the bulk density in each computational fluid         !
!       cell.                                                          !
!     - Calculate the volume average solids velocity in each           !
!       computational fluid cell.                                      !
!     - For parallel processing indices are altered                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PARTICLES_IN_CELL(pijk, iglobal_id, particle_state, des_pos_new, des_vel_new)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new
      INTEGER(KIND=1), DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: iglobal_id
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

! Local Variables
!---------------------------------------------------------------------//
! particle no.
      INTEGER L
! accounted for particles
      INTEGER PC
! ijk indices
      INTEGER I, J, K
! variables that count/store the number of particles in i, j, k cell
      INTEGER:: npic, pos
! The accumulated number of particles in each IJK.
      INTEGER, ALLOCATABLE :: PARTICLE_COUNT(:,:,:)
!......................................................................!


! following quantities are reset every call to particles_in_cell
      PINC(:,:,:) = 0


! Use an incremental approach to determine the new particle location.
!-----------------------------------------------------------------------
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE

         I = PIJK(L,1)
         IF(I <= ISTART3 .OR. I >= IEND3) THEN
            CALL PIC_SEARCH(I, DES_POS_NEW(L,1), XE,                   &
               DIMENSION_I, IMIN2, IMAX2)
         ELSE
            IF((DES_POS_NEW(L,1) >= XE(I-1)) .AND.                     &
               (DES_POS_NEW(L,1) <  XE(I))) THEN
               I = I
            ELSEIF((DES_POS_NEW(L,1) >= XE(I)) .AND.                   &
               (DES_POS_NEW(L,1) < XE(I+1))) THEN
              I = I+1
           ELSEIF((DES_POS_NEW(L,1) >= XE(I-2)) .AND.                 &
              (DES_POS_NEW(L,1) < XE(I-1))) THEN
               I = I-1
            ELSE
               CALL PIC_SEARCH(I, DES_POS_NEW(L,1), XE,                &
                  DIMENSION_I, IMIN2, IMAX2)
            ENDIF
         ENDIF

         J = PIJK(L,2)
         IF(J <= JSTART3 .OR. J >= JEND3) THEN
            CALL PIC_SEARCH(J, DES_POS_NEW(L,2), YN,                   &
               DIMENSION_J, JMIN2, JMAX2)
         ELSE
            IF((DES_POS_NEW(L,2) >= YN(J-1)) .AND.                     &
               (DES_POS_NEW(L,2) < YN(J))) THEN
               J = J
            ELSEIF((DES_POS_NEW(L,2) >= YN(J)) .AND.                   &
               (DES_POS_NEW(L,2) < YN(J+1))) THEN
               J = J+1
            ELSEIF((DES_POS_NEW(L,2) >= YN(J-2)) .AND.                 &
               (DES_POS_NEW(L,2) < YN(J-1)))THEN
               J = J-1
            ELSE
               CALL PIC_SEARCH(J, DES_POS_NEW(L,2), YN,                &
                  DIMENSION_J, JMIN2, JMAX2)
            ENDIF
         ENDIF


         K = PIJK(L,3)
         IF(K <= KSTART3 .OR. K >= KEND3) THEN
            CALL PIC_SEARCH(K, DES_POS_NEW(L,3), ZT,                &
               DIMENSION_K, KMIN2, KMAX2)
         ELSE
            IF((DES_POS_NEW(L,3) >= ZT(K-1)) .AND.                  &
               (DES_POS_NEW(L,3) < ZT(K))) THEN
               K = K
            ELSEIF((DES_POS_NEW(L,3) >= ZT(K)) .AND.               &
               (DES_POS_NEW(L,3) < ZT(K+1))) THEN
               K = K+1
            ELSEIF((DES_POS_NEW(L,3) >= ZT(K-2)) .AND.              &
               (DES_POS_NEW(L,3) < ZT(K-1))) THEN
               K = K-1
            ELSE
               CALL PIC_SEARCH(K, DES_POS_NEW(L,3), ZT,             &
                  DIMENSION_K, KMIN2, KMAX2)
            ENDIF
         ENDIF

! Assign PIJK(L,1:3)
         PIJK(L,1) = I
         PIJK(L,2) = J
         PIJK(L,3) = K

! Increment the number of particles in cell IJK
         IF(.NOT.NORMAL_GHOST==PARTICLE_STATE(L) .AND. .NOT.ENTERING_GHOST==PARTICLE_STATE(L) .AND. &
            .NOT.EXITING_GHOST==PARTICLE_STATE(L)) PINC(i,j,k) = PINC(i,j,k) + 1

      ENDDO

      CALL CHECK_CELL_MOVEMENT(pijk, iglobal_id, des_vel_new, des_pos_new)

! Assigning the variable PIC(IJK)%p(:). For each computational fluid
! cell compare the number of current particles in the cell to what was
! in the cell previously. If different reallocate. Store the particle
! ids
! ---------------------------------------------------------------->>>
        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

! checking all cells (including ghost cells); updating entering/exiting
! particle regions
         NPIC =  PINC(i,j,k)
         IF (ASSOCIATED(PIC(I,J,K)%p)) THEN
            IF (NPIC.NE.SIZE(PIC(I,J,K)%p)) THEN
               DEALLOCATE(PIC(I,J,K)%p)
               IF (NPIC.GT.0) ALLOCATE(PIC(I,J,K)%p(NPIC))
            ENDIF
         ELSE
            IF (NPIC.GT.0) ALLOCATE(PIC(I,J,K)%p(NPIC))
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      allocate( PARTICLE_COUNT(istart3:iend3, jstart3:jend3, kstart3:kend3))
      PARTICLE_COUNT(:,:,:) = 1
      PC = 1
      DO L = 1, MAX_PIP
! exiting loop if reached max number of particles in processor
         IF(PC.GT.PIP) exit
! skipping indices with no particles (non-existent particles)
         IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE
! incrementing particle account when particle exists
         PC = PC+1
! skipping ghost particles
         IF(NORMAL_GHOST==PARTICLE_STATE(L) .OR. ENTERING_GHOST==PARTICLE_STATE(L) .OR. EXITING_GHOST==PARTICLE_STATE(L)) CYCLE
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         POS = PARTICLE_COUNT(I,J,K)
         PIC(I,J,K)%P(POS) = L
         PARTICLE_COUNT(I,J,K) = PARTICLE_COUNT(I,J,K) + 1
      ENDDO

      deallocate(PARTICLE_COUNT)

      RETURN
      END SUBROUTINE PARTICLES_IN_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: INIT_PARTICLES_IN_CELL                                  !
!                                                                      !
!  Purpose:                                                            !
!     - For each particle find the computational fluid cell            !
!       containing the particle center.                                !
!     - Calculate the bulk density in each computational fluid         !
!       cell.                                                          !
!     - Calculate the volume average solids velocity in each           !
!       computational fluid cell.                                      !
!     - For parallel processing indices are altered                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_PARTICLES_IN_CELL(pijk, dg_pijk, dg_pijkprv, particle_state, des_pos_new)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: des_pos_new
      INTEGER(KIND=1), DIMENSION(:), INTENT(IN) :: particle_state
      INTEGER, DIMENSION(:), INTENT(INOUT) :: dg_pijk
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijkprv
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER :: L
! ijk indices
      INTEGER :: I, J, K

      CALL INIT_ERR_MSG("INIT_PARTICLES_IN_CELL")

! following quantities are reset every call to particles_in_cell
      PINC(:,:,:) = 0

! Bin the particles to the DES grid.
      CALL DESGRID_PIC(.TRUE., dg_pijkprv, dg_pijk, particle_state, des_pos_new)
! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
      CALL DES_PAR_EXCHANGE(des_pos_new, dg_pijk, dg_pijkprv, particle_state)

! Assigning PIJK(L,1), PIJK(L,2) and PIJK(L,3) the i, j, k indices
! of particle L (locating particle on fluid grid). Also determine
! composite ijk index. If first_pass, also assigning particle_phase(L,5) the
! solids phase index of particle.
! ---------------------------------------------------------------->>>
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE

! Use a brute force technique to determine the particle locations in
! the Eulerian fluid grid.

         CALL PIC_SEARCH(I, DES_POS_NEW(L,1), XE,                      &
            DIMENSION_I, IMIN2, IMAX2)
         PIJK(L,1) = I

         CALL PIC_SEARCH(J, DES_POS_NEW(L,2), YN,                      &
            DIMENSION_J, JMIN2, JMAX2)
         PIJK(L,2) = J

         CALL PIC_SEARCH(K, DES_POS_NEW(L,3), ZT,                   &
            DIMENSION_K, KMIN2, KMAX2)
         PIJK(L,3) = K

! Enumerate the number of 'real' particles in the ghost cell.
         IF(.NOT.NORMAL_GHOST==PARTICLE_STATE(L) .AND. &
            .NOT.ENTERING_GHOST==PARTICLE_STATE(L) .AND. &
            .NOT.EXITING_GHOST==PARTICLE_STATE(L)) PINC(i,j,k) = PINC(i,j,k) + 1
      ENDDO

! Bin the particles to the DES grid.
      CALL DESGRID_PIC(.TRUE., dg_pijkprv, dg_pijk, particle_state, des_pos_new)
! Calling exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
! unclear why this needs to be called again.
      CALL DES_PAR_EXCHANGE(des_pos_new, dg_pijk, dg_pijkprv, particle_state)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE INIT_PARTICLES_IN_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PIC_SEARCH                                              !
!                                                                      !
!  Purpose: Identify the I (or J or K) index of the fluid cell that    !
!  contains the particle centroid.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_SEARCH(IDX, lPOS, ENT_POS, lDIMN, lSTART, lEND)

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Index being searched for (I, J, or K)
      INTEGER, INTENT(OUT) :: IDX
! Particle x,y,z position
      DOUBLE PRECISION, INTENT(IN) :: lPOS
! Dimension of ENT_POS array
      INTEGER, INTENT(IN) :: lDIMN
! East, North, or Top cell face location
      DOUBLE PRECISION, INTENT(IN) :: ENT_POS(0:lDIMN)
! Search bounds (by rank)
      INTEGER, INTENT(IN) :: lSTART, lEND

      DO IDX = lSTART,lEND
         IF(lPOS >= ENT_POS(IDX-1) .AND. lPOS < ENT_POS(IDX)) EXIT
      ENDDO

      RETURN
      END SUBROUTINE PIC_SEARCH
END MODULE PARTICLES_IN_CELL_MODULE
