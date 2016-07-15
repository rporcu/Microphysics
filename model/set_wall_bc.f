!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_WALL_BC                                             C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!  Purpose: Set wall boundary conditions                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_WALL_BC()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE funits
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Local index for boundary condition
      INTEGER :: L
! indices
      INTEGER :: IJK, IPJK
! Starting & ending I index
      INTEGER :: I1, I2
! Starting & ending J index
      INTEGER :: J1, J2
! Starting and ending K index
      INTEGER :: K1, K2
!-----------------------------------------------


! Set the boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! The range of boundary cells
            I1 = BC_I_W(L)
            I2 = BC_I_E(L)
            J1 = BC_J_S(L)
            J2 = BC_J_N(L)
            K1 = BC_K_B(L)
            K2 = BC_K_T(L)

            SELECT CASE (TRIM(BC_TYPE(L)))
               CASE ('FREE_SLIP_WALL')
                  CALL SET_WALL_BC1 (I1, I2, J1, J2, K1, K2)

               CASE ('NO_SLIP_WALL')
                  CALL SET_WALL_BC1 (I1, I2, J1, J2, K1, K2)

               CASE ('PAR_SLIP_WALL')
! updating the boundary velocity may improve convergence
            END SELECT
         ENDIF
      ENDDO


! The above section did not address bc_type=undefined (which by default
! is either a ns wall, a fs wall) or
! bc_type='dummy' conditions. The section below will handle both events
! since default_wall_at will register as true
      K1 = 1
      DO J1 = JSTART3, JEND3
         DO I1 = ISTART3, IEND3
            IF(K1.NE.KSTART2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1,&
                J1, J1, K1, K1)
         ENDDO
      ENDDO

! top xy-plane
      K1 = KMAX2
      DO J1 = JSTART3, JEND3
         DO I1 = ISTART3, IEND3
            IF(K1.NE.KEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
               J1, J1, K1, K1)
         ENDDO
      ENDDO

! south xz-plane
      J1 = 1
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
            IF(J1.NE.JSTART2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
                J1, J1, K1, K1)
         ENDDO
      ENDDO

! north xz-plane
      J1 = JMAX2
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
            IF(J1.NE.JEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
               J1, J1, K1, K1)
         ENDDO
      ENDDO

! west zy-plane
      I1 = 1
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
            IF(I1.NE.ISTART2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
                J1, J1, K1, K1)
         ENDDO
      ENDDO

! east zy-plane
      I1 = IMAX2
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
            IF(I1.NE.IEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, &
                J1, J1, K1, K1)
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE SET_WALL_BC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_WALL_BC1                                            C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!                                                                      C
!  Purpose: Set U, V, and W components for the specified cells by      C
!           copying the same or negative values from near by fluid     C
!           cell                                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_WALL_BC1(II1, II2, JJ1, JJ2, KK1, KK2)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE funits
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Starting and ending I index
      INTEGER, INTENT(IN) :: II1, II2
! Starting and ending J index
      INTEGER, INTENT(IN) :: JJ1, JJ2
! Starting and ending K index
      INTEGER, INTENT(IN) :: KK1, KK2
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Sign with legal values +1 or -1
      DOUBLE PRECISION :: SIGN0
! Local indices near wall cell
      INTEGER :: I, J, K
      INTEGER :: IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP
      INTEGER :: I1, I2, J1, J2, K1, K2
! Local index for a fluid cell near the wall cell
      INTEGER :: LFLUID
!-----------------------------------------------

! Limit I1, I2 and all to local processor first ghost layer
      I1 = II1
      I2 = II2
      J1 = JJ1
      J2 = JJ2
      K1 = KK1
      K2 = KK2

      IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)
      IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)
      IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)
      IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)
      IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)
      IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

      DO K = K1, K2
         DO J = J1, J2
            DO I = I1, I2
               IJK = FUNIJK(I,J,K)

               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells


               IF(NS_WALL_AT(IJK))THEN
                  SIGN0 = -ONE
               ELSE
                  SIGN0 = ONE
               ENDIF

               IF (WALL_AT(IJK)) THEN
                  IMJK = IM_OF(IJK)
                  IJMK = JM_OF(IJK)
                  IJKM = KM_OF(IJK)
                  IPJK = IP_OF(IJK)
                  IJPK = JP_OF(IJK)
                  IJKP = KP_OF(IJK)

! Fluid cell at West
                  IF (.NOT.WALL_AT(IMJK)) THEN
                     LFLUID = IMJK
! Wall cell at North
                     IF (WALL_AT(IJPK)) THEN
                        V_G(IJK) = SIGN0*V_G(LFLUID)
                     ENDIF
! Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN
                        W_G(IJK) = SIGN0*W_G(LFLUID)
                     ENDIF
                  ENDIF

! Fluid cell at East
                  IF (.NOT.WALL_AT(IPJK)) THEN
                     LFLUID = IPJK
! Wall cell at North
                     IF (WALL_AT(IJPK)) THEN
                        V_G(IJK) = SIGN0*V_G(LFLUID)
                     ENDIF
! Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN
                        W_G(IJK) = SIGN0*W_G(LFLUID)
                     ENDIF
                  ENDIF


! Fluid cell at South
                  IF (.NOT.WALL_AT(IJMK)) THEN
                     LFLUID = IJMK
! Wall cell at East
                     IF (WALL_AT(IPJK)) THEN
                        U_G(IJK) = SIGN0*U_G(LFLUID)
                     ENDIF
! Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN
                        W_G(IJK) = SIGN0*W_G(LFLUID)
                     ENDIF
                  ENDIF

! Fluid cell at North
                  IF (.NOT.WALL_AT(IJPK)) THEN
                     LFLUID = IJPK
! Wall cell at East
                     IF (WALL_AT(IPJK)) THEN
                        U_G(IJK) = SIGN0*U_G(LFLUID)
                     ENDIF
! Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN
                        W_G(IJK) = SIGN0*W_G(LFLUID)
                     ENDIF
                  ENDIF


                  IF (DO_K) THEN
! Fluid cell at Bottom
                     IF (.NOT.WALL_AT(IJKM)) THEN
                        LFLUID = IJKM
! Wall cell at East
                        IF (WALL_AT(IPJK)) THEN
                           U_G(IJK) = SIGN0*U_G(LFLUID)
                        ENDIF
! Wall cell at North
                        IF (WALL_AT(IJPK)) THEN
                           V_G(IJK) = SIGN0*V_G(LFLUID)
                        ENDIF
                     ENDIF

! Fluid cell at Top
                     IF (.NOT.WALL_AT(IJKP)) THEN
                        LFLUID = IJKP
! Wall cell at East
                        IF (WALL_AT(IPJK)) THEN
                           U_G(IJK) = SIGN0*U_G(LFLUID)
                        ENDIF
! Wall cell at North
                        IF (WALL_AT(IJPK)) THEN
                           V_G(IJK) = SIGN0*V_G(LFLUID)
                        ENDIF
                     ENDIF
                  ENDIF   ! end if (do_k)

               ENDIF   ! end if (wall_at(ijk))
            ENDDO   ! end do loop (i = i1, i2)
         ENDDO   ! end do loop (j = j1, j2)
      ENDDO   ! end do loop (k = k1, k2)

      RETURN
      END SUBROUTINE SET_WALL_BC1
