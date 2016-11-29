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
      USE bc, only: bc_defined, bc_type, dimension_bc
      USE bc, only: bc_i_e, bc_i_w, bc_j_s, bc_j_n, bc_k_b, bc_k_t
      USE compar, only: istart2,iend2,jstart2,jend2,kstart2,kend2
      USE compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
      USE functions, only: dead_cell_at, default_wall_at
      USE functions, only: is_on_myPE_plus2layers
      USE geometry, only: imax2, jmax2, kmax2

      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Local index for boundary condition
      INTEGER :: L
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
            IF (DEFAULT_WALL_AT(i1,j1,k1)) CALL SET_WALL_BC1 &
                 (I1, I1, J1, J1, K1, K1)
         ENDDO
      ENDDO

! top xy-plane
      K1 = KMAX2
      DO J1 = JSTART3, JEND3
         DO I1 = ISTART3, IEND3
            IF(K1.NE.KEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IF (DEFAULT_WALL_AT(i1,j1,k1)) CALL SET_WALL_BC1 &
                 (I1, I1, J1, J1, K1, K1)
         ENDDO
      ENDDO

! south xz-plane
      J1 = 1
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
            IF(J1.NE.JSTART2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IF (DEFAULT_WALL_AT(i1,j1,k1)) CALL SET_WALL_BC1 &
                 (I1, I1, J1, J1, K1, K1)
         ENDDO
      ENDDO

! north xz-plane
      J1 = JMAX2
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
            IF(J1.NE.JEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IF (DEFAULT_WALL_AT(i1,j1,k1)) CALL SET_WALL_BC1 &
                 (I1, I1, J1, J1, K1, K1)
         ENDDO
      ENDDO

! west zy-plane
      I1 = 1
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
            IF(I1.NE.ISTART2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IF (DEFAULT_WALL_AT(i1,j1,k1)) CALL SET_WALL_BC1 &
                 (I1, I1, J1, J1, K1, K1)
         ENDDO
      ENDDO

! east zy-plane
      I1 = IMAX2
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
            IF(I1.NE.IEND2) EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IF (DEFAULT_WALL_AT(i1,j1,k1)) CALL SET_WALL_BC1 &
                 (I1, I1, J1, J1, K1, K1)
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
      USE bc
      USE compar, only: istart2,iend2,jstart2,jend2,kstart2,kend2
      USE fldvar, only: u_g, v_g, w_g
      USE functions, only: iplus, iminus, jplus, jminus, kplus, kminus
      USE functions, only: dead_cell_at, ns_wall_at, wall_at
      USE geometry , only: do_k
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
      INTEGER :: I1, I2, J1, J2, K1, K2
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

               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

               IF(NS_WALL_AT(i,j,k))THEN
                  SIGN0 = -ONE
               ELSE
                  SIGN0 = ONE
               ENDIF

               IF (wall_at(i,j,k)) THEN

! Fluid cell at West
                  IF (.NOT.wall_at(iminus(i,j,k),j,k)) THEN
! Wall cell at North
                     IF (wall_at(i,jplus(i,j,k),k)) THEN
                        V_G(I,J,K) = SIGN0*V_G(iminus(i,j,k),j,k)
                     ENDIF
! Wall cell at Top
                     IF (wall_at(i,j,kplus(i,j,k))) THEN
                        W_G(I,J,K) = SIGN0*W_G(iminus(i,j,k),j,k)
                     ENDIF
                  ENDIF

! Fluid cell at East
                  IF (.NOT.wall_at(iplus(i,j,k),j,k)) THEN
! Wall cell at North
                     IF (wall_at(i,jplus(i,j,k),k)) THEN
                        V_G(I,J,K) = SIGN0*V_G(iplus(i,j,k),j,k)
                     ENDIF
! Wall cell at Top
                     IF (wall_at(i,j,kplus(i,j,k))) THEN
                        W_G(I,J,K) = SIGN0*W_G(iplus(i,j,k),j,k)
                     ENDIF
                  ENDIF


! Fluid cell at South
                  IF (.NOT.wall_at(i,jminus(i,j,k),k)) THEN
! Wall cell at East
                     IF (wall_at(iplus(i,j,k),j,k)) THEN
                        U_G(I,J,K) = SIGN0*U_G(i,jminus(i,j,k),k)
                     ENDIF
! Wall cell at Top
                     IF (wall_at(i,j,kplus(i,j,k))) THEN
                        W_G(I,J,K) = SIGN0*W_G(i,jminus(i,j,k),k)
                     ENDIF
                  ENDIF

! Fluid cell at North
                  IF (.NOT.wall_at(i,jplus(i,j,k),k)) THEN
! Wall cell at East
                     IF (wall_at(iplus(i,j,k),j,k)) THEN
                        U_G(I,J,K) = SIGN0*U_G(i,jplus(i,j,k),k)
                     ENDIF
! Wall cell at Top
                     IF (wall_at(i,j,kplus(i,j,k))) THEN
                        W_G(I,J,K) = SIGN0*W_G(i,jplus(i,j,k),k)
                     ENDIF
                  ENDIF


                  IF (DO_K) THEN
! Fluid cell at Bottom
                     IF (.NOT.wall_at(i,j,kminus(i,j,k))) THEN
! Wall cell at East
                        IF (wall_at(iplus(i,j,k),j,k)) THEN
                           U_G(I,J,K) = SIGN0*U_G(i,j,kminus(i,j,k))
                        ENDIF
! Wall cell at North
                        IF (wall_at(i,jplus(i,j,k),k)) THEN
                           V_G(I,J,K) = SIGN0*V_G(i,j,kminus(i,j,k))
                        ENDIF
                     ENDIF

! Fluid cell at Top
                     IF (.NOT.wall_at(i,j,kplus(i,j,k))) THEN
! Wall cell at East
                        IF (wall_at(iplus(i,j,k),j,k)) THEN
                           U_G(I,J,K) = SIGN0*U_G(i,j,kplus(i,j,k))
                        ENDIF
! Wall cell at North
                        IF (wall_at(i,jplus(i,j,k),k)) THEN
                           V_G(I,J,K) = SIGN0*V_G(i,j,kplus(i,j,k))
                        ENDIF
                     ENDIF
                  ENDIF   ! end if (do_k)

               ENDIF   ! end if (wall_at(i,j,k))
            ENDDO   ! end do loop (i = i1, i2)
         ENDDO   ! end do loop (j = j1, j2)
      ENDDO   ! end do loop (k = k1, k2)

      RETURN
      END SUBROUTINE SET_WALL_BC1
