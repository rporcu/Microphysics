MODULE SET_WALL_BC_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_WALL_BC                                             C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!  Purpose: Set wall boundary conditions                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_WALL_BC(u_g,v_g,w_g, flag)&
         bind(C, name="set_wall_bc")

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc, only: bc_defined, bc_type, dimension_bc
      USE bc, only: bc_i_e, bc_i_w, bc_j_s, bc_j_n, bc_k_b, bc_k_t
      USE compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3

      implicit none

      real(c_real), intent(inout) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer(c_int), intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
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
                  CALL SET_WALL_BC1 (I1, I2, J1, J2, K1, K2, u_g, v_g, w_g, flag)

               CASE ('NO_SLIP_WALL')
                  CALL SET_WALL_BC1 (I1, I2, J1, J2, K1, K2, u_g, v_g, w_g, flag)

               CASE ('PAR_SLIP_WALL')
! updating the boundary velocity may improve convergence
            END SELECT
         ENDIF
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
      SUBROUTINE SET_WALL_BC1(II1, II2, JJ1, JJ2, KK1, KK2, &
         u_g, v_g, w_g, flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: one
      USE compar   , only: istart2,iend2,jstart2,jend2,kstart2,kend2
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE functions, only: iplus, iminus, jplus, jminus, kplus, kminus
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

      real(c_real), INTENT(INOUT) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(INOUT) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(INOUT) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Sign with legal values +1 or -1
      real(c_real) :: SIGN0
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

               IF(flag(i,j,k,1) == 100)THEN
                  SIGN0 = -ONE
               ELSE
                  SIGN0 = ONE
               ENDIF

               IF (flag(i,j,k,1) >= 100) THEN

! Fluid cell at West
                  IF (flag(iminus(i,j,k),j,k,1) < 100) THEN
! Wall cell at North
                     IF (flag(i,jplus(i,j,k),k,1) >= 100) THEN
                        V_G(I,J,K) = SIGN0*V_G(iminus(i,j,k),j,k)
                     ENDIF
! Wall cell at Top
                     IF (flag(i,j,kplus(i,j,k),1)>=100) THEN
                        W_G(I,J,K) = SIGN0*W_G(iminus(i,j,k),j,k)
                     ENDIF
                  ENDIF

! Fluid cell at East
                  IF (flag(iplus(i,j,k),j,k,1)<100) THEN
! Wall cell at North
                     IF (flag(i,jplus(i,j,k),k,1)>=100) THEN
                        V_G(I,J,K) = SIGN0*V_G(iplus(i,j,k),j,k)
                     ENDIF
! Wall cell at Top
                     IF (flag(i,j,kplus(i,j,k),1)>=100) THEN
                        W_G(I,J,K) = SIGN0*W_G(iplus(i,j,k),j,k)
                     ENDIF
                  ENDIF


! Fluid cell at South
                  IF (flag(i,jminus(i,j,k),k,1)<100) THEN
! Wall cell at East
                     IF (flag(iplus(i,j,k),j,k,1)>=100) THEN
                        U_G(I,J,K) = SIGN0*U_G(i,jminus(i,j,k),k)
                     ENDIF
! Wall cell at Top
                     IF (flag(i,j,kplus(i,j,k),1)>=100) THEN
                        W_G(I,J,K) = SIGN0*W_G(i,jminus(i,j,k),k)
                     ENDIF
                  ENDIF

! Fluid cell at North
                  IF (flag(i,jplus(i,j,k),k,1)<100) THEN
! Wall cell at East
                     IF (flag(iplus(i,j,k),j,k,1)>=100) THEN
                        U_G(I,J,K) = SIGN0*U_G(i,jplus(i,j,k),k)
                     ENDIF
! Wall cell at Top
                     IF (flag(i,j,kplus(i,j,k),1)>=100) THEN
                        W_G(I,J,K) = SIGN0*W_G(i,jplus(i,j,k),k)
                     ENDIF
                  ENDIF


! Fluid cell at Bottom
                  IF (flag(i,j,kminus(i,j,k),1)<100) THEN
! Wall cell at East
                     IF (flag(iplus(i,j,k),j,k,1)>=100) THEN
                        U_G(I,J,K) = SIGN0*U_G(i,j,kminus(i,j,k))
                     ENDIF
! Wall cell at North
                     IF (flag(i,jplus(i,j,k),k,1)>=100) THEN
                        V_G(I,J,K) = SIGN0*V_G(i,j,kminus(i,j,k))
                     ENDIF
                  ENDIF

! Fluid cell at Top
                  IF (flag(i,j,kplus(i,j,k),1)<100) THEN
! Wall cell at East
                     IF (flag(iplus(i,j,k),j,k,1)>=100) THEN
                        U_G(I,J,K) = SIGN0*U_G(i,j,kplus(i,j,k))
                     ENDIF
! Wall cell at North
                     IF (flag(i,jplus(i,j,k),k,1)>=100) THEN
                        V_G(I,J,K) = SIGN0*V_G(i,j,kplus(i,j,k))
                     ENDIF
                  ENDIF

               ENDIF   ! end if (wall_at(i,j,k))
            ENDDO   ! end do loop (i = i1, i2)
         ENDDO   ! end do loop (j = j1, j2)
      ENDDO   ! end do loop (k = k1, k2)

      RETURN
      END SUBROUTINE SET_WALL_BC1
END MODULE SET_WALL_BC_MODULE
