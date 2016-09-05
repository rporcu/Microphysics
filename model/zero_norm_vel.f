!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ZERO_NORM_VEL                                           C
!  Purpose: Set the velocity component normal to a wall to zero        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ZERO_NORM_VEL

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE physprop
      USE fldvar
      USE indices
      USE compar
      USE discretelement
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!
      INTEGER :: ISV
! Indicies
      INTEGER :: IJK, IMJK, IJMK, IJKM
      INTEGER :: M

!!$omp  parallel do private( IMJK, IJMK, IJKM)

      DO IJK = ijkstart3, ijkend3

         IF (.NOT.WALL_AT(IJK)) THEN
            IF (IP_AT_E(IJK)) U_G(IJK) = ZERO
            IF (IP_AT_N(IJK)) V_G(IJK) = ZERO
            IF (IP_AT_T(IJK)) W_G(IJK) = ZERO
         ELSE
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            U_G(IJK) = ZERO
            V_G(IJK) = ZERO
            W_G(IJK) = ZERO
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (I_OF(IJK)==IMAX2 .OR. &
                I_OF(IJK)==IMAX3))) U_G(IMJK) = ZERO
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (J_OF(IJK)==JMAX2 .OR. &
                J_OF(IJK)==JMAX3))) V_G(IJMK) = ZERO
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (K_OF(IJK)==KMAX2 .OR. &
                K_OF(IJK)==KMAX3))) W_G(IJKM) = ZERO
         ENDIF
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)


      RETURN

      END SUBROUTINE ZERO_NORM_VEL



