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
! Indices
      integer :: i,j,k
      INTEGER :: IJK, IMJK, IJMK, IJKM

      do k = kstart3, kend3
        do j = jstart3, jend3
          do i = istart3, iend3

          ijk = funijk(i,j,k)

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
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (I==IMAX2 .OR. &
                I==IMAX3))) U_G(IMJK) = ZERO
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (J==JMAX2 .OR. &
                J==JMAX3))) V_G(IJMK) = ZERO
            IF (.NOT.(CYCLIC_AT(IJK) .AND. (K==KMAX2 .OR. &
                K==KMAX3))) W_G(IJKM) = ZERO
         ENDIF

          end do
        end do
      end do

      END SUBROUTINE ZERO_NORM_VEL
