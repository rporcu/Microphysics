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
      USE compar
      USE discretelement
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      integer :: i,j,k,ijk

      do k = kstart3, kend3
        do j = jstart3, jend3
          do i = istart3, iend3


            IF (.NOT.WALL_AT(i,j,k)) THEN
               IF (ip_at_e(i,j,k)) U_G(I,J,K) = ZERO
               IF (ip_at_n(i,j,k)) V_G(I,J,K) = ZERO
               IF (ip_at_t(i,j,k)) W_G(I,J,K) = ZERO
            ELSE

               U_G(I,J,K) = ZERO
               V_G(I,J,K) = ZERO
               W_G(I,J,K) = ZERO

               IF (.NOT.(CYCLIC_AT(I,J,K) .AND. (I==IMAX2 .OR. &
                   I==IMAX3))) U_G(IMinus(i,j,k),J,K) = ZERO
               IF (.NOT.(CYCLIC_AT(I,J,K) .AND. (J==JMAX2 .OR. &
                   J==JMAX3))) V_G(I,JMinus(i,j,k),K) = ZERO
               IF (.NOT.(CYCLIC_AT(I,J,K) .AND. (K==KMAX2 .OR. &
                   K==KMAX3))) W_G(I,J,KMinus(i,j,k)) = ZERO
            ENDIF

          end do
        end do
      end do

      END SUBROUTINE ZERO_NORM_VEL
