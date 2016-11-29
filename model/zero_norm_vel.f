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

      SUBROUTINE ZERO_NORM_VEL(u_g,v_g,w_g)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1   , only: zero
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE geometry , only: imax2, jmax2, kmax2
      USE geometry , only: imax3, jmax3, kmax3
      USE functions, only: iminus, jminus, kminus
      USE functions, only: ip_at_e, ip_at_n, ip_at_t, cyclic_at, wall_at

      IMPLICIT NONE

      double precision, intent(inout) ::  u_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  v_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  w_g(istart3:iend3,jstart3:jend3,kstart3:kend3)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      integer :: i,j,k

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
