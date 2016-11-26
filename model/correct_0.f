!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CORRECT_0                                               C
!  Purpose: Correct the fluid pressure and gas velocities              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CORRECT_0()

      USE param
      USE param1
      USE geometry
      USE physprop
      USE compar
      USE functions, only: funijk, ieast, jnorth, ktop
      USE functions, only: fluidorp_flow_at

      use fldvar
      USE ur_facs

      IMPLICIT NONE


!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I,J,K,IJK
!-----------------------------------------------

! Underrelax pressure correction.  Velocity corrections should not be
! underrelaxed, so that the continuity eq. is satisfied.

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3

            IJK = FUNIJK(i,j,k)
            IF (FLUIDORP_FLOW_AT(i,j,k)) THEN

              P_G(IJK) = P_G(IJK) + UR_FAC(1)*PP_G(I,J,K)

              U_G(I,J,K) = U_G(I,J,K) - &
                 D_E(I,J,K)*(PP_G(ieast(i,j,k),j,k)-PP_G(I,J,K))
              V_G(I,J,K) = V_G(I,J,K) - &
                 D_N(I,J,K)*(PP_G(i,jnorth(i,j,k),k)-PP_G(I,J,K))
              IF (DO_K) W_G(I,J,K) = W_G(I,J,K) - &
                    D_T(I,J,K)*(PP_G(i,j,ktop(i,j,k)) - PP_G(I,J,K))

            ENDIF

          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CORRECT_0
