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
      INTEGER :: I,J,K,IJK, IJKE, IJKN, IJKT
!-----------------------------------------------

! Underrelax pressure correction.  Velocity corrections should not be
! underrelaxed, so that the continuity eq. is satisfied.

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3

            IJK = FUNIJK(i,j,k)
            IF (FLUIDORP_FLOW_AT(i,j,k)) THEN

              P_G(IJK) = P_G(IJK) + UR_FAC(1)*PP_G(IJK)

              IJKE = FUNIJK(ieast(i,j,k),j,k)
              IJKN = FUNIJK(i,jnorth(i,j,k),k)

              U_G(IJK) = U_G(IJK) - D_E(IJK)*(PP_G(IJKE)-PP_G(IJK))
              V_G(IJK) = V_G(IJK) - D_N(IJK)*(PP_G(IJKN)-PP_G(IJK))

              IF (DO_K) THEN
               IJKT = FUNIJK(i,j,ktop(i,j,k))
              ENDIF

            ENDIF

          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CORRECT_0
