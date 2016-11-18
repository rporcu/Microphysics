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
      USE indices
      USE physprop
      USE compar
      USE cutcell
      USE quadric
      USE functions

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

         IF (FLUIDORP_FLOW_AT(IJK)) THEN

            P_G(IJK) = P_G(IJK) + UR_FAC(1)*PP_G(IJK)

            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IF(.NOT.CARTESIAN_GRID) THEN
               U_G(IJK) = U_G(IJK) - D_E(IJK)*(PP_G(IJKE)-PP_G(IJK))
               V_G(IJK) = V_G(IJK) - D_N(IJK)*(PP_G(IJKN)-PP_G(IJK))
               IF (DO_K) THEN
                  IJKT = TOP_OF(IJK)
                  W_G(IJK) = W_G(IJK) - D_T(IJK)*(PP_G(IJKT)-PP_G(IJK))
               ENDIF
            ELSE
               U_G(IJK) = U_G(IJK) - D_E(IJK)*&
                  (PP_G(IJKE)*A_UPG_E(IJK) - PP_G(IJK)*A_UPG_W(IJK))
               V_G(IJK) = V_G(IJK) - D_N(IJK)*&
                  (PP_G(IJKN)*A_VPG_N(IJK) - PP_G(IJK)*A_VPG_S(IJK))
               IF (DO_K) THEN
                  IJKT = TOP_OF(IJK)
                  W_G(IJK) = W_G(IJK) - D_T(IJK)*&
                     (PP_G(IJKT)*A_WPG_T(IJK) - PP_G(IJK)*A_WPG_B(IJK))
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CORRECT_0
