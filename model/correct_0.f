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

      SUBROUTINE CORRECT_0(p_g,pp_g,u_g,v_g,w_g,d_e,d_n,d_t)

      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE functions, only: ieast, jnorth, ktop
      USE functions, only: fluidorp_flow_at
      USE ur_facs  , only: ur_fac

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: pp_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: d_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: d_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: d_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)


!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I,J,K
!-----------------------------------------------

! Underrelax pressure correction.  Velocity corrections should not be
! underrelaxed, so that the continuity eq. is satisfied.

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3

            IF (FLUIDORP_FLOW_AT(i,j,k)) THEN

              P_G(I,J,K) = P_G(I,J,K) + UR_FAC(1)*PP_G(I,J,K)

              U_G(I,J,K) = U_G(I,J,K) - &
                 D_E(I,J,K)*(PP_G(ieast(i,j,k),j,k)-PP_G(I,J,K))
              V_G(I,J,K) = V_G(I,J,K) - &
                 D_N(I,J,K)*(PP_G(i,jnorth(i,j,k),k)-PP_G(I,J,K))
              W_G(I,J,K) = W_G(I,J,K) - &
                 D_T(I,J,K)*(PP_G(i,j,ktop(i,j,k)) - PP_G(I,J,K))

            ENDIF

          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CORRECT_0
