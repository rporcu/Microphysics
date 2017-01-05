MODULE CORRECT_0_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CORRECT_0                                               C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!                                                                      C
!  Purpose: Correct the fluid pressure and gas velocities              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CORRECT_0(p_g,pp_g,u_g,v_g,w_g,d_e,d_n,d_t,flag)&
         bind(C, name="correct0")

      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE functions, only: ieast, jnorth, ktop
      USE ur_facs  , only: ur_fac

      IMPLICIT NONE

      real(c_real), intent(inout) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: pp_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: d_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: d_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: d_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer(c_int), intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)


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

            if(flag(i,j,k,1) <= 11) then

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
END MODULE CORRECT_0_MODULE
