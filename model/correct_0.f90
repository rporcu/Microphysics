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
      SUBROUTINE CORRECT_0(slo,shi,lo,hi,p_g,pp_g, &
           u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
           d_e,d_n,d_t,flag)&
           bind(C, name="correct0")

      USE functions, only: ieast, jnorth, ktop
      USE ur_facs  , only: ur_fac

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

      real(c_real), intent(in   ) :: pp_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I,J,K
!-----------------------------------------------

! Underrelax pressure correction.  Velocity corrections should not be
! underrelaxed, so that the continuity eq. is satisfied.

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)

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
