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
      SUBROUTINE CORRECT_0(slo,shi,lo,hi,&
           p_g,pp_g,u_g,v_g,w_g,d_e,d_n,d_t)&
           bind(C, name="correct0")

      USE ur_facs  , only: ur_fac

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      real(c_real), intent(in   ) :: pp_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I,J,K
!-----------------------------------------------

! Underrelax pressure correction.  Velocity corrections should not be
! underrelaxed, so that the continuity eq. is satisfied.

      do K = lo(3),hi(3)
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)
              P_G(I,J,K) = P_G(I,J,K) + UR_FAC(1)*PP_G(I,J,K)
          ENDDO
        ENDDO
      ENDDO

      do K = slo(3),shi(3)
        do J = slo(2),shi(2)
          do I = slo(1),hi(1)
            U_G(I,J,K) = U_G(I,J,K) - D_E(I,J,K)*(PP_G(i+1,j,k)-PP_G(I,J,K))
          ENDDO
        ENDDO
      ENDDO
      do K = slo(3),shi(3)
        do J = slo(2),hi(2)
          do I = slo(1),shi(1)
            V_G(I,J,K) = V_G(I,J,K) - D_N(I,J,K)*(PP_G(i,j+1,k)-PP_G(I,J,K))
          ENDDO
        ENDDO
      ENDDO

      do K = slo(3),hi(3)
        do J = slo(2),shi(2)
          do I = slo(1),shi(1)
            W_G(I,J,K) = W_G(I,J,K) - D_T(I,J,K)*(PP_G(i,j,k+1) - PP_G(I,J,K))
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE CORRECT_0

END MODULE CORRECT_0_MODULE
