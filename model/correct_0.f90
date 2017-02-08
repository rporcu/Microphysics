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
      SUBROUTINE CORRECT_0(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
           p_g,pp_g,u_g,v_g,w_g,d_e,d_n,d_t)&
           bind(C, name="correct_0")

      USE ur_facs  , only: ur_fac

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

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

      do K = lo(3),hi(3)
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)
              P_G(I,J,K) = P_G(I,J,K) + UR_FAC(1)*PP_G(I,J,K)
          ENDDO
        ENDDO
      ENDDO

      do K = lo(3),hi(3)
        do J = lo(2),hi(2)
          do I = slo(1),hi(1)
            u_g(i,j,k) = u_g(i,j,k) - d_e(i,j,k)*(pp_g(i+1,j,k)-pp_g(i,j,k))
          enddo
        enddo
     enddo

      do K = lo(3),hi(3)
        do J = slo(2),hi(2)
          do I = lo(1),hi(1)
            v_g(i,j,k) = v_g(i,j,k) - d_n(i,j,k)*(pp_g(i,j+1,k)-pp_g(i,j,k))
          enddo
        enddo
      enddo

      do K = slo(3),hi(3)
        do J = lo(2),hi(2)
          do I = lo(1),hi(1)
            w_g(i,j,k) = w_g(i,j,k) - d_t(i,j,k)*(pp_g(i,j,k+1) - pp_g(i,j,k))
          enddo
        enddo
      enddo

      END SUBROUTINE CORRECT_0

END MODULE CORRECT_0_MODULE
