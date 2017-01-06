MODULE CALC_TRD_G_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_trD_g                                             C
!  Purpose: Calculate the trace of gas phase rate of strain tensor     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TRD_G(slo,shi,lo,hi,trd_g,u_g,v_g,w_g,flag,dx,dy,dz)

      USE functions, only: iminus, jminus, kminus
      USE compar, only: istart, iend, jstart, jend, kstart, kend

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
      real(c_real), intent(inout) :: trd_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dx, dy, dz

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! Indices
      INTEGER :: I, J, K

      real(c_real) :: odx, ody, odz

      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)
            IF (flag(i,j,k,1)<100) THEN

              TRD_G(i,j,k) = &
                 (U_G(I,J,K)-U_G(iminus(i,j,k),j,k))*ODX + &
                 (V_G(I,J,K)-V_G(i,jminus(i,j,k),k))*ODY + &
                 (W_G(I,J,K)-W_G(i,j,kminus(i,j,k)))*ODZ

             ENDIF
          END DO
        END DO
      END DO

      END SUBROUTINE CALC_TRD_G
END MODULE CALC_TRD_G_MODULE
