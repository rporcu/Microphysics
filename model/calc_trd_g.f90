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
      SUBROUTINE CALC_TRD_G(trd_g,u_g,v_g,w_g,flag,dx,dy,dz)

      USE functions, only: iminus, jminus, kminus
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE compar, only: istart, iend, jstart, jend, kstart, kend

      IMPLICIT NONE

!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
     double precision, intent(inout) :: trd_g(istart:iend,jstart:jend,kstart:kend)

      real(c_real), intent(in   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(in   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(in   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      integer(c_int), intent(in   ) :: flag&
            (istart3:iend3,jstart3:jend3,kstart3:kend3,4)
      real(c_real), intent(in   ) :: dx, dy, dz

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! Indices
      INTEGER :: I, J, K

      double precision :: odx, ody, odz

      odx = 1.d0 / dx
      ody = 1.d0 / dy
      odz = 1.d0 / dz
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
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
