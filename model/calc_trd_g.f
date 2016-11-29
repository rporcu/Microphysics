!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_trD_g(IER)g, IER)                                 C
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
      SUBROUTINE CALC_TRD_G(trd_g)
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE fldvar, only: u_g, v_g, w_g
!     USE compar
!     USE bc
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
     double precision, intent(inout) :: trd_g(istart:iend,jstart:jend,kstart:kend)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! Indices
      INTEGER :: I, J, K
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
            IF (.NOT.WALL_AT(i,j,k)) THEN

              TRD_G(i,j,k) = &
                 (U_G(I,J,K)-U_G(iminus(i,j,k),j,k))*ODX + &
                 (V_G(I,J,K)-V_G(i,jminus(i,j,k),k))*ODY + &
                 (W_G(I,J,K)-W_G(i,j,kminus(i,j,k)))*ODZ


         ENDIF
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE CALC_TRD_G
