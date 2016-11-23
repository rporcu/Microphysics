!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_BC_AREA                                            C
!  Purpose: Compute area of boundary surfaces                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_BC_AREA
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE bc
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
! BC number
      INTEGER :: BCV
!-----------------------------------------------
!

      DO BCV = 1, DIMENSION_BC
         IF (BC_DEFINED(BCV)) THEN
            BC_AREA(BCV) = ZERO
            IF (BC_PLANE(BCV) == 'W' .OR. BC_PLANE(BCV) == 'E') THEN
               BC_AREA(BCV) = &
                  DY*dble(BC_J_N(BCV)-BC_J_S(BCV)+1)* &
                  DZ*dble(BC_K_T(BCV)-BC_K_B(BCV)+1)
            ELSE IF (BC_PLANE(BCV)=='S' .OR. BC_PLANE(BCV)=='N') THEN
               BC_AREA(BCV) = &
                  DX*dble(BC_I_E(BCV)-BC_I_W(BCV)+1)* &
                  DZ*dble(BC_K_T(BCV)-BC_K_B(BCV)+1)
            ELSE IF (BC_PLANE(BCV)=='B' .OR. BC_PLANE(BCV)=='T') THEN
               BC_AREA(BCV) = &
                  DX*dble(BC_I_E(BCV)-BC_I_W(BCV)+1)* &
                  DY*dble(BC_J_N(BCV)-BC_J_S(BCV)+1)
            ENDIF
         ENDIF
      ENDDO


      RETURN
      END SUBROUTINE GET_BC_AREA
