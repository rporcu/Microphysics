MODULE GET_BC_AREA_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS
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
      SUBROUTINE GET_BC_AREA(dx,dy,dz)

      use param, only: dimension_bc
      use param1, only: zero
      use bc, only: bc_defined, bc_area, bc_plane, bc_k_b, bc_k_t, bc_defined, bc_j_n, bc_j_s, bc_i_e, bc_i_w

      IMPLICIT NONE

      real(c_real), intent(in) :: dx, dy, dz
!
! BC number
      integer :: BCV
!-----------------------------------------------
!

      DO BCV = 1, DIMENSION_BC
         IF (BC_DEFINED(BCV)) THEN
            BC_AREA(BCV) = ZERO
            IF (bc_plane(BCV) == 'W' .OR. bc_plane(BCV) == 'E') THEN
               BC_AREA(BCV) = &
                  DY*dble(BC_J_N(BCV)-BC_J_S(BCV)+1)* &
                  DZ*dble(BC_K_T(BCV)-BC_K_B(BCV)+1)
            ELSE IF (bc_plane(BCV)=='S' .OR. bc_plane(BCV)=='N') THEN
               BC_AREA(BCV) = &
                  DX*dble(BC_I_E(BCV)-BC_I_W(BCV)+1)* &
                  DZ*dble(BC_K_T(BCV)-BC_K_B(BCV)+1)
            ELSE IF (bc_plane(BCV)=='B' .OR. bc_plane(BCV)=='T') THEN
               BC_AREA(BCV) = &
                  DX*dble(BC_I_E(BCV)-BC_I_W(BCV)+1)* &
                  DY*dble(BC_J_N(BCV)-BC_J_S(BCV)+1)
            ENDIF
         ENDIF
      ENDDO


      RETURN
      END SUBROUTINE GET_BC_AREA
END MODULE GET_BC_AREA_MODULE
