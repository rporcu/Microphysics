!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RESET_NEW                                              C
!  Purpose: Reset the new variables with the stored previous-time-step C
!           values of field variables.                                 C
!    *****Remember to modify update_old also
!                                                                      C
!  Author: M. Syamlal                                 Date: FEB-6-97   C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: ROP_g, EP_g, ROP_s, IJKMAX2, MMAX, U_s, V_s,  C
!                        W_s                                           C
!                                                                      C
!  Variables modified: ROP_go, ROP_so, IJK, M, U_so, V_so, W_so C
!                                                                      C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE RESET_NEW

!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run

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
!                    Indices
      INTEGER :: M
!
!                    error index
      INTEGER :: IER

!-----------------------------------------------

      EP_G(:) = EP_GO(:)
      P_G(:) = P_GO(:)
      RO_G(:) = RO_GO(:)
      ROP_G(:) = ROP_GO(:)
      U_G(:) = U_GO(:)
      V_G(:) = V_GO(:)
      W_G(:) = W_GO(:)

!     Recalculate all coefficients
      CALL CALC_COEFF_ALL (0, IER)

      RETURN
      END SUBROUTINE RESET_NEW

!// Comments on the modifications for DMP version implementation
!// 120 Replaced the index for initialization: (:IJKMAX2) to just (:)
