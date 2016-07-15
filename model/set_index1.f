!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_INDEX1                                             !
!  Author: M. Syamlal, W. A. Rogers                   Date: 17-Dec-91  C

!  Purpose: Set the indices of the first neighbors of cell ijk         C
!           This version adds 'increments' stored in STORE_INCREMENTS  C
!           to IJK to find indices of the neighbors.                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_INDEX1(IJK, I, J, K, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
         IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, IM, JM, KM)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE fldvar
      USE geometry
      USE constant
      USE indices
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
                      IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, &
                      IM, JM, KM
!
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
!
      IM = IM1(I)
      JM = JM1(J)
      KM = KM1(K)

!     Determine the true indices of neighboring cells
      IJKW = WEST_OF(IJK)
      IJKE = EAST_OF(IJK)
      IJKS = SOUTH_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKB = BOTTOM_OF(IJK)
      IJKT = TOP_OF(IJK)
      IMJK = IM_OF(IJK)
      IPJK = IP_OF(IJK)
      IJMK = JM_OF(IJK)
      IJPK = JP_OF(IJK)
      IJKM = KM_OF(IJK)
      IJKP = KP_OF(IJK)
!
      RETURN
      END SUBROUTINE SET_INDEX1
