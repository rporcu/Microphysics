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
      USE functions, only: funijk, ieast, iwest, jnorth, jsouth, ktop, kbot
      USE functions, only: iplus, iminus, jplus, jminus, kplus, kminus
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
      IJKW = FUNIJK(iwest(i,j,k),j,k)
      IJKE = FUNIJK(ieast(i,j,k),j,k)
      IJKS = FUNIJK(i,jsouth(i,j,k),k)
      IJKN = FUNIJK(i,jnorth(i,j,k),k)
      IJKB = FUNIJK(i,j,kbot(i,j,k))
      IJKT = FUNIJK(i,j,ktop(i,j,k))

      IMJK = FUNIJK(iminus(i,j,k),j,k)
      IPJK = FUNIJK(iplus(i,j,k),j,k)
      IJMK = FUNIJK(i,jminus(i,j,k),k)
      IJPK = FUNIJK(i,jplus(i,j,k),k)
      IJKM = FUNIJK(i,j,kminus(i,j,k))
      IJKP = FUNIJK(i,j,kplus(i,j,k))
 
      END SUBROUTINE SET_INDEX1
