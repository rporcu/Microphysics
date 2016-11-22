!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_trD_g(trD_g, IER)                                 C
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
      SUBROUTINE CALC_TRD_G
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE fldvar
      USE compar
      USE sendrecv
      USE bc
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
      INTEGER          I, J, K, IJK, IMJK, IJMK, IJKM, IM
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      DOUBLE PRECISION :: dudx,dvdy,dwdz
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: UW_g,VW_g,WW_g

      LOGICAL :: U_NODE_AT_E, U_NODE_AT_W
      LOGICAL :: V_NODE_AT_N, V_NODE_AT_S
      LOGICAL :: W_NODE_AT_T, W_NODE_AT_B
      INTEGER :: BCV
      INTEGER :: BCT


      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
         IF (.NOT.WALL_AT(IJK)) THEN
            IM = IM1(I)
            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJKM = FUNIJK(i,j,kminus(i,j,k))

            TRD_G(IJK) = (U_G(IJK)-U_G(IMJK))*ODX(I) + (&
              V_G(IJK)-V_G(IJMK))*ODY(J) + (W_G(IJK)-W_G(IJKM))*(ODZ(K))


! Original term:
!            TRD_G(IJK) = (X_E(I)*U_G(IJK)-X_E(IM)*U_G(IMJK))*OX(I)*ODX(I) + (&
!               V_G(IJK)-V_G(IJMK))*ODY(J) + (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K))
         ENDIF
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE CALC_TRD_G



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_CALC_VEL_G_GRAD(IJK,DELV, IER)                      C
!  Purpose: Calculate velocity derivatives in scalar cut-cell          C
!           Gas phase                                                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 25-JAN-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CG_CALC_VEL_G_GRAD(IJK,DELV)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE fldvar
      USE compar
      USE sendrecv
      USE bc
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
      INTEGER          I, J, K, IJK, IMJK, IJMK, IJKM
!
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      DOUBLE PRECISION :: dudx,dudy,dudz
      DOUBLE PRECISION :: dvdx,dvdy,dvdz
      DOUBLE PRECISION :: dwdx,dwdy,dwdz
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: UW_g,VW_g,WW_g
      DOUBLE PRECISION, DIMENSION (3,3) :: DELV

!              |  du/dx    du/dy   du/dz  |
!      DELV =  |  dv/dx    dv/dy   dv/dz  |  =  dUi/dxj
!              |  dw/dx    dw/dy   dw/dz  |

      LOGICAL :: U_NODE_AT_E, U_NODE_AT_W
      LOGICAL :: V_NODE_AT_N, V_NODE_AT_S
      LOGICAL :: W_NODE_AT_T, W_NODE_AT_B
      INTEGER :: BCV
      INTEGER :: BCT
!-----------------------------------------------
!
!
      DELV = ZERO

      RETURN
      END SUBROUTINE CG_CALC_VEL_G_GRAD
