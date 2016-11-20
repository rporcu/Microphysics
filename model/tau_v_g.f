!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_V_g                                            C
!  Purpose: Cross terms in the gradient of stress in V_g momentum      c
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
!                                                                      C
!                                                                      C
!  Comments: This routine calculates the components of the gas phase   C
!  viscous stress tensor of the v-momentum equation that cannot be     C
!  cast in the form: mu.grad(v). These components are stored in the    C
!  passed variable, which is then used as a source of the v-momentum   C
!  equation.                                                           C
!                                                                      C
!  The following v component viscous stress tensor terms are           C
!  calculated here:                                                    C
!  > part of d/dy (tau_yy) xdxdydz =>                                  C
!            d/dy (lambda.trcD) xdxdydz =>                             C
!    delta (lambda.trcD)Ap|N-S                                         C
!  > part of 1/x d/dx(x.tau_xy) xdxdydz =>                             C
!            1/x d/dx (x.mu.du/dy) xdxdydz =>                          C
!    delta (x.mu.du/dy)Ayz |E-W                                        C
!  > part of d/dy (tau_xy) xdxdydz =>                                  C
!           d/dy (mu.dv/dy) xdxdydz =>                                 C
!    delta (mu.dv/dx)Axz |N-S                                          C
!  > part of 1/x d/dz (tau_xz) xdxdydz =>                              C
!            1/x d/dz (mu.dw/dy) xdxdydz =>                            C
!    delta (mu.dw/dx)Axy |T-B                                          C
!                                                                      C
!  To reconstitute the full v-momentum gas phase viscous stress        C
!  tensor would require including the the 'diffusional' components     C
!  (i.e., those of the form mu.grad(v)                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_V_G(lTAU_V_G)

! Modules
!---------------------------------------------------------------------//
      USE param
      USE param1
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE sendrecv
      USE compar
      USE functions
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_V_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_V_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JP, KM
      INTEGER :: IJK, IJKN, IJKE, IJKW, IJKT, IJKB
      INTEGER :: IJKTN, IJKBN, IJKNE, IJKNW
      INTEGER :: IJPK, IJMK, IMJK, IJKM
      INTEGER :: IMJPK, IJPKM
      INTEGER :: itmp, jtmp, ktmp
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
!---------------------------------------------------------------------//
!     NOTE -- triply nested functions seem to break things -- hence the
!             use of the *tmp variables below
!---------------------------------------------------------------------//

        DO K = kstart3, kend3
         DO J = jstart3, jend3
          DO I = istart3, iend3

            IJK  = FUNIJK(i,j,k)
            IJKN = FUNIJK(i,jnorth(i,j,k),k)

            EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)
            IF ( .NOT.IP_AT_N(IJK) .AND. EPGA>DIL_EP_S) THEN

               JP = JP1(J)
               IM = IM1(I)
               KM = KM1(K)

               IJMK = FUNIJK(i,jminus(i,j,k),k)
               IMJK = FUNIJK(iminus(i,j,k),j,k)

               ktmp = kminus(i,j,k)
               IJKM  = FUNIJK(i,j,ktmp)
               IJPKM = FUNIJK(i,jplus(i,j,ktmp),ktmp)

               jtmp = jplus(i,j,k)
               IJPK  = FUNIJK(i,jtmp,k)
               IMJPK = FUNIJK(iminus(i,jtmp,k),jtmp,k)

               IJKE = FUNIJK(ieast(i,j,k),j,k)
               IJKW = FUNIJK(iwest(i,j,k),j,k)
               IJKT = FUNIJK(i,j,ktop(i,j,k))
               IJKB = FUNIJK(i,j,kbot(i,j,k))

               jtmp = jnorth(i,j,k)
               IJKNE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

               itmp = iwest(i,j,k)
               IJKNW = FUNIJK(itmp,jnorth(itmp,j,k),k)

               ktmp = ktop(i,j,k)
               IJKTN = FUNIJK(i,jnorth(i,j,ktmp),ktmp)

               ktmp = kbot(i,j,k)
               IJKBN = FUNIJK(i,jnorth(i,j,ktmp),ktmp)

! Surface forces at i, j+1/2, k
! bulk viscosity term
! part of d/dy (tau_yy) xdxdydz =>
!         d/dy (lambda.trcD) xdxdydz =>
! delta (lambda.trcD)Ap|N-S  : at (i, j+1 - j-1, k)
               SBV = (LAMBDA_G(IJKN)*TRD_G(IJKN)-&
                      LAMBDA_G(IJK)*TRD_G(IJK))*AXZ(IJK)

! shear stress terms at i, j+1/2, k

! part of 1/x d/dx(x.tau_xy) xdxdydz =>
!         1/x d/dx (x.mu.du/dy) xdxdydz =>
! delta (x.mu.du/dy)Ayz |E-W : at (i+1/2 - i-1/2, j+1/2, k)
               SSX = AVG_Y_H(AVG_X_H(MU_G(IJK),MU_G(IJKE),I),&
                             AVG_X_H(MU_G(IJKN),MU_G(IJKNE),I),J)*&
                     (U_G(IJPK)-U_G(IJK))*ODY_N(J)*AYZ_V(IJK) - &
                     AVG_Y_H(AVG_X_H(MU_G(IJKW),MU_G(IJK),IM),&
                             AVG_X_H(MU_G(IJKNW),MU_G(IJKN),IM),J)*&
                     (U_G(IMJPK)-U_G(IMJK))*ODY_N(J)*AYZ_V(IMJK)

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dy) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at (i, j+1 - j-1, k)
               SSY = MU_G(IJKN)*(V_G(IJPK)-V_G(IJK))*ODY(JP)*&
                        AXZ_V(IJK) - &
                     MU_G(IJK)*(V_G(IJK)-V_G(IJMK))*ODY(J)*&
                        AXZ_V(IJMK)

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dy) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at (i, j+1/2, k+1/2 - k-1/2)
               SSZ = AVG_Y_H(AVG_Z_H(MU_G(IJK),MU_G(IJKT),K),&
                             AVG_Z_H(MU_G(IJKN),MU_G(IJKTN),K),J)*&
                     (W_G(IJPK)-W_G(IJK))*ODY_N(J)*AXY_V(IJK) - &
                     AVG_Y_H(AVG_Z_H(MU_G(IJKB),MU_G(IJK),KM),&
                             AVG_Z_H(MU_G(IJKBN),MU_G(IJKN),KM),J)*&
                     (W_G(IJPKM)-W_G(IJKM))*ODY_N(J)*AXY_V(IJKM)

! Add the terms
               lTAU_V_G(IJK) = SBV + SSX + SSY + SSZ

            ELSE
               lTAU_V_G(IJK) = ZERO
            ENDIF   ! end if (.NOT. IP_AT_N(IJK) .AND. EPGA>DIL_EP_S)
         ENDDO
         ENDDO
         ENDDO

      call send_recv(ltau_v_g,2)
      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE CALC_TAU_V_G
