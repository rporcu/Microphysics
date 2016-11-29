!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_U_g                                            C
!  Purpose: Cross terms in the gradient of stress in U_g momentum      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
!                                                                      C
!                                                                      C
!  Comments: This routine calculates the components of the gas phase   C
!  viscous stress tensor of the u-momentum equation that cannot be     C
!  cast in the form: mu.grad(u). These components are stored in the    C
!  passed variable, which is then used as a source of the u-momentum   C
!  equation.                                                           C
!                                                                      C
!  The following u component viscous stress tensor terms are           C
!  calculated here:                                                    C
!  > Combined part of 1/x d/dx (x.tau_xx) xdxdydz and                  C
!                              -tau_zz/x xdxdydz =>                    C
!    1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>     C
!        d/dx (lambda.trcD) xdxdydz =>                                 C
!    delta (lambda.trcD)Ap|E-W                                         C
!  > part of 1/x d/dx(x.tau_xx) xdxdydz =>                             C
!            1/x d/dx (x.mu.du/dx) xdxdydz =>                          C
!    delta (mu du/dx)Ayz |E-W                                          C
!  > part of d/dy (tau_xy) xdxdydz =>                                  C
!           d/dy (mu.dv/dx) xdxdydz =>                                 C
!    delta (mu.dv/dx)Axz |N-S                                          C
!  > part of 1/x d/dz (tau_xz) xdxdydz =>                              C
!            1/x d/dz (mu.dw/dx) xdxdydz =>                            C
!    delta (mu.dw/dx)Axy |T-B                                          C
!                                                                      C
!  To reconstitute the complete u-momentum gas phase viscous stress    C
!  tensor would require including the term calculated in source_u_g    C
!  and the 'diffusional' components (i.e., those of the form           C
!  mu.grad(u)                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_U_G(lTAU_U_G)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3
      USE param1, only: zero, half
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc, only: dil_ep_s
      USE geometry
      USE compar
      USE functions
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_U_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_U_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IP, JM, KM
      INTEGER :: IJK, IJKE, IJKN, IJKS, IJKT, IJKB
      INTEGER :: IJKNE, IJKSE, IJKTE, IJKBE
      INTEGER :: IPJK, IMJK, IJMK, IJKM
      INTEGER :: IPJMK, IPJKM
      INTEGER :: itmp, jtmp, ktmp
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average viscosity
      DOUBLE PRECISION :: MU_Ge, MU_gbe
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
            IJKE = FUNIJK(ieast(i,j,k),j,k)

            EPGA = AVG(EP_G(I,J,K),EP_G(ieast(i,j,k),j,k))

            IF (.NOT.ip_at_e(i,j,k) .AND. EPGA>DIL_EP_S) THEN

               IP = IP1(I)
               JM = JM1(J)
               KM = KM1(K)

               IPJK = FUNIJK(iplus(i,j,k),j,k)
               IMJK = FUNIJK(iminus(i,j,k),j,k)
               IJMK = FUNIJK(i,jminus(i,j,k),k)
               IJKM = FUNIJK(i,j,kminus(i,j,k))

               itmp = iplus(i,j,k)
               IPJMK = FUNIJK(itmp,jminus(itmp,j,k),k)


               IPJMK = FUNIJK(iplus(i,j,k),jminus(i,j,k),k)

               ktmp = kminus(i,j,k)
               IPJKM = FUNIJK(iplus(i,j,ktmp),j,ktmp)

               IJKE = FUNIJK(ieast(i,j,k),j,k)

               jtmp = jnorth(i,j,k)
               IJKN  = FUNIJK(i,jtmp,k)
               IJKNE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

               jtmp = jsouth(i,j,k)
               IJKS  = FUNIJK(i,jtmp,k)
               IJKSE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

               ktmp = ktop(i,j,k)
               IJKT  = FUNIJK(i,j,ktmp)
               IJKTE = FUNIJK(ieast(i,j,ktmp),j,ktmp)

               ktmp = kbot(i,j,k)
               IJKB  = FUNIJK(i,j,ktmp)
               IJKBE = FUNIJK(ieast(i,j,ktmp),j,ktmp)


! Surface forces at i+1/2, j, k
! bulk viscosity term
! combines part of 1/x d/dx (x.tau_xx) xdxdydz and -tau_zz/x xdxdydz =>
! combines 1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>
!              d/dx (lambda.trcD) xdxdydz
! delta (lambda.trcD)Ap |E-W : at (i+1 - i-1), j, k
               SBV = (LAMBDA_G(ieast(i,j,k),j,k)*TRD_G(ieast(i,j,k),j,k)-&
                      LAMBDA_G(i,j,k)*TRD_G(i,j,k))*AYZ

! shear stress terms at i+1/2, j, k
! part of 1/x d/dx(x.tau_xx) xdxdydz =>
!         1/x d/dx (x.mu.du/dx) xdxdydz =>
! delta (mu du/dx)Ayz |E-W : at (i+1 - i-1), j, k
               SSX = MU_G(ieast(i,j,k),j,k)*(U_G(iplus(i,j,k),j,k)-U_G(I,J,K))*ODX*AYZ - &
                  MU_G(i,j,k)*(U_G(I,J,K)-U_G(iminus(i,j,k),j,k))*ODX*AYZ

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dx) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at i+1/2, (j+1/2 - j-1/2), k
               SSY = AVG_H(AVG_H(MU_G(i,j,k),MU_G(i,jtmp,k)),&
                           AVG_H(MU_G(ieast(i,j,k),j,k),MU_G(ieast(i,jtmp,k),jtmp,k)))*&
                        (V_G(iplus(i,j,k),j,k)-V_G(I,J,K))*ODX*AXZ - &
                     AVG_H(AVG_H(MU_G(i,jtmp,k),MU_G(i,j,k)),&
                           AVG_H(MU_G(ieast(i,jtmp,k),jtmp,k),MU_G(ieast(i,j,k),j,k)))*&
                           (V_G(iplus(i,j,k),jminus(i,j,k),k)-&
                           V_G(i,jminus(i,j,k),k))*ODX*AXZ

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dx) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at i+1/2, j, (k+1/2 - k-1/2)
               MU_GE = AVG_H(AVG_H(MU_G(i,j,k),MU_G(i,j,ktmp)),&
                             AVG_H(MU_G(ieast(i,j,k),j,k),MU_G(ieast(i,j,ktmp),j,ktmp)))
               MU_GBE = AVG_H(AVG_H(MU_G(i,j,ktmp),MU_G(i,j,k)),&
                              AVG_H(MU_G(ieast(i,j,ktmp),j,ktmp),MU_G(ieast(i,j,k),j,k)))
               SSZ = MU_GE*(W_G(iplus(i,j,k),j,k)-W_G(I,J,K))*ODX*AXY - &
                  MU_GBE*(W_G(IPlus(i,j,k),J,KMinus(i,j,k))-&
                  W_G(I,J,KMinus(i,j,k)))*ODX*AXY


! Add the terms
               lTAU_U_G(IJK) = SBV + SSX + SSY + SSZ

            ELSE
               lTAU_U_G(IJK) = ZERO
            ENDIF   ! end if (.NOT.ip_at_e(i,j,k) .AND. EPGA>DIL_EP_S)
         ENDDO
         ENDDO
         ENDDO

      ! call send_recv(ltau_u_g,2)

      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE CALC_TAU_U_G
