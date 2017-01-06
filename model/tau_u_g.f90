MODULE CALC_TAU_U_G_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
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
      SUBROUTINE CALC_TAU_U_G(slo,shi,lo,hi,&
                              lTAU_U_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dx,dy,dz)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero
      USE toleranc, only: dil_ep_s
      IMPLICIT NONE

      integer(c_int), intent(in ) :: slo(3),shi(3),lo(3),hi(3)

      ! TAU_U_g
      real(c_real), INTENT(INOUT) :: trd_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(OUT) :: lTAU_U_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), INTENT(IN   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), INTENT(IN   ) :: dx,dy,dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IP, JM, KM
! Average volume fraction
      real(c_real) :: EPGA
! Average viscosity
      real(c_real) :: MU_Ge, MU_gbe
! Source terms (Surface)
      real(c_real) :: Sbv, Ssx, Ssy, Ssz

      real(c_real) :: odx, axy, axz, ayz

      odx = 1.d0/dx
      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz
!---------------------------------------------------------------------//
!     NOTE -- triply nested functions seem to break things -- hence the
!             use of the *tmp variables below
!---------------------------------------------------------------------//

        DO K = slo(3),shi(3)
          DO J = slo(2),shi(2)
            DO I = slo(1),shi(1)

            EPGA = AVG(EP_G(I,J,K),EP_G(ieast(i,j,k),j,k))

            IF (flag(i,j,k,2) > 1000 .AND. EPGA>DIL_EP_S) THEN

               IP = IP1(I)
               JM = JM1(J)
               KM = KM1(K)

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
               SSY = AVG_H(AVG_H(MU_G(i,j,k), &
                                 MU_G(i,jnorth(i,j,k),k)), &
                           AVG_H(MU_G(ieast(i,j,k),j,k), &
                                 MU_G(ieast(i,jnorth(i,j,k),k),jnorth(i,j,k),k))) &
                     *(V_G(iplus(i,j,k),j,k)-V_G(I,J,K))*ODX*AXZ &
                     - AVG_H(AVG_H(MU_G(i,jsouth(i,j,k),k), &
                                   MU_G(i,j,k)), &
                             AVG_H(MU_G(ieast(i,jsouth(i,j,k),k),jsouth(i,j,k),k), &
                                   MU_G(ieast(i,j,k),j,k))) &
                          *(V_G(iplus(i,j,k),jminus(i,j,k),k) &
                         - V_G(i,jminus(i,j,k),k))*ODX*AXZ

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dx) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at i+1/2, j, (k+1/2 - k-1/2)
               MU_GE = AVG_H(AVG_H(MU_G(i,j,k), &
                                   MU_G(i,j,ktop(i,j,k))), &
                             AVG_H(MU_G(ieast(i,j,k),j,k), &
                                   MU_G(ieast(i,j,ktop(i,j,k)),j,ktop(i,j,k))))
               MU_GBE = AVG_H(AVG_H(MU_G(i,j,kbot(i,j,k)), &
                                    MU_G(i,j,k)), &
                              AVG_H(MU_G(ieast(i,j,kbot(i,j,k)),j,kbot(i,j,k)), &
                                    MU_G(ieast(i,j,k),j,k)))
               SSZ = MU_GE*(W_G(iplus(i,j,k),j,k)-W_G(I,J,K))*ODX*AXY - &
                  MU_GBE*(W_G(IPlus(i,j,k),J,KMinus(i,j,k))-&
                  W_G(I,J,KMinus(i,j,k)))*ODX*AXY


! Add the terms
               lTAU_U_G(i,j,k) = SBV + SSX + SSY + SSZ

            ELSE
               lTAU_U_G(i,j,k) = ZERO
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ! call send_recv(ltau_u_g,2)

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

    END SUBROUTINE CALC_TAU_U_G
END MODULE CALC_TAU_U_G_MODULE
