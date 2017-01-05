MODULE CALC_TAU_V_G_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   contains

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
      SUBROUTINE CALC_TAU_V_G(lTAU_V_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dy)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE param1, only: zero
      USE toleranc, only: dil_ep_s
      USE geometry, only: axy, axz, ayz
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_V_g
      real(c_real), INTENT(INOUT) :: trd_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), INTENT(OUT) :: lTAU_V_g&
          (istart3:iend3,jstart3:jend3,kstart3:kend3)

      real(c_real), INTENT(IN   ) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), INTENT(IN   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), INTENT(IN   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), INTENT(IN   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), INTENT(IN   ) :: lambda_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), INTENT(IN   ) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
            (istart3:iend3,jstart3:jend3,kstart3:kend3,4)
      real(c_real), INTENT(IN   ) :: dy

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JP, KM
! Average volume fraction
      real(c_real) :: EPGA
! Source terms (Surface)
      real(c_real) :: Sbv, Ssx, Ssy, Ssz

      real(c_real) :: ody

      ody = 1.d0/dy

!---------------------------------------------------------------------//
!     NOTE -- triply nested functions seem to break things -- hence the
!             use of the *tmp variables below
!---------------------------------------------------------------------//

        DO K = kstart3, kend3
         DO J = jstart3, jend3
          DO I = istart3, iend3

            EPGA = AVG(EP_G(I,J,K),EP_G(i,jnorth(i,j,k),k))
            IF (flag(i,j,k,3) > 1000 .AND. EPGA>DIL_EP_S) THEN

               JP = JP1(J)
               IM = IM1(I)
               KM = KM1(K)

! Surface forces at i, j+1/2, k
! bulk viscosity term
! part of d/dy (tau_yy) xdxdydz =>
!         d/dy (lambda.trcD) xdxdydz =>
! delta (lambda.trcD)Ap|N-S  : at (i, j+1 - j-1, k)
               SBV = (LAMBDA_G(i,jnorth(i,j,k),k)*TRD_G(i,jnorth(i,j,k),k)-&
                      LAMBDA_G(i,j,k)*TRD_G(i,j,k))*AXZ

! shear stress terms at i, j+1/2, k

! part of 1/x d/dx(x.tau_xy) xdxdydz =>
!         1/x d/dx (x.mu.du/dy) xdxdydz =>
! delta (x.mu.du/dy)Ayz |E-W : at (i+1/2 - i-1/2, j+1/2, k)
               SSX = AVG_H(AVG_H(MU_G(i,j,k), &
                                 MU_G(ieast(i,j,k),j,k)),&
                           AVG_H(MU_G(i,jnorth(i,j,k),k), &
                                 MU_G(ieast(i,jnorth(i,j,k),k),jnorth(i,j,k),k))) &
                    *(U_G(I,JPlus(i,j,k),K)-U_G(I,J,K))*ODY*AYZ &
                    - AVG_H(AVG_H(MU_G(iwest(i,j,k),j,k), &
                                  MU_G(i,j,k)), &
                            AVG_H(MU_G(iwest(i,j,k),jnorth(iwest(i,j,k),j,k),k), &
                                  MU_G(i,jnorth(i,j,k),k))) &
                    *(U_G(IMinus(i,j,k),JPlus(i,j,k),K)-U_G(IMinus(i,j,k),J,K))*ODY*AYZ

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dy) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at (i, j+1 - j-1, k)
               SSY = MU_G(i,jnorth(i,j,k),k)*(V_G(I,JPlus(i,j,k),K)-V_G(I,J,K))*ODY*AXZ - &
                     MU_G(i,j,k)*(V_G(I,J,K)-V_G(I,JMinus(i,j,k),K))*ODY*AXZ

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dy) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at (i, j+1/2, k+1/2 - k-1/2)
               SSZ = AVG_H(AVG_H(MU_G(i,j,k), &
                                 MU_G(i,j,ktop(i,j,k))), &
                           AVG_H(MU_G(i,jnorth(i,j,k),k), &
                                 MU_G(i,jnorth(i,j,ktop(i,j,k)),ktop(i,j,k)))) &
                    *(W_G(I,JPlus(i,j,k),K)-W_G(I,J,K))*ODY*AXY &
                    - AVG_H(AVG_H(MU_G(i,j,kbot(i,j,k)), &
                                  MU_G(i,j,k)), &
                            AVG_H(MU_G(i,jnorth(i,j,kbot(i,j,k)),kbot(i,j,k)), &
                                  MU_G(i,jnorth(i,j,k),k))) &
                    *(W_G(I,JPlus(i,j,k),KMinus(i,j,k))-W_G(I,J,KMinus(i,j,k)))*ODY*AXY

! Add the terms
               lTAU_V_G(i,j,k) = SBV + SSX + SSY + SSZ

            ELSE
               lTAU_V_G(i,j,k) = ZERO
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ! call send_recv(ltau_v_g,2)
      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

    END SUBROUTINE CALC_TAU_V_G
END MODULE CALC_TAU_V_G_MODULE
