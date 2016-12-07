!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_W_g                                            C
!  Purpose: Cross terms in the gradient of stress in W_g momentum      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
!                                                                      C
!                                                                      C
!  Comments: This routine calculates the components of the gas phase   C
!  viscous stress tensor of the w-momentum equation that cannot be     C
!  cast in the form: mu.grad(w). These components are stored in the    C
!  passed variable, which is then used as a source of the w-momentum   C
!  equation.                                                           C
!                                                                      C
!  The following w component viscous stress tensor terms are           C
!  calculated here:                                                    C
!  > part of 1/x d/dz (tau_zz) xdxdydz =>                              C
!            1/x d/dz (lambda.trcD) xdxdydz=>                          C
!    delta (lambda.trcD)Ap |T-B                                        C
!  > part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently        C
!    part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>              C
!            1/x d/dx(mu.du/dz) xdxdydz =>                             C
!    delta (mu/x du/dz)Ayz |E-W                                        C
!  > part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently        C
!    part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>              C
!            1/x d/dx(mu.du/dz) xdxdydz =>                             C
!    delta (mu/x du/dz)Ayz |E-W                                        C
!  > part of 1/x d/dz (tau_zz) xdxdydz =>                              C
!            1/x d/dz (mu/x dw/dz) xdxdydz =>                          C
!    delta (mu/x dw/dz)Axy |T-B                                        C
!                                                                      C
!  To reconstitute the complete w-momentum gas phase viscous stress    C
!  tensor requires including the terms calculated in source_w_g        C
!  and the 'diffusional' components (i.e., those of the form           C
!  mu.grad(w)                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_W_G(lTAU_W_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE param1, only: zero
      USE toleranc, only: dil_ep_s
      USE geometry, only: axy, axz, ayz, odz, dy, dz

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_W_g
      DOUBLE PRECISION, INTENT(INOUT) :: trd_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)

      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: lambda_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JM, KP
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


            EPGA = AVG(EP_G(I,J,K),EP_G(i,j,ktop(i,j,k)))

            IF ( .NOT.ip_at_t(i,j,k) .AND. EPGA>DIL_EP_S) THEN

               IM = IM1(I)
               JM = JM1(J)
               KP = KP1(K)

! Surface forces
! Bulk viscosity term
! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (lambda.trcD) xdxdydz=>
! delta (lambda.trcD)Ap |T-B : at (i, j, k+1 - k-1)
               SBV = (LAMBDA_G(i,j,ktop(i,j,k))*TRD_G(i,j,ktop(i,j,k))-&
                      LAMBDA_G(i,j,k)*TRD_G(i,j,k))*AXY

! shear stress terms
! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>
!         1/x d/dx(mu.du/dz) xdxdydz =>
! delta (mu/x du/dz)Ayz |E-W : at (i+1/2-i-1/2, j, k+1/2)
               SSX = AVG_H(AVG_H(MU_G(i,j,k), &
                                 MU_G(ieast(i,j,k),j,k)), &
                           AVG_H(MU_G(i,j,ktop(i,j,k)), &
                                 MU_G(ieast(i,j,ktop(i,j,k)),j,ktop(i,j,k)))) &
                    *(U_G(I,J,KPlus(i,j,k))-U_G(I,J,K))*ODZ*AYZ &
                   - AVG_H(AVG_H(MU_G(iwest(i,j,k),j,k), &
                                 MU_G(i,j,k)), &
                           AVG_H(MU_G(iwest(i,j,ktop(i,j,k)),j,ktop(i,j,k)), &
                                 MU_G(i,j,ktop(i,j,k)))) &
                     *(U_G(IMinus(i,j,k),J,KPlus(i,j,k))-U_G(IMinus(i,j,k),J,K))*ODZ*DY*DZ
! DY(J)*HALF(DZ(k)+DZ(kp)) = oX_E(IM)*AYZ_W(IMJK), but avoids singularity

! part of d/dy (tau_zy) xdxdydz =>
!         d/dy (mu/x dv/dz) xdxdydz =>
! delta (mu/x dv/dz)Axz |N-S : at (i, j+1/2 - j-1/2, k+1/2)
               SSY = AVG_H(AVG_H(MU_G(i,j,k),MU_G(i,jnorth(i,j,k),k)),&
                           AVG_H(MU_G(i,j,ktop(i,j,k)), &
                           MU_G(i,jnorth(i,j,k),ktop(i,jnorth(i,j,k),k)))) &
                     *(V_G(I,J,KPlus(i,j,k))-V_G(I,J,K))*ODZ*AXZ &
                       - AVG_H(AVG_H(MU_G(i,jsouth(i,j,k),k), &
                                     MU_G(i,j,k)), &
                               AVG_H(MU_G(i,jsouth(i,j,k),ktop(i,jsouth(i,j,k),k)), &
                                     MU_G(i,j,ktop(i,j,k)))) &
                         *(V_G(I,JMinus(i,j,k),KPlus(i,j,k))-V_G(I,JMinus(i,j,k),K))*ODZ*AXZ

! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (mu/x dw/dz) xdxdydz =>
! delta (mu/x dw/dz)Axy |T-B : at (i, j, k+1 - k-1)
               SSZ = MU_G(i,j,ktop(i,j,k))*(W_G(I,J,KPlus(i,j,k))-W_G(I,J,K))*ODZ*AXY - &
                     MU_G(i,j,k)*(W_G(I,J,K)-W_G(I,J,KMinus(i,j,k)))*ODZ*AXY

! Add the terms
               lTAU_W_G(i,j,k) =  SBV + SSX + SSY + SSZ


            ELSE
               lTAU_W_G(i,j,k) = ZERO
            ENDIF   ! end if (.NOT. ip_at_t(i,j,k) .AND. EPGA>DIL_EP_S)
         ENDDO
         ENDDO
         ENDDO


      ! call send_recv(ltau_w_g,2)

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

    END SUBROUTINE CALC_TAU_W_G
