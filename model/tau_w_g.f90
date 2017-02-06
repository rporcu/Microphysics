module calc_tau_w_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use geometry      , only: domlo, domhi

   contains
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
    subroutine calc_tau_w_g(slo,shi,lo,hi,&
                            ltau_w_g,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dx,dy,dz)

! Modules
!---------------------------------------------------------------------//

      use functions, only: avg, avg_h
      use param1   , only: zero
      use toleranc , only: dil_ep_s

      IMPLICIT NONE

      integer(c_int), intent(in ) :: slo(3),shi(3),lo(3),hi(3)

      ! tau_w_g
      real(c_real), INTENT(INOUT) :: trd_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(OUT) :: ltau_w_g&
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
      INTEGER :: I, J, K, IM, JM, KP
! Average volume fraction
      real(c_real) :: EPGA
! Source terms (Surface)
      real(c_real) :: Sbv, Ssx, Ssy, Ssz

      real(c_real) :: odz, axy, ayz, axz

      odz = 1.d0/dz
      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

!---------------------------------------------------------------------//
!     NOTE -- triply nested functions seem to break things -- hence the
!             use of the *tmp variables below
!---------------------------------------------------------------------//

        DO K = lo(3),hi(3)-1
         DO J = lo(2),hi(2)
          DO I = lo(1),hi(1)

            EPGA = AVG(EP_G(I,J,K),EP_G(i,j,k+1))

            IF (flag(i,j,k,4) > 1000 .AND. EPGA>DIL_EP_S) THEN

               IM = max(domlo(1)-1, i-1)
               JM = max(domlo(2)-1, j-1)
               KP = min(domhi(3)+1, k+1)

! Surface forces
! Bulk viscosity term
! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (lambda.trcD) xdxdydz=>
! delta (lambda.trcD)Ap |T-B : at (i, j, k+1 - k-1)
               SBV = (LAMBDA_G(i,j,k+1)*TRD_G(i,j,k+1)-&
                      LAMBDA_G(i,j,k  )*TRD_G(i,j,k  ))*AXY

! shear stress terms
! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>
!         1/x d/dx(mu.du/dz) xdxdydz =>
! delta (mu/x du/dz)Ayz |E-W : at (i+1/2-i-1/2, j, k+1/2)
               SSX = AVG_H(AVG_H(mu_g(i  ,j,k), &
                                 mu_g(i+1,j,k)), &
                           AVG_H(mu_g(i  ,j,k+1), &
                                 mu_g(i+1,j,k+1))) &
                    *(u_g(i,j,k+1)-u_g(i,j,k))*ODZ*AYZ &
                   - AVG_H(AVG_H(mu_g(i-1,j,k), &
                                 mu_g(i  ,j,k)), &
                           AVG_H(mu_g(i-1,j,k+1), &
                                 mu_g(i  ,j,k+1))) &
                     *(u_g(i-1,j,k+1)-u_g(i-1,j,k))*ODZ*AXZ
! DY(J)*HALF(DZ(k)+DZ(kp)) = oX_E(IM)*AYZ_W(IMJK), but avoids singularity

! part of d/dy (tau_zy) xdxdydz =>
!         d/dy (mu/x dv/dz) xdxdydz =>
! delta (mu/x dv/dz)Axz |N-S : at (i, j+1/2 - j-1/2, k+1/2)
               SSY = AVG_H(AVG_H(mu_g(i,j,k  ),mu_g(i,j+1,k  )),&
                           AVG_H(mu_g(i,j,k+1),mu_g(i,j+1,k+1))) &
                     *(v_g(i,j,k+1)-v_g(i,j,k))*ODZ*AXZ &
                       - AVG_H(AVG_H(mu_g(i,j-1,k  ),mu_g(i,j  ,k  )), &
                               AVG_H(mu_g(i,j-1,k+1),mu_g(i,j  ,k+1))) &
                         *(v_g(I,j-1,k+1)-v_g(i,j-1,k))*ODZ*AXZ

! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (mu/x dw/dz) xdxdydz =>
! delta (mu/x dw/dz)Axy |T-B : at (i, j, k+1 - k-1)
               SSZ = mu_g(i,j,k+1)*(w_g(i,j,k+1)-w_g(i,j,k  ))*ODZ*AXY - &
                     mu_g(i,j,k  )*(w_g(i,j,k  )-w_g(i,j,k-1))*ODZ*AXY

               ! Add the terms
               ltau_w_g(i,j,k) =  SBV + SSX + SSY + SSZ

            ELSE

               ltau_w_g(i,j,k) = ZERO

            ENDIF
         ENDDO
         ENDDO
         ENDDO

    end subroutine calc_tau_w_g

end module calc_tau_w_g_module
