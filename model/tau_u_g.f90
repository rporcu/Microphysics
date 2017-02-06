module calc_tau_u_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use geometry      , only: domlo, domhi

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_Tau_U_g                                            C
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
      subroutine calc_tau_u_g(slo,shi,lo,hi,&
                              ltau_u_g,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dx,dy,dz)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero
      USE toleranc, only: dil_ep_s

      use functions, only: avg, avg_h

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
      real(c_real) :: mu_ge, MU_gbe
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

        DO K = lo(3),hi(3)
          DO J = lo(2),hi(2)
            DO I = lo(1),hi(1)-1

            EPGA = AVG(EP_G(I,J,K),EP_G(i+1,j,k))

            IF (flag(i,j,k,2) > 1000 .AND. EPGA>DIL_EP_S) THEN

               IP = min(domhi(1)+1, i+1)
               JM = max(domlo(2)-1, j-1)
               KM = max(domlo(3)-1, k-1)

! Surface forces at i+1/2, j, k
! bulk viscosity term
! combines part of 1/x d/dx (x.tau_xx) xdxdydz and -tau_zz/x xdxdydz =>
! combines 1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>
!              d/dx (lambda.trcD) xdxdydz
! delta (lambda.trcD)Ap |E-W : at (i+1 - i-1), j, k
               SBV = (LAMBDA_G(i+1,j,k)*TRD_G(i+1,j,k)-&
                      LAMBDA_G(i  ,j,k)*TRD_G(i,j,k))*AYZ

! shear stress terms at i+1/2, j, k
! part of 1/x d/dx(x.tau_xx) xdxdydz =>
!         1/x d/dx (x.mu.du/dx) xdxdydz =>
! delta (mu du/dx)Ayz |E-W : at (i+1 - i-1), j, k
               SSX = mu_g(i+1,j,k)*(u_g(i+1,j,k)-u_g(I,J,K))*odx*AYZ - &
                     mu_g(i  ,j,k)*(u_g(i  ,j,k)-u_g(i-1,j,k))*odx*AYZ

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dx) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at i+1/2, (j+1/2 - j-1/2), k
               SSY = AVG_H(AVG_H(mu_g(i,j,k), &
                                 mu_g(i,j+1,k)), &
                           AVG_H(mu_g(i+1,j,k), &
                                 mu_g(i+1,j+1,k))) &
                     *(v_g(i+1,j,k)-v_g(I,J,K))*odx*axz &
                     - AVG_H(AVG_H(mu_g(i,j-1,k), &
                                   mu_g(i,j,k)), &
                             AVG_H(mu_g(i+1,j-1,k), &
                                   mu_g(i+1,j  ,k))) &
                          *(v_g(i+1,j-1,k) &
                          - v_g(i  ,j-1,k))*odx*axz

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dx) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at i+1/2, j, (k+1/2 - k-1/2)
               mu_gE  = AVG_H(AVG_H(mu_g(i  ,j,k), &
                                    mu_g(i  ,j,k+1)), &
                              AVG_H(mu_g(i+1,j,k  ), &
                                    mu_g(i+1,j,k+1)))
               mu_gBE = AVG_H(AVG_H(mu_g(i  ,j,k-1), &
                                    mu_g(i  ,j,k  )), &
                              AVG_H(mu_g(i+1,j,k-1), &
                                    mu_g(i+1,j,k  )))
               SSZ = mu_gE *(w_g(i+1,j,k  )-w_g(i,j,k  ))*odx*axy - &
                     mu_gBE*(w_g(i+1,j,k-1)-w_g(i,j,k-1))*odx*axy


               ! Add the terms
               ltau_u_g(i,j,k) = SBV + SSX + SSY + SSZ

            ELSE

               ltau_u_g(i,j,k) = ZERO

            ENDIF
         ENDDO
         ENDDO
         ENDDO

      RETURN

    end subroutine calc_tau_u_g

end module calc_tau_u_g_module
