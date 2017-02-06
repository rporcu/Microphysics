module calc_tau_v_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use geometry      , only: domlo, domhi

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
      subroutine calc_tau_v_g(slo,shi,lo,hi,&
                              lTAU_v_g,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dx,dy,dz)

      use functions, only: avg, avg_h 
      use toleranc , only: dil_ep_s
      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! TAU_V_g
      real(c_real), INTENT(INOUT) :: trd_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(OUT) :: lTAU_V_g&
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
      INTEGER :: i, j, k, im, jp, km
! Average volume fraction
      real(c_real) :: EPGA
! Source terms (Surface)
      real(c_real) :: Sbv, Ssx, Ssy, Ssz

      real(c_real) :: ody, axy, axz, ayz

      ody = 1.d0/dy
      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

!---------------------------------------------------------------------//
!     NOTE -- triply nested functions seem to break things -- hence the
!             use of the *tmp variables below
!---------------------------------------------------------------------//

        DO k = lo(3),hi(3)
         DO j = lo(2),hi(2)-1
          DO i = lo(1),hi(1)

            EPGA = AVG(EP_G(i,j,k),EP_G(i,j+1,k))
            IF (flag(i,j,k,3) > 1000 .AND. EPGA>DIL_EP_S) THEN

               im = max(domlo(1)-1, i-1)
               jp = min(domhi(2)+1, j+1)
               km = max(domlo(3)-1, k-1)

! Surface forces at i, j+1/2, k
! bulk viscosity term
! part of d/dy (tau_yy) xdxdydz =>
!         d/dy (lambda.trcD) xdxdydz =>
! delta (lambda.trcD)Ap|N-S  : at (i, j+1 - j-1, k)
               SBV = (LAMBDA_G(i,j+1,k)*TRD_G(i,j+1,k)-&
                      LAMBDA_G(i,j  ,k)*TRD_G(i,j  ,k))*AXZ

! shear stress terms at i, j+1/2, k

! part of 1/x d/dx(x.tau_xy) xdxdydz =>
!         1/x d/dx (x.mu.du/dy) xdxdydz =>
! delta (x.mu.du/dy)Ayz |E-W : at (i+1/2 - i-1/2, j+1/2, k)
               SSX = AVG_H(AVG_H(mu_g(i  ,j  ,k), &
                                 mu_g(i+1,j  ,k)),&
                           AVG_H(mu_g(i  ,j+1,k), &
                                 mu_g(i+1,j+1,k))) &
                    *(u_g(i,j+1,k)-u_g(i,j,k))*ODY*AYZ &
                    - AVG_H(AVG_H(mu_g(i-1,j  ,k), &
                                  mu_g(i  ,j  ,k)), &
                            AVG_H(mu_g(i-1,j+1,k), &
                                  mu_g(i  ,j+1,k))) &
                    *(u_g(i-1,j+1,k)-u_g(i-1,j,k))*ODY*AYZ

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dy) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at (i, j+1 - j-1, k)
               SSY = mu_g(i,j+1,k)*(v_g(i,j+1,k)-v_g(i,j,k))*ODY*AXZ - &
                     mu_g(i,j,k)*(v_g(i,j,k)-v_g(i,j-1,k))*ODY*AXZ

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dy) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at (i, j+1/2, k+1/2 - k-1/2)
               SSZ = AVG_H(AVG_H(mu_g(i,j  ,k), &
                                 mu_g(i,j  ,k+1)), &
                           AVG_H(mu_g(i,j+1,k  ), &
                                 mu_g(i,j+1,k+1))) &
                    *(w_g(i,j+1,k)-w_g(i,j,k))*ODY*AXY &
                    - AVG_H(AVG_H(mu_g(i,j  ,k-1), &
                                  mu_g(i,j  ,k  )), &
                            AVG_H(mu_g(i,j+1,k-1), &
                                  mu_g(i,j+1,k  ))) &
                    *(w_g(i,j+1,k-1)-w_g(i,j,k-1))*ODY*AXY

               ! Add the terms
               lTAU_v_g(i,j,k) = SBV + SSX + SSY + SSZ

            ELSE

               lTAU_v_g(i,j,k) = 0.d0

            ENDIF
         ENDDO
         ENDDO
         ENDDO

    end subroutine calc_tau_v_g

end module calc_tau_v_g_module
