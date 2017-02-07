module calc_tau_v_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use geometry      , only: domlo, domhi

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_Tau_V_g                                            !
!  Purpose: Cross terms in the gradient of stress in V_g momentum      !
!                                                                      !
!  Comments: This routine calculates the components of the gas phase   !
!  viscous stress tensor of the v-momentum equation that cannot be     !
!  cast in the form: mu.grad(v). These components are stored in the    !
!  passed variable, which is then used as a source of the v-momentum   !
!  equation.                                                           !
!                                                                      !
!  The following v component viscous stress tensor terms are           !
!  calculated here:                                                    !
!  > part of d/dy (tau_yy) xdxdydz =>                                  !
!            d/dy (lambda.trcD) xdxdydz =>                             !
!    delta (lambda.trcD)Ap|N-S                                         !
!  > part of 1/x d/dx(x.tau_xy) xdxdydz =>                             !
!            1/x d/dx (x.mu.du/dy) xdxdydz =>                          !
!    delta (x.mu.du/dy)Ayz |E-W                                        !
!  > part of d/dy (tau_xy) xdxdydz =>                                  !
!           d/dy (mu.dv/dy) xdxdydz =>                                 !
!    delta (mu.dv/dx)Axz |N-S                                          !
!  > part of 1/x d/dz (tau_xz) xdxdydz =>                              !
!            1/x d/dz (mu.dw/dy) xdxdydz =>                            !
!    delta (mu.dw/dx)Axy |T-B                                          !
!                                                                      !
!  To reconstitute the full v-momentum gas phase viscous stress        !
!  tensor would require including the the 'diffusional' components     !
!  (i.e., those of the form mu.grad(v)                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_tau_v_g(slo,shi,lo,hi,&
         lTAU_v_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)

      use functions, only: avg, avg_h
      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) ::  lo(3), hi(3)
      real(c_real), intent(in   ) :: dx,dy,dz

      ! TAU_V_g
      real(c_real), intent(  out) :: ltau_v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: trd_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables
!---------------------------------------------------------------------//
! indices
      integer :: i, j, k
! Source terms (Surface)
      real(c_real) :: Sbv, Ssx, Ssy, Ssz
      real(c_real) :: ody, axy, axz, ayz

      ody = 1.d0/dy
      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)-1
            do i = lo(1),hi(1)

! Surface forces at i, j+1/2, k
! bulk viscosity term
! part of d/dy (tau_yy) xdxdydz =>
!         d/dy (lambda.trcD) xdxdydz =>
! delta (lambda.trcD)Ap|N-S  : at (i, j+1 - j-1, k)
               sbv = (lambda_g(i,j+1,k)*trd_g(i,j+1,k)-&
                      lambda_g(i,j  ,k)*trd_g(i,j  ,k))*axz

! shear stress terms at i, j+1/2, k

! part of 1/x d/dx(x.tau_xy) xdxdydz =>
!         1/x d/dx (x.mu.du/dy) xdxdydz =>
! delta (x.mu.du/dy)Ayz |E-W : at (i+1/2 - i-1/2, j+1/2, k)
               ssx = avg_h(avg_h(mu_g(i  ,j  ,k), &
                                 mu_g(i+1,j  ,k)),&
                           avg_h(mu_g(i  ,j+1,k), &
                                 mu_g(i+1,j+1,k))) &
                    *(u_g(i,j+1,k)-u_g(i,j,k))*ody*ayz &
                    - avg_h(avg_h(mu_g(i-1,j  ,k), &
                                  mu_g(i  ,j  ,k)), &
                            avg_h(mu_g(i-1,j+1,k), &
                                  mu_g(i  ,j+1,k))) &
                    *(u_g(i-1,j+1,k)-u_g(i-1,j,k))*ody*ayz

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

               ltau_v_g(i,j,k) = sbv + ssx + ssy + ssz

            enddo
         enddo
      enddo

    end subroutine calc_tau_v_g

end module calc_tau_v_g_module
