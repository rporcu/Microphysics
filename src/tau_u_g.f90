module calc_tau_u_g_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_Tau_U_g                                            !
!  Purpose: Cross terms in the gradient of stress in U_g momentum      !
!                                                                      !
!  Comments: This routine calculates the components of the gas phase   !
!  viscous stress tensor of the u-momentum equation that cannot be     !
!  cast in the form: mu.grad(u). These components are stored in the    !
!  passed variable, which is then used as a source of the u-momentum   !
!  equation.                                                           !
!                                                                      !
!  The following u component viscous stress tensor terms are           !
!  calculated here:                                                    !
!  > Combined part of 1/x d/dx (x.tau_xx) xdxdydz and                  !
!                              -tau_zz/x xdxdydz =>                    !
!    1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>     !
!        d/dx (lambda.trcD) xdxdydz =>                                 !
!    delta (lambda.trcD)Ap|E-W                                         !
!  > part of 1/x d/dx(x.tau_xx) xdxdydz =>                             !
!            1/x d/dx (x.mu.du/dx) xdxdydz =>                          !
!    delta (mu du/dx)Ayz |E-W                                          !
!  > part of d/dy (tau_xy) xdxdydz =>                                  !
!           d/dy (mu.dv/dx) xdxdydz =>                                 !
!    delta (mu.dv/dx)Axz |N-S                                          !
!  > part of 1/x d/dz (tau_xz) xdxdydz =>                              !
!            1/x d/dz (mu.dw/dx) xdxdydz =>                            !
!    delta (mu.dw/dx)Axy |T-B                                          !
!                                                                      !
!  To reconstitute the complete u-momentum gas phase viscous stress    !
!  tensor would require including the term calculated in source_u_g    !
!  and the 'diffusional' components (i.e., those of the form           !
!  mu.grad(u)                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_tau_u_g(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
         ltau_u_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)

! Modules
!---------------------------------------------------------------------//
      use functions, only: avg, avg_h

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      real(c_real)  , intent(in   ) :: dx,dy,dz

      real(c_real), intent(  out) :: ltau_u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):shi(3))

      real(c_real), intent(in   ) :: trd_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables
!---------------------------------------------------------------------//
! Indices
      integer :: I, J, K
! Source terms (Surface)
      real(c_real) :: Sbv, Ssx, Ssy, Ssz
      real(c_real) :: odx, axy, axz, ayz

      odx = 1.d0/dx
      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1)-1,hi(1)

! Surface forces at i+1/2, j, k
! bulk viscosity term
! combines part of 1/x d/dx (x.tau_xx) xdxdydz and -tau_zz/x xdxdydz =>
! combines 1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>
!              d/dx (lambda.trcD) xdxdydz
! delta (lambda.trcD)Ap |E-W : at (i+1 - i-1), j, k
               sbv = (lambda_g(i+1,j,k)*trd_g(i+1,j,k)-&
                      lambda_g(i  ,j,k)*trd_g(i,j,k))*ayz

! shear stress terms at i+1/2, j, k
! part of 1/x d/dx(x.tau_xx) xdxdydz =>
!         1/x d/dx (x.mu.du/dx) xdxdydz =>
! delta (mu du/dx)Ayz |E-W : at (i+1 - i-1), j, k
               ssx = mu_g(i+1,j,k)*(u_g(i+1,j,k)-u_g(i  ,j,k))*odx*ayz - &
                     mu_g(i  ,j,k)*(u_g(i  ,j,k)-u_g(i-1,j,k))*odx*ayz

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dx) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at i+1/2, (j+1/2 - j-1/2), k
               ssy = avg_h(avg_h(mu_g(i  ,j,k), mu_g(i,  j+1,k)),    &
                           avg_h(mu_g(i+1,j,k), mu_g(i+1,j+1,k)))*   &
                           (v_g(i+1,j  ,k)-v_g(i,j  ,k))*odx*axz -   &
                     avg_h(avg_h(mu_g(i,  j-1,k), mu_g(i,  j,k)),    &
                           avg_h(mu_g(i+1,j-1,k), mu_g(i+1,j,k)))*   &
                           (v_g(i+1,j-1,k)-v_g(i,j-1,k))*odx*axz

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dx) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at i+1/2, j, (k+1/2 - k-1/2)
               ssz = avg_h(avg_h(mu_g(i  ,j,k), mu_g(i  ,j,k+1)),    &
                           avg_h(mu_g(i+1,j,k), mu_g(i+1,j,k+1))) *  &
                           (w_g(i+1,j,k  )-w_g(i,j,k  ))*odx*axy -   &
                     avg_h(avg_h(mu_g(i  ,j,k-1), mu_g(i  ,j,k  )),  &
                           avg_h(mu_g(i+1,j,k-1), mu_g(i+1,j,k  )))* &
                           (w_g(i+1,j,k-1)-w_g(i,j,k-1))*odx*axy


               ltau_u_g(i,j,k) = SBV + SSX + SSY + SSZ

            enddo
         enddo
      enddo

      RETURN

    end subroutine calc_tau_u_g

end module calc_tau_u_g_module
