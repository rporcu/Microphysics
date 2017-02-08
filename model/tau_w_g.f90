module calc_tau_w_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use geometry      , only: domlo, domhi

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_Tau_W_g                                            !
!  Purpose: Cross terms in the gradient of stress in W_g momentum      !
!                                                                      !
!  Comments: This routine calculates the components of the gas phase   !
!  viscous stress tensor of the w-momentum equation that cannot be     !
!  cast in the form: mu.grad(w). These components are stored in the    !
!  passed variable, which is then used as a source of the w-momentum   !
!  equation.                                                           !
!                                                                      !
!  The following w component viscous stress tensor terms are           !
!  calculated here:                                                    !
!  > part of 1/x d/dz (tau_zz) xdxdydz =>                              !
!            1/x d/dz (lambda.trcD) xdxdydz=>                          !
!    delta (lambda.trcD)Ap |T-B                                        !
!  > part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently        !
!    part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>              !
!            1/x d/dx(mu.du/dz) xdxdydz =>                             !
!    delta (mu/x du/dz)Ayz |E-W                                        !
!  > part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently        !
!    part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>              !
!            1/x d/dx(mu.du/dz) xdxdydz =>                             !
!    delta (mu/x du/dz)Ayz |E-W                                        !
!  > part of 1/x d/dz (tau_zz) xdxdydz =>                              !
!            1/x d/dz (mu/x dw/dz) xdxdydz =>                          !
!    delta (mu/x dw/dz)Axy |T-B                                        !
!                                                                      !
!  To reconstitute the complete w-momentum gas phase viscous stress    !
!  tensor requires including the terms calculated in source_w_g        !
!  and the 'diffusional' components (i.e., those of the form           !
!  mu.grad(w)                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_tau_w_g(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
       ltau_w_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)

! Modules
!---------------------------------------------------------------------//

      use functions, only: avg, avg_h

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      real(c_real)  , intent(in   ) :: dx,dy,dz

      real(c_real), intent(inout) :: ltau_w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

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
      integer :: i, j, k
! Source terms (Surface)
      real(c_real) :: Sbv, Ssx, Ssy, Ssz
      real(c_real) :: odz, axy, ayz, axz

      odz = 1.d0/dz
      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz


      do k = lo(3),hi(3)-1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

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

               ltau_w_g(i,j,k) =  SBV + SSX + SSY + SSZ

            enddo
         enddo
      enddo

    end subroutine calc_tau_w_g

end module calc_tau_w_g_module
