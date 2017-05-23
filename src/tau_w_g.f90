module calc_tau_w_g_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

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
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

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

      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

! Surface forces
! Bulk viscosity term
! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (lambda.trcD) xdxdydz=>
! delta (lambda.trcD)Ap |T-B : at (i, j, k+1 - k-1)
               sbv = (lambda_g(i,j,k  )*trd_g(i,j,k  )-&
                      lambda_g(i,j,k-1)*trd_g(i,j,k-1))*axy

! shear stress terms
! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>
!         1/x d/dx(mu.du/dz) xdxdydz =>
! delta (mu/x du/dz)Ayz |E-W : at (i+1/2-i-1/2, j, k+1/2)
               ssx = avg_h(avg_h(mu_g(i,j,k-1), mu_g(i+1,j,k-1)),    &
                           avg_h(mu_g(i,j,k  ), mu_g(i+1,j,k  )))*   &
                           (u_g(i+1,j,k) - u_g(i+1,j,k-1))*odz*ayz - &
                     avg_h(avg_h(mu_g(i-1,j,k-1), mu_g(i,j,k-1)),    &
                           avg_h(mu_g(i-1,j,k  ), mu_g(i,j,k  )))*   &
                           (u_g(i  ,j,k) - u_g(i  ,j,k-1))*odz*ayz

! part of d/dy (tau_zy) xdxdydz =>
!         d/dy (mu/x dv/dz) xdxdydz =>
! delta (mu/x dv/dz)Axz |N-S : at (i, j+1/2 - j-1/2, k+1/2)
               ssy = avg_h(avg_h(mu_g(i,j,k-1), mu_g(i,j+1,k-1)),    &
                           avg_h(mu_g(i,j,k  ), mu_g(i,j+1,k  )))*   &
                           (v_g(i,j+1,k) - v_g(i,j+1,k-1))*odz*axz - &
                     avg_h(avg_h(mu_g(i,j-1,k-1), mu_g(i,j,k-1)),    &
                           avg_h(mu_g(i,j-1,k  ), mu_g(i,j,k  )))*   &
                           (v_g(i,j  ,k) - v_g(i,j  ,k-1))*odz*axz

! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (mu/x dw/dz) xdxdydz =>
! delta (mu/x dw/dz)Axy |T-B : at (i, j, k+1 - k-1)
               ssz = mu_g(i,j,k  )*(w_g(i,j,k+1)-w_g(i,j,k  ))*odz*axy - &
                     mu_g(i,j,k-1)*(w_g(i,j,k  )-w_g(i,j,k-1))*odz*axy

               ltau_w_g(i,j,k) =  sbv + ssx + ssy + ssz

            enddo
         enddo
      enddo

    end subroutine calc_tau_w_g

end module calc_tau_w_g_module
