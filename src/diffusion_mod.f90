! 
!              
!  This module contains the subroutines to compute the three components
!  of the diffusion term div(tau) where
!
!      tau = mu ( grad(u) + grad(u)^T ) + div ( lambda div(u) I )
!  
!  Author: Michele Rosso
! 
!  Date: October 12, 2017
!
! 
module diffusion_mod
   
   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   
   implicit none
   private

contains

   
   !
   ! Computes d(txx)/dx + d(txy)/dy + d(txz)/dz at the
   ! x-edges ( u-component location ), where:
   !
   !  txx = 2*mu*(du/dx) + lambda div(u)
   !  txy =  mu * ( du/dy + dv/dx ) 
   !  txz =  mu * ( du/dz + dw/dx )
   ! 
   subroutine compute_divtau_x ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, mu, slo, shi, lambda, divu, divtau_x, dx )    &
        & bind(C, name="compute_divtau_x")

      integer(c_int),  intent(in)    :: lo(3),  hi(3)   ! Tile indeces for mf associated to ug
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in)    :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      real(ar),        intent(in   ) :: dx(3)
      
      real(ar),        intent(inout) ::                           &
           & divtau_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), & 
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)),       &
           & mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),       &
           & lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),   &
           & divu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
           
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: txx_e, txx_w
      real(ar)                       :: txy_n, txy_s, mu_n, mu_s
      real(ar)                       :: txz_t, txz_b, mu_t, mu_b
      real(ar)                       :: dudy, dudz, dvdx, dwdx
      real(ar), parameter            :: i4  = one / 4.0_ar
      real(ar), parameter            :: two = 2.0_ar
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! Compute txx_e
               txx_e = two * mu(i,j,k)   * ( ug(i+1,j,k) - ug(i,j,k)   ) * idx + & 
                    & lambda(i,j,k)   * divu(i,j,k)


               ! Compute txx_w
               txx_w = two * mu(i-1,j,k) * ( ug(i,j,k)   - ug(i-1,j,k) ) * idx + &
                    & lambda(i-1,j,k) * divu(i-1,j,k)
               
               ! Compute txy_n
               mu_n = ( mu(i,j,k) + mu(i-1,j,k) + mu(i,j+1,k) + mu(i-1,j+1,k) ) * i4
               
               dudy = ( ug(i,j+1,k) - ug(i,j,k) ) * idy
               dvdx = ( vg(i,j+1,k) - vg(i-1,j+1,k) ) * idx
               
               txy_n = mu_n * ( dudy + dvdx ) 

               ! Compute txy_s
               mu_s = ( mu(i,j-1,k) + mu(i-1,j-1,k) + mu(i,j,k) + mu(i-1,j,k) ) * i4
               
               dudy = ( ug(i,j,k) - ug(i,j-1,k) ) * idy
               dvdx = ( vg(i,j,k) - vg(i-1,j,k) ) * idx
               
               txy_s = mu_s * ( dudy + dvdx )
               
               ! Compute txz_t
               mu_t = ( mu(i,j,k) + mu(i-1,j,k) + mu(i,j,k+1) + mu(i-1,j,k+1) ) * i4

               dudz = ( ug(i,j,k+1) - ug(i,j,k) ) * idz
               dwdx = ( wg(i,j,k+1) - wg(i-1,j,k+1) ) * idx

               txz_t = mu_t * ( dudz + dwdx )

               ! Compute txz_b
               mu_b = ( mu(i,j,k-1) + mu(i-1,j,k-1) + mu(i,j,k) + mu(i-1,j,k) ) * i4

               dudz = ( ug(i,j,k) - ug(i,j,k-1) ) * idz
               dwdx = ( wg(i,j,k) - wg(i-1,j,k) ) * idx

               txz_b = mu_b * ( dudz + dwdx )

               ! Assemble divtau_x
               divtau_x(i,j,k) = ( txx_e - txx_w ) * idx + &
                    &            ( txy_n - txy_s ) * idy + &
                    &            ( txz_t - txz_b ) * idz

            end do
         end do
      end do

   end subroutine compute_divtau_x



   
   !
   ! Computes d(tyx)/dx + d(tyy)/dy + d(tyz)/dz at the
   ! y-edges ( v-component location ), where:
   !
   !  tyx =  mu * ( du/dy + dv/dx ) 
   !  tyy =  2*mu*(dv/dy) + lambda div(u)
   !  txz =  mu * ( dv/dz + dw/dy )
   ! 
   subroutine compute_divtau_y ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, mu, slo, shi, lambda, divu, divtau_y, dx )    &
        & bind(C, name="compute_divtau_y")

      integer(c_int),  intent(in)    :: lo(3),  hi(3)   ! Tile indeces for mf associated to ug
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in)    :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      real(ar),        intent(in   ) :: dx(3)
      
      real(ar),        intent(inout) ::                           &
           & divtau_y(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), & 
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)),       &
           & mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),       &
           & lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),   &
           & divu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
           
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: tyx_e, tyx_w, mu_e, mu_w
      real(ar)                       :: tyy_n, tyy_s
      real(ar)                       :: tyz_t, tyz_b, mu_t, mu_b
      real(ar)                       :: dudy, dvdx, dvdz, dwdy
      real(ar), parameter            :: i4  = one / 4.0_ar
      real(ar), parameter            :: two = 2.0_ar
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               
               ! Compute tyx_e
               mu_e = ( mu(i,j,k) + mu(i,j-1,k) + mu(i+1,j,k) + mu(i+1,j-1,k) ) * i4
               
               dudy = ( ug(i+1,j,k) - ug(i+1,j-1,k) ) * idy
               dvdx = ( vg(i+1,j,k) - vg(i,j,k)     ) * idx
               
               tyx_e = mu_e * ( dudy + dvdx ) 

               ! Compute txy_s
               mu_w = ( mu(i-1,j,k) + mu(i-1,j-1,k) + mu(i,j,k) + mu(i,j-1,k) ) * i4
               
               dudy = ( ug(i,j,k) - ug(i,j-1,k) ) * idy
               dvdx = ( vg(i,j,k) - vg(i-1,j,k) ) * idx
               
               tyx_w = mu_w * ( dudy + dvdx )

               ! Compute tyy_n
               tyy_n = two * mu(i,j,k)   * ( vg(i,j+1,k) - vg(i,j,k)   ) * idy + & 
                    & lambda(i,j,k)   * divu(i,j,k)

               ! Compute tyy_s
               tyy_s = two * mu(i,j-1,k) * ( vg(i,j,k)   - vg(i,j-1,k) ) * idy + &
                    & lambda(i,j-1,k) * divu(i,j-1,k)

               
               ! Compute tyz_t
               mu_t = ( mu(i,j,k) + mu(i,j-1,k) + mu(i,j,k+1) + mu(i,j-1,k+1) ) * i4

               dvdz = ( vg(i,j,k+1) - vg(i,j,k) )     * idz
               dwdy = ( wg(i,j,k+1) - wg(i,j-1,k+1) ) * idy

               tyz_t = mu_t * ( dvdz + dwdy )

               ! Compute tyz_b
               mu_b = ( mu(i,j,k-1) + mu(i,j-1,k-1) + mu(i,j,k) + mu(i,j-1,k) ) * i4

               dvdz = ( vg(i,j,k) - vg(i,j,k-1) ) * idz
               dwdy = ( wg(i,j,k) - wg(i,j-1,k) ) * idy

               tyz_b = mu_b * ( dvdz + dwdy )

               ! Assemble divtau_x
               divtau_y(i,j,k) = ( tyx_e - tyx_w ) * idx + &
                    &            ( tyy_n - tyy_s ) * idy + &
                    &            ( tyz_t - tyz_b ) * idz

            end do
         end do
      end do

   end subroutine compute_divtau_y




   
   
   !
   ! Computes d(tzx)/dx + d(tzy)/dy + d(tzz)/dz at the
   ! y-edges ( v-component location ), where:
   !
   !  tzx =  mu * ( du/dz + dw/dx ) 
   !  tzy =  mu * ( dv/dz + dw/dy )
   !  tzz =  2*mu*(dw/dy) + lambda div(u)
   ! 
   subroutine compute_divtau_z ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, mu, slo, shi, lambda, divu, divtau_z, dx )    &
        & bind(C, name="compute_divtau_z")

      integer(c_int),  intent(in)    :: lo(3),  hi(3)   ! Tile indeces for mf associated to ug
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in)    :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      real(ar),        intent(in   ) :: dx(3)
      
      real(ar),        intent(inout) ::                           &
           & divtau_z(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), & 
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)),       &
           & mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),       &
           & lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),   &
           & divu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
           
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: tzx_e, tzx_w, mu_e, mu_w
      real(ar)                       :: tzy_n, tzy_s, mu_n, mu_s
      real(ar)                       :: tzz_t, tzz_b
      real(ar)                       :: dudz, dvdz, dwdx, dwdy
      real(ar), parameter            :: i4  = one / 4.0_ar
      real(ar), parameter            :: two = 2.0_ar
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               
               ! Compute tzx_e
               mu_e = ( mu(i,j,k) + mu(i,j,k-1) + mu(i+1,j,k) + mu(i+1,j,k-1) ) * i4
               
               dudz = ( ug(i+1,j,k) - ug(i+1,j,k-1) ) * idz
               dwdx = ( wg(i+1,j,k) - wg(i,j,k)     ) * idx
               
               tzx_e = mu_e * ( dudz + dwdx ) 

               ! Compute tzx_w
               mu_w = ( mu(i-1,j,k) + mu(i-1,j,k-1) + mu(i,j,k) + mu(i,j,k-1) ) * i4
               
               dudz = ( ug(i,j,k) - ug(i,j,k-1) ) * idz
               dwdx = ( wg(i,j,k) - wg(i-1,j,k) ) * idx
               
               tzx_w = mu_w * ( dudz + dwdx )
               
               ! Compute tzy_n
               mu_n = ( mu(i,j,k) + mu(i,j,k-1) + mu(i,j+1,k) + mu(i,j+1,k-1) ) * i4

               dvdz = ( vg(i,j+1,k) - vg(i,j+1,k-1) ) * idz
               dwdy = ( wg(i,j+1,k) - wg(i,j,k) )     * idy

               tzy_n = mu_n * ( dvdz + dwdy )

               ! Compute tzy_s
               mu_s = ( mu(i,j-1,k) + mu(i,j-1,k-1) + mu(i,j,k) + mu(i,j,k-1) ) * i4

               dvdz = ( vg(i,j,k) - vg(i,j,k-1) ) * idz
               dwdy = ( wg(i,j,k) - wg(i,j-1,k) ) * idy

               tzy_s = mu_s * ( dvdz + dwdy )


               ! Compute tzz_t
               tzz_t = two * mu(i,j,k)   * ( wg(i,j,k+1) - wg(i,j,k)   ) * idz + & 
                    & lambda(i,j,k)   * divu(i,j,k)

               ! Compute tzz_s
               tzz_b = two * mu(i,j,k-1) * ( wg(i,j,k)   - wg(i,j,k-1) ) * idz + &
                    & lambda(i,j,k-1) * divu(i,j,k-1)
              
               ! Assemble divtau_x
               divtau_z(i,j,k) = ( tzx_e - tzx_w ) * idx + &
                    &            ( tzy_n - tzy_s ) * idy + &
                    &            ( tzz_t - tzz_b ) * idz

            end do
         end do
      end do

   end subroutine compute_divtau_z


   
end module diffusion_mod
