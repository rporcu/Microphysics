! 
!              
!  This module contains the subroutines to compute the three components
!  of the convetion term (U.grad)U
!  
!  Author: Michele Rosso
! 
!  Date: October 25, 2017
!
! 
module convection_mod

   
   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   
   implicit none
   private

contains

   !
   ! Computes udu/dx + vdu/dy + wdu/dz at the
   ! x-edges ( u-component location )
   ! 
   subroutine compute_ugradu_x ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, ugradu_x, dx )  bind(C, name="compute_ugradu_x")

      implicit none

      integer(c_int),  intent(in)    :: lo(3),  hi(3)   ! Tile indeces for mf associated to ug
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in)    :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)
      
      real(ar),        intent(inout) ::                           &
           & ugradu_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), & 
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
           
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: ududx, vdudy, wdudz
      real(ar)                       :: u_e, u_w, dudy, dudz
      real(ar)                       :: v_n, v_s
      real(ar)                       :: w_t, w_b
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! udu/dx
               u_e = edge_velocity ( ug(i,j,k),   ug(i+1,j,k) )
               u_w = edge_velocity ( ug(i-1,j,k), ug(i,j,k)   )

               ududx = ug(i,j,k) * ( u_e - u_w ) * idx

               ! vdu/dy
               v_n = half * ( vg(i,j+1,k) + vg(i-1,j+1,k) )
               v_s = half * ( vg(i,j,k)   + vg(i-1,j,k)   )

               dudy  = half * ( ug(i,j+1,k) - ug(i,j-1,k) ) * idy
               
               vdudy = half * ( v_n + v_s ) * dudy

               ! wdu/dz
               w_t = half * ( wg(i,j,k+1) + wg(i-1,j,k+1) )
               w_b = half * ( wg(i,j,k)   + wg(i-1,j,k)   )

               dudz  = half * ( ug(i,j,k+1) - ug(i,j,k-1) ) * idz
               
               wdudz = half * ( w_t + w_b ) * dudz

               ! Assemble terms
               ugradu_x(i,j,k) = ududx + vdudy + wdudz
               
            end do
         end do
      end do

      
   end subroutine compute_ugradu_x



   !
   ! Computes udv/dx + vdv/dy + wdv/dz at the
   ! y-edges ( v-component location )
   ! 
   subroutine compute_ugradu_y ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, ugradu_y, dx )  bind(C, name="compute_ugradu_y")

      implicit none

      integer(c_int),  intent(in)    :: lo(3),  hi(3)   ! Tile indeces for mf associated to vg
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in)    :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)
      
      real(ar),        intent(inout) ::                           &
           & ugradu_y(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), & 
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
           
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: udvdx, vdvdy, wdvdz
      real(ar)                       :: u_e, u_w
      real(ar)                       :: v_n, v_s, dvdx, dvdz
      real(ar)                       :: w_t, w_b
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! udvdx
               u_e = half * ( ug(i+1,j,k) + ug(i+1,j-1,k) )
               u_w = half * ( ug(i,j,k)   + ug(i,j-1,k)   )

               dvdx  = half * ( vg(i+1,j,k) - vg(i-1,j,k) ) * idx
               
               udvdx = half * ( u_e + u_w ) * dvdx
               
               ! vdv/dy
               v_n = edge_velocity ( vg(i,j,k),   vg(i,j+1,k) )
               v_s = edge_velocity ( vg(i,j-1,k), vg(i,j,k)   )

               vdvdy = vg(i,j,k) * ( v_n - v_s ) * idy

               ! wdv/dz
               w_t = half * ( wg(i,j,k+1) + wg(i,j-1,k+1) )
               w_b = half * ( wg(i,j,k)   + wg(i,j-1,k)   )

               dvdz  = half * ( vg(i,j,k+1) - vg(i,j,k-1) ) * idz
               
               wdvdz = half * ( w_t + w_b ) * dvdz

               ! Assemble terms
               ugradu_y(i,j,k) = udvdx + vdvdy + wdvdz
               
            end do
         end do
      end do

      
   end subroutine compute_ugradu_y




   !
   ! Computes udw/dx + vdw/dy + wdw/dz at the
   ! z-edges ( w-component location )
   ! 
   subroutine compute_ugradu_z ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, ugradu_z, dx )  bind(C, name="compute_ugradu_z")

      implicit none

      integer(c_int),  intent(in)    :: lo(3),  hi(3)   ! Tile indeces for mf associated to vg
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in)    :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)
      
      real(ar),        intent(inout) ::                           &
           & ugradu_z(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), & 
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
           
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: udwdx, vdwdy, wdwdz
      real(ar)                       :: u_e, u_w
      real(ar)                       :: v_n, v_s
      real(ar)                       :: w_t, w_b, dwdx, dwdy
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! udwdx
               u_e = half * ( ug(i+1,j,k) + ug(i+1,j,k-1) )
               u_w = half * ( ug(i,j,k)   + ug(i,j,k-1)   )

               dwdx  = half * ( wg(i+1,j,k) - wg(i-1,j,k) ) * idx
               
               udwdx = half * ( u_e + u_w ) * dwdx

               ! vdw/dy
               v_n = half * ( vg(i,j+1,k) + vg(i,j+1,k-1) )
               v_s = half * ( vg(i,j,k)   + vg(i,j,k-1)   )

               dwdy  = half * ( wg(i,j+1,k) - wg(i,j-1,k) ) * idy
               
               vdwdy = half * ( v_n + v_s ) * dwdy
              
               ! wdw/dz
               w_t = edge_velocity ( wg(i,j,k),   wg(i,j,k+1) )
               w_b = edge_velocity ( wg(i,j,k-1), wg(i,j,k)   )

               wdwdz = wg(i,j,k) * ( w_t - w_b ) * idz
               
               ! Assemble terms
               ugradu_z(i,j,k) = udwdx + vdwdy + wdwdz
               
            end do
         end do
      end do

      
   end subroutine compute_ugradu_z

   
   ! Upwind along the same direction as the velocity component
   function edge_velocity (u_minus, u_plus) result (ev)
      
      real(ar), intent(in) :: u_minus, u_plus
      real(ar)             :: ev
      
      if ( ( u_minus >= zero ) .and. ( u_plus >= zero ) ) then
         ev = u_minus
      else if ( ( u_minus <= zero ) .and. ( u_plus <= zero ) ) then
         ev = u_plus
      else if ( ( u_minus > zero ) .and. ( u_plus < zero ) ) then
         ev = half * ( u_minus + u_plus ) ! "Shock"
      else
         ev = zero              ! "Expantion fan"
      end if
      
   end function edge_velocity

   
   ! ! Upwind along direction normal to velocity component
   ! function edge_velocity_normal (u_minus, u_plus, u_trasport) result (ev)
      
   !    real(ar), intent(in) :: u_minus, u_plus, u_trasport
   !    real(ar)             :: ev
      
   !    if ( u_transport > zero ) then
   !       ev = u_minus  
   !    else ( u_trasport < zero ) then
   !       ev = u_plus
   !    else
   !       ev = zero
   !    end if
      
   ! end function edge_velocity_normal
   
end module convection_mod
