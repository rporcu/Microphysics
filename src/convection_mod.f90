! 
!              
!  This module contains the subroutines to compute the three components
!  of the convetion term (U.grad)U
!  
!  Author: Michele Rosso
! 
!  Date: October 11, 2017
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

      ! Tile bounds
      integer(c_int),  intent(in)    :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in)    :: ulo(3), uhi(3)
      integer(c_int),  intent(in)    :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)

      ! Grid
      real(ar),        intent(in   ) :: dx(3)

      ! Arrays
      real(ar),        intent(inout) ::                           &
           & ugradu_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), & 
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: ududx, vdudy, wdudz
      real(ar)                       :: u, v, w
      real(ar),        parameter     :: over4 = half * half
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! udu/dx
               u     = ug(i,j,k) 
               ududx = ( max ( u, zero ) * ( u - ug(i-1,j,k) )  + &
                     &   min ( u, zero ) * ( ug(i+1,j,k) - u ) ) * idx
               
               ! vdu/dy
               v     = over4 * ( vg(i,j+1,k) + vg(i-1,j+1,k) + &
                     &           vg(i,j,k)   + vg(i-1,j,k)   )

               vdudy = ( max ( v, zero ) * ( ug(i,j,k)   - ug(i,j-1,k) ) + &
                     &   min ( v, zero ) * ( ug(i,j+1,k) - ug(i,j,k)   ) ) * idy  
                

               ! wdu/dz
               w     = over4 * ( wg(i,j,k+1) + wg(i-1,j,k+1) + &
                     &           wg(i,j,k)   + wg(i-1,j,k)   )

               wdudz = ( max ( w, zero ) * ( ug(i,j,k)   - ug(i,j,k-1) ) + &
                     &   min ( w, zero ) * ( ug(i,j,k+1) - ug(i,j,k)   ) ) * idz  

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
      real(ar)                       :: u, v, w
      real(ar)                       :: udvdx, vdvdy, wdvdz
      real(ar),        parameter     :: over4 = half * half
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! udvdx
               u     = over4 * ( ug(i+1,j,k) + ug(i+1,j-1,k) + &
                     &           ug(i,j,k)   + ug(i,j-1,k)   )
              
               udvdx = ( max ( u, zero ) * ( vg(i,j,k) - vg(i-1,j,k) ) + &
                     &   min ( u, zero ) * ( vg(i+1,j,k) - vg(i,j,k) ) ) * idx
               
               ! vdv/dy
               v     = vg(i,j,k)
               vdvdy = ( max ( v, zero ) * ( v - vg(i,j-1,k) ) + &
                     &   min ( v, zero ) * ( vg(i,j+1,k) - v ) ) * idy
               
               ! wdv/dz
               w     = over4 * ( wg(i,j,k+1) + wg(i,j-1,k+1) + &
                     &           wg(i,j,k)   + wg(i,j-1,k)   )
                           
               wdvdz = ( max ( w, zero ) * ( vg(i,j,k) - vg(i,j,k-1) ) + &
                     &   min ( w, zero ) * ( vg(i,j,k+1) - vg(i,j,k) ) ) * idz

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
      real(ar)                       :: u, v, w
      real(ar),        parameter     :: over4 = half * half
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! udwdx
               u     = over4 * ( ug(i+1,j,k) + ug(i+1,j,k-1) + &
                     &           ug(i,j,k)   + ug(i,j,k-1)   )
               
               udwdx = ( max ( u, zero ) * ( wg(i,j,k) - wg(i-1,j,k) ) + &
                     &   min ( u, zero ) * ( wg(i+1,j,k) - wg(i,j,k) ) ) * idx

               ! vdw/dy
               v   = over4 * ( vg(i,j+1,k) + vg(i,j+1,k-1) + &
                               vg(i,j,k)   + vg(i,j,k-1)   )
               
               vdwdy = ( max ( v, zero ) * ( wg(i,j,k) - wg(i,j-1,k) ) + &
                     &   min ( v, zero ) * ( wg(i,j+1,k) - wg(i,j,k) ) ) * idy 
              
               ! wdw/dz
               w     = wg(i,j,k)
               wdwdz = ( max ( w, zero ) * ( w - wg(i,j,k-1) ) + &
                     &   min ( w, zero ) * ( wg(i,j,k+1) - w  ) ) * idz  
               
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
