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

   ! Public members
   public compute_ugradu_x
   public compute_ugradu_y
   public compute_ugradu_z

   ! Small value to protect against tiny velocities used in upwinding
   real(ar),        parameter     :: small_vel = 1.0d-10

contains

   !
   ! Computes udu/dx + vdu/dy + wdu/dz at the
   ! x-edges ( u-component location )
   ! 
   subroutine compute_ugradu_x ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, slopes, ugradu_x, dx )

      ! Tile bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)

      ! Grid
      real(ar),        intent(in   ) :: dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                           &
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)),       &
           & slopes(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3) 

      real(ar),        intent(  out) ::                           &
           & ugradu_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: udu, vdu, wdu
      real(ar)                       :: dpls, dmns  
      real(ar)                       :: u, v, w
      real(ar),        parameter     :: over4 = half * half
      real(ar),        parameter     :: w1 = 1.5_ar
      real(ar),        parameter     :: w2 = 2.0_ar
      real(ar),        parameter     :: w3 = half      
      real(ar)                       :: u_e, u_w, u_s, u_n, u_b, u_t
      real(ar)                       :: v_s, v_n, w_b, w_t

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)


      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! East face
               dpls  = ug(i+1,j,k) - half * slopes(i+1,j,k,1)
               dmns  = ug(i,j,k)   + half * slopes(i,j,k,1)
               u_e   = edge_velocity ( dmns, dpls )

               ! West face
               dpls  = ug(i,j,k)   - half * slopes(i,j,k,1)
               dmns  = ug(i-1,j,k) + half * slopes(i-1,j,k,1)
               u_w   = edge_velocity ( dmns, dpls )

               udu   = ug(i,j,k) * (u_e - u_w)

               ! North face
               v_n   = half * ( vg(i,j+1,k) + vg(i-1,j+1,k) )
               dpls  = ug(i,j+1,k) - half * slopes(i,j+1,k,2)
               dmns  = ug(i,j,k)   + half * slopes(i,j,k,2)
               if (abs(v_n) .gt. small_vel) then
                  u_n   = merge ( dmns, dpls, v_n > zero )   
               else
                  u_n   = half * (dmns + dpls)
               end if

               ! South face
               v_s   = half * ( vg(i,j,k)   + vg(i-1,j,k)   )
               dpls  = ug(i,j,k)   - half * slopes(i,j,k,2)
               dmns  = ug(i,j-1,k) + half * slopes(i,j-1,k,2)
               if (abs(v_s) .gt. small_vel) then
                  u_s   = merge ( dmns, dpls, v_s > zero )   
               else
                  u_s   = half * (dmns + dpls)
               end if

               vdu   = 0.5d0 * (v_n + v_s) * (u_n - u_s)

               ! Top face
               w_t   = half * ( wg(i,j,k+1) + wg(i-1,j,k+1) )
               dpls  = ug(i,j,k+1) - half * slopes(i,j,k+1,3)
               dmns  = ug(i,j,k)   + half * slopes(i,j,k,3)
               if (abs(w_t) .gt. small_vel) then
                  u_t   = merge ( dmns, dpls, w_t > zero ) 
               else
                  u_t   = half * (dmns + dpls)
               end if

               ! Bottom face
               w_b   = half * ( wg(i,j,k)   + wg(i-1,j,k)   )
               dpls  = ug(i,j,k)   - half * slopes(i,j,k,3)
               dmns  = ug(i,j,k-1) + half * slopes(i,j,k-1,3)
               if (abs(w_b) .gt. small_vel) then
                  u_b   = merge ( dmns, dpls, w_b > zero ) 
               else
                  u_b   = half * (dmns + dpls)
               end if

               wdu   = 0.5d0 * (w_t + w_b) * (u_t - u_b)

               ! Assemble terms
               ugradu_x(i,j,k) = udu*idx + vdu*idy + wdu*idz

            end do
         end do
      end do

   end subroutine compute_ugradu_x

   !
   ! Computes udv/dx + vdv/dy + wdv/dz at the
   ! y-edges ( v-component location )
   ! 
   subroutine compute_ugradu_y ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, slopes, ugradu_y, dx )

      integer(c_int),  intent(in   ) :: lo(3),  hi(3)   ! Tile indeces for mf associated to vg
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                           &
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)),       &
           & slopes(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3) 

      real(ar),        intent(  out) ::                           &
           & ugradu_y(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: u, v, w
      real(ar)                       :: dpls, dmns
      real(ar)                       :: udv, vdv, wdv
      real(ar),        parameter     :: over4 = half * half
      real(ar),        parameter     :: w1 = 1.5_ar
      real(ar),        parameter     :: w2 = 2.0_ar
      real(ar),        parameter     :: w3 = half      
      real(ar)                       :: v_e, v_w, v_b, v_t
      real(ar)                       :: u_e, u_w, v_n, v_s, w_b, w_t

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)


      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! East face
               u_e   = half * ( ug(i+1,j,k) + ug(i+1,j-1,k) )
               dpls  = vg(i+1,j,k) - half * slopes(i+1,j,k,1)
               dmns  = vg(i  ,j,k) + half * slopes(i  ,j,k,1)
               if (abs(u_e) .gt. small_vel) then
                  v_e   = merge ( dmns, dpls, u_e > zero )
               else
                  v_e   = half * (dmns + dpls)
               end if

               ! West face
               u_w   = half * ( ug(i,j,k) + ug(i,j-1,k) )
               dpls  = vg(i  ,j,k) - half * slopes(i  ,j,k,1)
               dmns  = vg(i-1,j,k) + half * slopes(i-1,j,k,1)
               v_w   = merge ( dmns, dpls, u_w > zero )
               if (abs(u_w) .gt. small_vel) then
                  v_w   = merge ( dmns, dpls, u_w > zero )
               else
                  v_w   = half * (dmns + dpls)
               end if

               udv   = half * (u_e + u_w) * (v_e - v_w)

               ! North face
               dpls  = vg(i,j+1,k) - half * slopes(i,j+1,k,2)
               dmns  = vg(i,j,k)   + half * slopes(i,j,k,2)
               v_n   = edge_velocity ( dmns, dpls ) 

               ! South face
               dpls  = vg(i,j,k)   - half * slopes(i,j,k,2)
               dmns  = vg(i,j-1,k) + half * slopes(i,j-1,k,2)
               v_s   = edge_velocity ( dmns, dpls )

               vdv   = vg(i,j,k) * (v_n - v_s)

               ! Top face
               w_t   = half * ( wg(i,j,k+1) + wg(i,j-1,k+1) )
               dpls  = vg(i,j,k+1) - half * slopes(i,j,k+1,3)
               dmns  = vg(i,j,k)   + half * slopes(i,j,k,3)
               v_t   = merge ( dmns, dpls, w_t > zero ) 
               if (abs(w_t) .gt. small_vel) then
                  v_t   = merge ( dmns, dpls, w_t > zero ) 
               else
                  v_t   = half * (dmns + dpls)
               end if

               ! Bottom face
               w_b   = half * ( wg(i,j,k)   + wg(i,j-1,k)   )
               dpls  = vg(i,j,k)   - half * slopes(i,j,k,3)
               dmns  = vg(i,j,k-1) + half * slopes(i,j,k-1,3)
               if (abs(w_b) .gt. small_vel) then
                  v_b   = merge ( dmns, dpls, w_b > zero ) 
               else
                  v_b   = half * (dmns + dpls)
               end if

               wdv   = half * (w_t + w_b) * (v_t - v_b)

               ! Assemble terms
               ugradu_y(i,j,k) = udv*idx + vdv*idy + wdv*idz

            end do
         end do
      end do


   end subroutine compute_ugradu_y

   !
   ! Computes udw/dx + vdw/dy + wdw/dz at the
   ! z-edges ( w-component location )
   ! 
   subroutine compute_ugradu_z ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, slopes, ugradu_z, dx )

      integer(c_int),  intent(in   ) :: lo(3),  hi(3)   ! Tile indeces for mf associated to vg
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                           &
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)),       &
           & slopes(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),3) 

      real(ar),        intent(  out) ::                           &
           & ugradu_z(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: udw, vdw, wdw
      real(ar)                       :: dpls, dmns 
      real(ar)                       :: u, v, w
      real(ar),        parameter     :: over4 = half * half
      real(ar),        parameter     :: w1 = 1.5_ar
      real(ar),        parameter     :: w2 = 2.0_ar
      real(ar),        parameter     :: w3 = half      
      real(ar)                       :: u_e, u_w, v_n, v_s
      real(ar)                       :: w_e, w_w, w_n, w_s, w_b, w_t

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! East face
               u_e   = half * ( ug(i+1,j,k) + ug(i+1,j,k-1) )
               dpls  = wg(i+1,j,k) - half * slopes(i+1,j,k,1)
               dmns  = wg(i,j,k)   + half * slopes(i,j,k,1)
               if (abs(u_e) .gt. small_vel) then
                  w_e   = merge ( dmns, dpls, u_e > zero )
               else
                  w_e   = half * (dmns + dpls)
               end if

               ! West face
               u_w   = half * ( ug(i,j,k) + ug(i,j,k-1) )
               dpls  = wg(i,j,k)   - half * slopes(i,j,k,1)
               dmns  = wg(i-1,j,k) + half * slopes(i-1,j,k,1)
               if (abs(u_w) .gt. small_vel) then
                  w_w   = merge ( dmns, dpls, u_w > zero )
               else
                  w_w   = half * (dmns + dpls)
               end if

               udw   = 0.5d0 * (u_e + u_w) * (w_e - w_w)

               ! North face
               v_n   = half * ( vg(i,j+1,k) + vg(i,j+1,k-1) )
               dpls  = wg(i,j+1,k) - half * slopes(i,j+1,k,2)
               dmns  = wg(i,j,k)   + half * slopes(i,j,k,2)
               if (abs(v_n) .gt. small_vel) then
                  w_n   = merge ( dmns, dpls, v_n > zero )   
               else
                  w_n   = half * (dmns + dpls)
               end if

               ! South face
               v_s   = half * ( vg(i,j,k)   + vg(i,j,k-1)   )
               dpls  = wg(i,j,k)   - half * slopes(i,j,k,2)
               dmns  = wg(i,j-1,k) + half * slopes(i,j-1,k,2)
               if (abs(v_s) .gt. small_vel) then
                  w_s   = merge ( dmns, dpls, v_s > zero )   
               else
                  w_s   = half * (dmns + dpls)
               end if

               vdw   = 0.5d0 * (v_n + v_s) * (w_n - w_s)

               ! Top face
               dpls  = wg(i,j,k+1) - half * slopes(i,j,k+1,3)
               dmns  = wg(i,j,k)   + half * slopes(i,j,k,3)
               w_t   = edge_velocity ( dmns, dpls )

               ! Bottom face
               dpls  = wg(i,j,k)   - half * slopes(i,j,k,3)
               dmns  = wg(i,j,k-1) + half * slopes(i,j,k-1,3)
               w_b   = edge_velocity ( dmns, dpls )

               wdw   =  wg(i,j,k) * (w_t - w_b)

               ! Assemble terms
               ugradu_z(i,j,k) = udw*idx + vdw*idy + wdw*idz

            end do
         end do
      end do


   end subroutine compute_ugradu_z


   ! Upwind along direction normal to velocity component
   function edge_velocity ( umns, upls ) result (ev)

      real(ar), intent(in) :: umns, upls
      real(ar)             :: ev, avg

      if ( umns < zero .and. upls > zero ) then
         ev = zero
      else 
         avg = half * ( upls + umns )
         ev = merge ( umns, upls, avg >= zero ) 
      end if

   end function edge_velocity

end module convection_mod
