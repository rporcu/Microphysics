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
   public compute_divuu_x
   public compute_divuu_y
   public compute_divuu_z

   ! Small value to protect against tiny velocities used in upwinding
   real(ar),        parameter     :: small_vel = 1.0d-10
   
contains
   
   !
   ! Computes udu/dx + vdu/dy + wdu/dz at the
   ! x-edges ( u-component location )
   ! 
   subroutine compute_ugradu_x ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, slopes, ugradu_x, dx, order )

      ! Tile bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)

      ! Grid
      real(ar),        intent(in   ) :: dx(3)

      ! Order of the method
      integer(c_int),  intent(in   ) :: order

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
      
      if ( order == 1 ) then 

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  ! udu/dx
                  u     = ug(i,j,k)
                  dmns  = u - ug(i-1,j,k)
                  dpls  = ug(i+1,j,k) - u 
                  udu   = u * merge ( dmns, dpls, u > zero ) 

                  ! vdu/dy
                  v     = over4 * ( vg(i,j+1,k) + vg(i-1,j+1,k) + &
                       &            vg(i,j,k)   + vg(i-1,j,k)   )
                  dmns  = ug(i,j,k)   - ug(i,j-1,k)
                  dpls  = ug(i,j+1,k) - ug(i,j,k)
                  vdu   = v * merge ( dmns, dpls, v > zero ) 

                  ! wdu/dz
                  w     = over4 * ( wg(i,j,k+1) + wg(i-1,j,k+1) + &
                       &           wg(i,j,k)   + wg(i-1,j,k)   )
                  dmns  = ug(i,j,k)   - ug(i,j,k-1)
                  dpls  = ug(i,j,k+1) - ug(i,j,k)
                  wdu   = w * merge ( dmns, dpls, w > zero )

                  ! Assemble terms
                  ugradu_x(i,j,k) = udu*idx + vdu*idy + wdu*idz

               end do
            end do
         end do

      else if ( order == 2 ) then
         
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
 
                  ! ************************************************************

                  ! udu/dx
                  ! u     = ug(i,j,k)
                  ! dmns  =  w1*ug(i,j,k)   - w2*ug(i-1,j,k) + w3*ug(i-2,j,k)
                  ! dpls  = -w3*ug(i+2,j,k) + w2*ug(i+1,j,k) - w1*ug(i,j,k)   
                  ! udu   = u * merge ( dmns, dpls, u > zero ) 

                  ! East face
                  dpls  = ug(i+1,j,k) - half * slopes(i+1,j,k,1)
                  dmns  = ug(i,j,k)   + half * slopes(i,j,k,1)
                  u_e   = edge_velocity ( dmns, dpls )
               
                  ! West face
                  dpls  = ug(i,j,k)   - half * slopes(i,j,k,1)
                  dmns  = ug(i-1,j,k) + half * slopes(i-1,j,k,1)
                  u_w   = edge_velocity ( dmns, dpls )
               
                  udu   = ug(i,j,k) * (u_e - u_w)
 
                  ! ************************************************************

                  ! vdu/dy
                  ! v     = over4 * ( vg(i,j+1,k) + vg(i-1,j+1,k) + & 
                  !                   vg(i,j,k)   + vg(i-1,j,k)   )
                  ! dmns  =  w1*ug(i,j,k)   - w2*ug(i,j-1,k) + w3*ug(i,j-2,k)
                  ! dpls  = -w3*ug(i,j+2,k) + w2*ug(i,j+1,k) - w1*ug(i,j,k)   
                  ! vdu   = v * merge ( dmns, dpls, v > zero ) 

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
 
                  ! ************************************************************

                  ! wdu/dz
                  ! w     = over4 * ( wg(i,j,k+1) + wg(i-1,j,k+1) + &
                  !      &           wg(i,j,k)   + wg(i-1,j,k)   )
                  ! dmns  =  w1*ug(i,j,k)   - w2*ug(i,j,k-1) + w3*ug(i,j,k-2)
                  ! dpls  = -w3*ug(i,j,k+2) + w2*ug(i,j,k+1) - w1*ug(i,j,k)
                  ! wdu   = w * merge ( dmns, dpls, w > zero )

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
 
                  ! ************************************************************
                  ! Assemble terms
                  ugradu_x(i,j,k) = udu*idx + vdu*idy + wdu*idz
                   
               end do
            end do
         end do

      else

         print*, " order can only be 1 or 2 ! " 
         stop
         
      end if

   end subroutine compute_ugradu_x

   !
   ! Computes udv/dx + vdv/dy + wdv/dz at the
   ! y-edges ( v-component location )
   ! 
   subroutine compute_ugradu_y ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, slopes, ugradu_y, dx, order )

      integer(c_int),  intent(in   ) :: lo(3),  hi(3)   ! Tile indeces for mf associated to vg
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)

      ! Order of the method
      integer(c_int),  intent(in   ) :: order
      
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

      if ( order == 1 ) then 

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  ! udvdx
                  u     = over4 * ( ug(i+1,j,k) + ug(i+1,j-1,k) + &
                       &            ug(i,j,k)   + ug(i,j-1,k)   )
                  dmns  = vg(i,j,k) - vg(i-1,j,k)
                  dpls  = vg(i+1,j,k) - vg(i,j,k)
                  udv   = u * merge ( dmns, dpls, u > zero )

                  ! vdv/dy
                  v     = vg(i,j,k)
                  dmns  = vg(i,j,k) - vg(i,j-1,k)
                  dpls  = vg(i,j+1,k) - vg(i,j,k)
                  vdv   = v * merge ( dmns, dpls, v > zero )

                  ! wdv/dz
                  w     = over4 * ( wg(i,j,k+1) + wg(i,j-1,k+1) + &
                       &           wg(i,j,k)   + wg(i,j-1,k)   )
                  dmns  = vg(i,j,k) - vg(i,j,k-1)
                  dpls  = vg(i,j,k+1) - vg(i,j,k)
                  wdv   = w * merge ( dmns, dpls, w > zero )

                  ! Assemble terms
                  ugradu_y(i,j,k) = udv*idx + vdv*idy + wdv*idz

               end do
            end do
         end do
         
      else if ( order == 2 ) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
 
                  ! ************************************************************
                  ! udvdx
                  ! ************************************************************
                  ! u     = over4 * ( ug(i+1,j,k) + ug(i+1,j-1,k) + &
                  !                   ug(i,j,k)   + ug(i,j-1,k)   )
                  ! dmns  =  w1*vg(i,j,k)   - w2*vg(i-1,j,k) + w3*vg(i-2,j,k)
                  ! dpls  = -w3*vg(i+2,j,k) + w2*vg(i+1,j,k) - w1*vg(i,j,k)
                  ! udv   = u * merge ( dmns, dpls, u > zero )

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

                  ! ************************************************************
                  ! vdv/dy
                  ! ************************************************************
                  ! v     = vg(i,j,k)
                  ! dmns  =  w1*vg(i,j,k)   - w2*vg(i,j,k-1) + w3*vg(i,j,k-2)
                  ! dpls  = -w3*vg(i,j+2,k) + w2*vg(i,j,k+1) - w1*vg(i,j,k)
                  ! vdv   = v * merge ( dmns, dpls, v > zero )

                  ! North face
                  dpls  = vg(i,j+1,k) - half * slopes(i,j+1,k,2)
                  dmns  = vg(i,j,k)   + half * slopes(i,j,k,2)
                  v_n   = edge_velocity ( dmns, dpls ) 
   
                  ! South face
                  dpls  = vg(i,j,k)   - half * slopes(i,j,k,2)
                  dmns  = vg(i,j-1,k) + half * slopes(i,j-1,k,2)
                  v_s   = edge_velocity ( dmns, dpls )
                                 
                  vdv   = vg(i,j,k) * (v_n - v_s)
 
                  ! ************************************************************
                  ! wdv/dz
                  ! ************************************************************
                  ! w     = over4 * ( wg(i,j,k+1) + wg(i,j-1,k+1) + &
                  !                   wg(i,j,k)   + wg(i,j-1,k)   )
                  ! dmns  =  w1*vg(i,j,k)   - w2*vg(i,j,k-1) + w3*vg(i,j,k-2)
                  ! dpls  = -w3*vg(i,j,k+2) + w2*vg(i,j,k+1) - w1*vg(i,j,k)
                  ! wdv   = w * merge ( dmns, dpls, w > zero )

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

                  ! ************************************************************

                  ! Assemble terms
                  ugradu_y(i,j,k) = udv*idx + vdv*idy + wdv*idz

               end do
            end do
         end do
         
      else
         print*, " order can only be 1 or 2 ! " 
         stop
      end if

      
   end subroutine compute_ugradu_y

   !
   ! Computes udw/dx + vdw/dy + wdw/dz at the
   ! z-edges ( w-component location )
   ! 
   subroutine compute_ugradu_z ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, slopes, ugradu_z, dx, order )

      integer(c_int),  intent(in   ) :: lo(3),  hi(3)   ! Tile indeces for mf associated to vg
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)

      ! Order of the method
      integer(c_int),  intent(in   ) :: order
      
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

      if ( order == 1 ) then
         
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  ! udwdx
                  u     = over4 * ( ug(i+1,j,k) + ug(i+1,j,k-1) + &
                       &            ug(i,j,k)   + ug(i,j,k-1)   )
                  dmns  = wg(i,j,k) - wg(i-1,j,k)
                  dpls  = wg(i+1,j,k) - wg(i,j,k)
                  udw   = u * merge ( dmns, dpls, u > zero )

                  ! vdw/dy
                  v   = over4 * ( vg(i,j+1,k) + vg(i,j+1,k-1) + &
                       &          vg(i,j,k)   + vg(i,j,k-1)   )
                  dmns  = wg(i,j,k) - wg(i,j-1,k)
                  dpls  = wg(i,j+1,k) - wg(i,j,k)
                  vdw   = v * merge ( dmns, dpls, v > zero )

                  ! wdw/dz
                  w     = wg(i,j,k)
                  dmns  = wg(i,j,k) - wg(i,j,k-1)
                  dpls  = wg(i,j,k+1) - wg(i,j,k)
                  wdw   = w * merge ( dmns, dpls, w > zero )

                  ! Assemble terms
                  ugradu_z(i,j,k) = udw*idx + vdw*idy + wdw*idz


               end do
            end do
         end do

      else if ( order == 2 ) then 

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
 
                  ! ************************************************************

                  ! udwdx
                  ! u     = over4 * ( ug(i+1,j,k) + ug(i+1,j,k-1) + &
                  !      &            ug(i,j,k)   + ug(i,j,k-1)   )
                  ! dmns  =  w1*wg(i,j,k)   - w2*wg(i-1,j,k) + w3*wg(i-2,j,k)
                  ! dpls  = -w3*wg(i+2,j,k) + w2*wg(i+1,j,k) - w2*wg(i,j,k)
                  ! udw   = u * merge ( dmns, dpls, u > zero )

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
 
                  ! ************************************************************

                  ! vdw/dy
                  ! v   = over4 * ( vg(i,j+1,k) + vg(i,j+1,k-1) + &
                  !      &          vg(i,j,k)   + vg(i,j,k-1)   )
                  ! dmns  =  w1*wg(i,j,k)   - w2*wg(i,j-1,k) + w3*wg(i,j-2,k)
                  ! dpls  = -w3*wg(i,j+2,k) + w2*wg(i,j+1,k) - w2*wg(i,j,k)
                  ! vdw   = v * merge ( dmns, dpls, v > zero )

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
 
                  ! ************************************************************

                  ! wdw/dz
                  ! w     = wg(i,j,k)
                  ! dmns  =  w1*wg(i,j,k)   - w2*wg(i,j,k-1) + w3*wg(i,j,k-2)
                  ! dpls  = -w3*wg(i,j,k+2) + w2*wg(i,j,k+1) - w2*wg(i,j,k)
                  ! wdw   = w * merge ( dmns, dpls, w > zero )

                  ! Top face
                  dpls  = wg(i,j,k+1) - half * slopes(i,j,k+1,3)
                  dmns  = wg(i,j,k)   + half * slopes(i,j,k,3)
                  w_t   = edge_velocity ( dmns, dpls )
   
                  ! Bottom face
                  dpls  = wg(i,j,k)   - half * slopes(i,j,k,3)
                  dmns  = wg(i,j,k-1) + half * slopes(i,j,k-1,3)
                  w_b   = edge_velocity ( dmns, dpls )
   
                  wdw   =  wg(i,j,k) * (w_t - w_b)
 
                  ! ************************************************************

                  ! Assemble terms
                  ugradu_z(i,j,k) = udw*idx + vdw*idy + wdw*idz

               end do
            end do
         end do
      else
         print*, " order can only be 1 or 2 ! " 
         stop
      end if

      
   end subroutine compute_ugradu_z


   ! 
   !
   ! Computes d(uu)/dx + d(uv)/dy + d(uw)/dz at the
   ! x-edges ( u-component location )
   !
   ! 
   subroutine compute_divuu_x ( lo, hi, ug, ulo, uhi, slopes, vg, vlo, vhi, &
        & wg, wlo, whi, divuu_x, dx )

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
           & divuu_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
          
      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: dmns, dpls, u, v, w
      real(ar)                       :: duu, duv, duw
      real(ar)                       :: upls, umns ! Plus, minus
      real(ar)                       :: u_e1, u_w1, u_e2, u_w2, u_e, u_w
      real(ar)                       :: u_n, u_s, v_n, v_s
      real(ar)                       :: u_t, u_b, w_t, w_b
      real(ar),        parameter     :: over4 = half * half
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! 
               ! d(uu)/dx
               !

               ! East face
               upls  = ug(i+1,j,k) - half * slopes(i+1,j,k,1)
               umns  = ug(i,j,k)   + half * slopes(i,j,k,1)
               u_e   = edge_velocity ( umns, upls )
               
               ! West face
               upls  = ug(i,j,k)   - half * slopes(i,j,k,1)
               umns  = ug(i-1,j,k) + half * slopes(i-1,j,k,1)
               u_w   = edge_velocity ( umns, upls )
               
               ! duu   = u_e*u_e - u_w*u_w
               ! duu   = 0.5d0 * (u_e + u_w) * (u_e - u_w) ! CONV 2
               duu   = ug(i,j,k) * (u_e - u_w)

               ! ***************************************************************
               ! CONV udu/dx
               ! u     = ug(i,j,k)
               ! dmns  = u - ug(i-1,j,k)
               ! dpls  = ug(i+1,j,k) - u
               ! duu   = u * merge ( dmns, dpls, u > zero )
               ! ***************************************************************
                              
               !
               ! d(uv)/dy
               !

               ! North face
               v_n   = half * ( vg(i,j+1,k) + vg(i-1,j+1,k) )
               upls  = ug(i,j+1,k) - half * slopes(i,j+1,k,2)
               umns  = ug(i,j,k)   + half * slopes(i,j,k,2)
               u_n   = merge ( umns, upls, v_n > zero )   

               ! South face
               v_s   = half * ( vg(i,j,k)   + vg(i-1,j,k)   )
               upls  = ug(i,j,k)   - half * slopes(i,j,k,2)
               umns  = ug(i,j-1,k) + half * slopes(i,j-1,k,2)
               u_s   = merge ( umns, upls, v_s > zero )   
                              
               ! duv   = u_n*v_n - u_s*v_s   
               duv   = 0.5d0 * (v_n + v_s) * (u_n - u_s)

               ! ***************************************************************
               ! CONV vdu/dy
               ! v     = over4 * ( vg(i,j+1,k) + vg(i-1,j+1,k) + &
               !      &            vg(i,j,k)   + vg(i-1,j,k)   )
               ! dmns  = ug(i,j,k)   - ug(i,j-1,k)
               ! dpls  = ug(i,j+1,k) - ug(i,j,k)
               ! duv   = v * merge ( dmns, dpls, v > zero )
               ! ***************************************************************
                
               !
               ! d(uw)/dz
               !

               ! Top face
               w_t   = half * ( wg(i,j,k+1) + wg(i-1,j,k+1) )
               upls  = ug(i,j,k+1) - half * slopes(i,j,k+1,3)
               umns  = ug(i,j,k)   + half * slopes(i,j,k,3)
               u_t   = merge ( umns, upls, w_t > zero ) 

               ! Bottom face
               w_b   = half * ( wg(i,j,k)   + wg(i-1,j,k)   )
               upls  = ug(i,j,k)   - half * slopes(i,j,k,3)
               umns  = ug(i,j,k-1) + half * slopes(i,j,k-1,3)
               u_b   = merge ( umns, upls, w_b > zero ) 

               ! duw   = u_t*w_t - u_b*w_b
               duw   = 0.5d0 * (w_t + w_b) * (u_t - u_b)

               ! ***************************************************************
               ! CONV wdu/dz
               ! w     = over4 * ( wg(i,j,k+1) + wg(i-1,j,k+1) + &
               !                   wg(i,j,k)   + wg(i-1,j,k)   )
               ! dmns  = ug(i,j,k)   - ug(i,j,k-1)
               ! dpls  = ug(i,j,k+1) - ug(i,j,k)
               ! duw   = w * merge ( dmns, dpls, w > zero )
               ! ***************************************************************
               
               ! Assemble terms
               divuu_x(i,j,k) = duu*idx + duv*idy + duw*idz
               
            end do
         end do
      end do

   end subroutine compute_divuu_x

   ! 
   !
   ! Computes d(vu)/dx + d(vv)/dy + d(vw)/dz at the
   ! y-edges ( v-component location )
   !
   ! 
   subroutine compute_divuu_y ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, slopes, &
        & wg, wlo, whi, divuu_y, dx )

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
           & slopes(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3) 

      real(ar),        intent(  out) ::                           &
           & divuu_y(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
          
      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: dmns, dpls, u, v, w
      real(ar)                       :: dvu, dvv, dvw
      real(ar)                       :: vpls, vmns ! Plus, minus
      real(ar)                       :: v_n1, v_s1, v_n2, v_s2, v_n, v_s
      real(ar)                       :: u_e, u_w, v_e, v_w
      real(ar)                       :: v_t, v_b, w_t, w_b
      real(ar),        parameter     :: over4 = half * half
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! 
               ! d(vu)/dx
               !

               ! East face
               u_e   = half * ( ug(i+1,j,k) + ug(i+1,j-1,k) )
               vpls  = vg(i+1,j,k) - half * slopes(i+1,j,k,1)
               vmns  = vg(i,j,k)   + half * slopes(i,j,k,1)
               v_e   = merge ( vmns, vpls, u_e > zero )

               ! West face
               u_w   = half * ( ug(i,j,k) + ug(i,j-1,k) )
               vpls  = vg(i,j,k)   - half * slopes(i,j,k,1)
               vmns  = vg(i-1,j,k) + half * slopes(i-1,j,k,1)
               v_w   = merge ( vmns, vpls, u_w > zero )

               ! dvu   = u_e*v_e - u_w*v_w
               dvu   = 0.5d0 * (u_e + u_w) * (v_e - v_w)

               ! ***************************************************************
               ! CONV udv/dx
               ! u     = over4 * ( ug(i+1,j,k) + ug(i+1,j-1,k) + &
               !                   ug(i,j,k)   + ug(i,j-1,k)   )
               ! dmns  = vg(i,j,k) - vg(i-1,j,k)
               ! dpls  = vg(i+1,j,k) - vg(i,j,k)
               ! dvu   = u * merge ( dmns, dpls, u > zero )
               ! ***************************************************************
               
               !
               ! d(vv)/dy
               !

               ! North face
               vpls  = vg(i,j+1,k) - half * slopes(i,j+1,k,2)
               vmns  = vg(i,j,k)   + half * slopes(i,j,k,2)
               v_n   = edge_velocity ( vmns, vpls ) 

               ! South face
               vpls  = vg(i,j,k)   - half * slopes(i,j,k,2)
               vmns  = vg(i,j-1,k) + half * slopes(i,j-1,k,2)
               v_s   = edge_velocity ( vmns, vpls )
                              
               ! dvv   = v_n*v_n - v_s*v_s   
               ! dvv   = 0.5d0 * (v_n + v_s) * (v_n - v_s) ! CONV 2
               dvv   = vg(i,j,k) * (v_n - v_s)

               ! ***************************************************************
               ! CONV vdv/dy
               ! dmns  = vg(i,j  ,k) - vg(i,j-1,k)
               ! dpls  = vg(i,j+1,k) - vg(i,j  ,k) 
               ! dvv   = vg(i,j  ,k) * merge ( dmns, dpls, vg(i,j,k) > zero ) 
               ! ***************************************************************
                
               !
               ! d(uw)/dz
               !

               ! Top face
               w_t   = half * ( wg(i,j,k+1) + wg(i,j-1,k+1) )
               vpls  = vg(i,j,k+1) - half * slopes(i,j,k+1,3)
               vmns  = vg(i,j,k)   + half * slopes(i,j,k,3)
               v_t   = merge ( vmns, vpls, w_t > zero ) 

               ! Bottom face
               w_b   = half * ( wg(i,j,k)   + wg(i,j-1,k)   )
               vpls  = vg(i,j,k)   - half * slopes(i,j,k,3)
               vmns  = vg(i,j,k-1) + half * slopes(i,j,k-1,3)
               v_b   = merge ( vmns, vpls, w_b > zero ) 

               ! dvw   = v_t*w_t - v_b*w_b
               dvw   = 0.5d0 * (w_t + w_b) * (v_t - v_b)

               ! ***************************************************************
               ! CONV wdv/dz
               ! w     = over4 * ( wg(i,j,k+1) + wg(i,j-1,k+1) + &
               !                   wg(i,j,k)   + wg(i,j-1,k)   )
               ! dmns  = vg(i,j,k) - vg(i,j,k-1)
               ! dpls  = vg(i,j,k+1) - vg(i,j,k)
               ! dvw   = w * merge ( dmns, dpls, w > zero )
               ! ***************************************************************
               
               ! Assemble terms
               divuu_y(i,j,k) = dvu*idx + dvv*idy + dvw*idz
               
            end do
         end do
      end do

   end subroutine compute_divuu_y


   ! 
   !
   ! Computes d(wu)/dx + d(wv)/dy + d(ww)/dz at the
   ! z-edges ( w-component location )
   !
   ! 
   subroutine compute_divuu_z ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, slopes, divuu_z, dx )

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
           & slopes(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),3) 

      real(ar),        intent(  out) ::                           &
           & divuu_z(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
          
      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: dmns, dpls, u, v, w
      real(ar)                       :: dwu, dwv, dww
      real(ar)                       :: wpls, wmns ! Plus, minus
      real(ar)                       :: w_t1, w_b1, w_t2, w_b2, w_t, w_b
      real(ar)                       :: u_e, u_w, w_e, w_w
      real(ar)                       :: v_n, v_s, w_n, w_s
      real(ar),        parameter     :: over4 = half * half
      
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! 
               ! d(wu)/dx
               !

               ! East face
               u_e   = half * ( ug(i+1,j,k) + ug(i+1,j,k-1) )
               wpls  = wg(i+1,j,k) - half * slopes(i+1,j,k,1)
               wmns  = wg(i,j,k)   + half * slopes(i,j,k,1)
               w_e   = merge ( wmns, wpls, u_e > zero )

               ! West face
               u_w   = half * ( ug(i,j,k) + ug(i,j,k-1) )
               wpls  = wg(i,j,k)   - half * slopes(i,j,k,1)
               wmns  = wg(i-1,j,k) + half * slopes(i-1,j,k,1)
               w_w   = merge ( wmns, wpls, u_w > zero )

               ! dwu   = u_e*w_e - u_w*w_w
               ! dwu   = 0.5d0 * (u_e + u_w) * (w_e - w_w)

               ! ***************************************************************
               ! CONV udwdx
               u     = over4 * ( ug(i+1,j,k) + ug(i+1,j,k-1) + &
                                 ug(i,j,k)   + ug(i,j,k-1)   )
               dmns  = wg(i,j,k) - wg(i-1,j,k)
               dpls  = wg(i+1,j,k) - wg(i,j,k)
               dwu   = u * merge ( dmns, dpls, u > zero )
               ! ***************************************************************
               
               !
               ! d(wv)/dy
               !

               ! North face
               v_n   = half * ( vg(i,j+1,k) + vg(i,j+1,k-1) )
               wpls  = wg(i,j+1,k) - half * slopes(i,j+1,k,2)
               wmns  = wg(i,j,k)   + half * slopes(i,j,k,2)
               w_n   = merge ( wmns, wpls, v_n > zero )   

               ! South face
               v_s   = half * ( vg(i,j,k)   + vg(i,j,k-1)   )
               wpls  = wg(i,j,k)   - half * slopes(i,j,k,2)
               wmns  = wg(i,j-1,k) + half * slopes(i,j-1,k,2)
               w_s   = merge ( wmns, wpls, v_s > zero )   
                              
               ! dwv   = w_n*v_n - w_s*v_s   
               ! dwv   = 0.5d0 * (v_n + v_s) * (w_n - w_s)

               ! ***************************************************************
               ! CONV vdw/dy
               v   = over4 * ( vg(i,j+1,k) + vg(i,j+1,k-1) + &
                    &          vg(i,j,k)   + vg(i,j,k-1)   )
               dmns  = wg(i,j,k) - wg(i,j-1,k)
               dpls  = wg(i,j+1,k) - wg(i,j,k)
               dwv   = v * merge ( dmns, dpls, v > zero )
               ! ***************************************************************
                
               !
               ! d(ww)/dz
               !

               ! Top face
               wpls  = wg(i,j,k+1) - half * slopes(i,j,k+1,3)
               wmns  = wg(i,j,k)   + half * slopes(i,j,k,3)
               w_t   = edge_velocity ( wmns, wpls )

               ! Bottom face
               wpls  = wg(i,j,k)   - half * slopes(i,j,k,3)
               wmns  = wg(i,j,k-1) + half * slopes(i,j,k-1,3)
               w_b   = edge_velocity ( wmns, wpls )

               ! dww   = w_t*w_t - w_b*w_b
               ! dww   = 0.5d0 * (w_t + w_b) * (w_t - w_b) ! CONV 2
               ! dww   =  wg(i,j,k) * (w_t - w_b)

               dmns  = wg(i,j,k  ) - wg(i,j,k-1)
               dpls  = wg(i,j,k+1) - wg(i,j,k  ) 
               ! dww   = wg(i,j,k  ) * merge ( dmns, dpls, wg(i,j,k) > zero ) 

               ! ***************************************************************
               ! CONV wdw/dz
               w     = wg(i,j,k)
               dmns  = wg(i,j,k) - wg(i,j,k-1)
               dpls  = wg(i,j,k+1) - wg(i,j,k)
               dww   = w * merge ( dmns, dpls, w > zero )
               ! ***************************************************************
               
               ! Assemble terms
               divuu_z(i,j,k) = dwu*idx + dwv*idy + dww*idz
               
            end do
         end do
      end do

   end subroutine compute_divuu_z

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

   
   
   ! ! Upwind along the same direction as the velocity component
   ! function edge_velocity (u_minus, u_plus) result (ev)
      
   !    real(ar), intent(in) :: u_minus, u_plus
   !    real(ar)             :: ev
      
   !    if ( ( u_minus >= zero ) .and. ( u_plus >= zero ) ) then
   !       ev = u_minus
   !    else if ( ( u_minus <= zero ) .and. ( u_plus <= zero ) ) then
   !       ev = u_plus
   !    else if ( ( u_minus > zero ) .and. ( u_plus < zero ) ) then
   !       ev = half * ( u_minus + u_plus ) ! "Shock"
   !    else
   !       ev = zero              ! "Expantion fan"
   !    end if
      
   ! end function edge_velocity




   
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
