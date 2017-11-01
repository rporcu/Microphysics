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
   public compute_divuu_x
   public compute_divuu_y
   public compute_divuu_z
   public compute_ugradu_x
   public compute_ugradu_y
   public compute_ugradu_z
   
contains

   ! 
   !
   ! Computes d(uu)/dx + d(uv)/dy + d(uw)/dz at the
   ! x-edges ( u-component location )
   !
   ! 
   subroutine compute_divuu_x ( lo, hi, ug, ulo, uhi, slopes, vg, vlo, vhi, &
        & wg, wlo, whi, divuu_x, dx )  bind(C)

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
      real(ar)                       :: duu, duv, duw
      real(ar)                       :: upls, umns ! Plus, minus
      real(ar)                       :: u_e1, u_w1, u_e2, u_w2
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
               u_e1  = half * ( ug(i+1,j,k) + ug(i,j,k) )
               upls  = ug(i+1,j,k) - half * slopes(i+1,j,k,1)
               umns  = ug(i,j,k)   + half * slopes(i,j,k,1)
               u_e2  = merge ( umns, upls, u_e1 > zero )

               ! West face
               u_w1  = half * ( ug(i,j,k) + ug(i-1,j,k) )
               upls  = ug(i,j,k)   - half * slopes(i,j,k,1)
               umns  = ug(i-1,j,k) + half * slopes(i-1,j,k,1)
               u_w2  = merge ( umns, upls, u_w1 > zero )

               duu   = u_e1*u_e2 - u_w1*u_w2
               
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
                              
               duv   = u_n*v_n - u_s*v_s   
                
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

               duw   = u_t*w_t - u_b*w_b
               
               ! Assemble terms
               divuu_x(i,j,k) = duu*idx + duv*idy + duw*idz
               
            end do
         end do
      end do

   end subroutine compute_divuu_x





   ! 
   !
   ! Computes d(uv)/dx + d(vv)/dy + d(vw)/dz at the
   ! y-edges ( v-component location )
   !
   ! 
   subroutine compute_divuu_y ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, slopes, &
        & wg, wlo, whi, divuu_y, dx )  bind(C)

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
      real(ar)                       :: dvu, dvv, dvw
      real(ar)                       :: vpls, vmns ! Plus, minus
      real(ar)                       :: v_n1, v_s1, v_n2, v_s2
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

               dvu   = u_e*v_e - u_w*v_w
               
               !
               ! d(vv)/dy
               !

               ! North face
               v_n1  = half * ( vg(i,j+1,k) + vg(i,j,k) )
               vpls  = vg(i,j+1,k) - half * slopes(i,j+1,k,2)
               vmns  = vg(i,j,k)   + half * slopes(i,j,k,2)
               v_n2  = merge ( vmns, vpls, v_n1 > zero )   

               ! South face
               v_s1  = half * ( vg(i,j,k)   + vg(i,j-1,k)   )
               vpls  = vg(i,j,k)   - half * slopes(i,j,k,2)
               vmns  = vg(i,j-1,k) + half * slopes(i,j-1,k,2)
               v_s2  = merge ( vmns, vpls, v_s1 > zero )   
                              
               dvv   = v_n1*v_n2 - v_s1*v_s2   
                
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

               dvw   = v_t*w_t - v_b*w_b
               
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
        & wg, wlo, whi, slopes, divuu_z, dx )  bind(C)

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
      real(ar)                       :: dwu, dwv, dww
      real(ar)                       :: wpls, wmns ! Plus, minus
      real(ar)                       :: w_t1, w_b1, w_t2, w_b2
      real(ar)                       :: u_e, u_w, w_e, w_w
      real(ar)                       :: v_n, v_s, w_n, w_s
      
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
               u_e   = half * ( ug(i+1,j,k) + ug(i,j,k-1) )
               wpls  = wg(i+1,j,k) - half * slopes(i+1,j,k,1)
               wmns  = wg(i,j,k)   + half * slopes(i,j,k,1)
               w_e   = merge ( wmns, wpls, u_e > zero )

               ! West face
               u_w   = half * ( ug(i,j,k) + ug(i,j,k-1) )
               wpls  = wg(i,j,k)   - half * slopes(i,j,k,1)
               wmns  = wg(i-1,j,k) + half * slopes(i-1,j,k,1)
               w_w   = merge ( wmns, wpls, u_w > zero )

               dwu   = u_e*w_e - u_w*w_w
               
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
                              
               dwv   = w_n*v_n - w_s*v_s   
                
               !
               ! d(ww)/dz
               !

               ! Top face
               w_t1  = half * ( wg(i,j,k+1) + wg(i,j,k) )
               wpls  = wg(i,j,k+1) - half * slopes(i,j,k+1,3)
               wmns  = wg(i,j,k)   + half * slopes(i,j,k,3)
               w_t2  = merge ( wmns, wpls, w_t1 > zero ) 

               ! Bottom face
               w_b1  = half * ( wg(i,j,k)   + wg(i-1,j,k)   )
               wpls  = wg(i,j,k)   - half * slopes(i,j,k,3)
               wmns  = wg(i,j,k-1) + half * slopes(i,j,k-1,3)
               w_b2   = merge ( wmns, wpls, w_b1 > zero ) 

               dww   = w_t1*w_t2 - w_b1*w_b2
               
               ! Assemble terms
               divuu_z(i,j,k) = dwu*idx + dwv*idy + dww*idz
               
            end do
         end do
      end do

   end subroutine compute_divuu_z

   
   !
   ! Computes udu/dx + vdu/dy + wdu/dz at the
   ! x-edges ( u-component location )
   ! 
   subroutine compute_ugradu_x ( lo, hi, ug, ulo, uhi, vg, vlo, vhi, &
        & wg, wlo, whi, ugradu_x, dx )  bind(C, name="compute_ugradu_x")

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
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(ar),        intent(  out) ::                           &
           & ugradu_x(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
           
      
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

      integer(c_int),  intent(in   ) :: lo(3),  hi(3)   ! Tile indeces for mf associated to vg
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)
      
      real(ar),        intent(in   ) ::                           &
           & ug(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),       &
           & vg(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),       &
           & wg(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(ar),        intent(  out) ::                           &
           & ugradu_y(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      
      
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

      integer(c_int),  intent(in   ) :: lo(3),  hi(3)   ! Tile indeces for mf associated to vg
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      real(ar),        intent(in   ) :: dx(3)
      
      real(ar),        intent(  out) ::                           &
           & ugradu_z(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
          
      real(ar),        intent(in   ) ::                           &
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
