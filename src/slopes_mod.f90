!   
!  This module contains the subroutines to compute the slopes in the
!  three directions for each velocity component.
!  The x,y, and z slopes for each velocity component are calculated at the
!  velocity component location via the second order Monotonized Central (MC)
!  limiter (van Leer, 1977). The scheme is described below for the u-velocity.
!
!
!                      |--x--|--x--|--x--|
!                        i-1    i    i+1
!
!  In the sketch above, the x represents the u-velocities while the vertical 
!  bars | | enclose a u-cell (NOT a scalar cell!).
!  The MC limiter computes the slope at cell "i" by combining the left, central
!  and right u-variation "du":
!
!       du_l = u(i) - u(i-1)               = left variation
!       du_c = 0.5 * ( u(i+1) - u(i-1) )   = central (umlimited) variation
!       du_r = u(i+1) - u(i)               = right variation
!
!  Finally, the u-variation at cell "i" is given by :
!
!       du(i) = sign(du_c) min(2|du_l|, |du_c|, 2|du_r|)) if du_l*du_r > 0
!       du(i) = 0                                         otherwise 
!
!  The above procedure is applied direction by direction.
!
!  BOUNDARY CONDITIONS
!  When periodic or Neumann's BCs are imposed, the scheme can be applied
!  without any change since the ghost cells at the boundary are filled
!  by either periodicity or by extrapolation.
!  For Dirichlet's BCs in the transversal direction, the scheme can again
!  be applied as is since the velocity is known at the first ghost cell
!  out of the domain.
!  However, for Dirichlet's BCs in the longitudinal direction, the velocity
!  is not known outside the domain since the BC is applied directly at the first
!  valid node which lies on the boundary itself. Therefore, the scheme must be
!  arranged as follows to use ONLY values from inside the domain.
!  For a left boundary (i=0), the u-variations are:
!
!       du_l = 0                             Don't use values on the left
!       du_c = -1.5*u(0) + 2*u(1) -0.5*u(2)  2nd order right-biased  
!       du_r = u(1) - u(0)                   Right variation
!

!  
!  Author: Michele Rosso
! 
!  Date: October 30, 2017
!
! 
module slopes_mod
   
   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   use bc,                only: minf_, nsw_, fsw_
   
   implicit none
   private

   public compute_u_slopes
   public compute_v_slopes
   public compute_w_slopes
   

contains


   !
   ! Compute u-velocity slopes
   ! 
   subroutine compute_u_slopes ( lo, hi, u, ulo, uhi, slopes, &
        & domlo, domhi, bc_ilo_type, bc_ihi_type ) bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)

      ! Grid bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! BC types 
      integer(c_int), intent(in   ) ::  &
           & bc_ilo_type(domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2), &
           & bc_ihi_type(domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      
      ! Arrays
      real(ar),       intent(in   ) ::                      &
           &   u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      
      real(ar),       intent(  out) ::                           &
           & slopes(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer                       :: i, j, k, slope_order
      real(ar)                      :: du_l, du_c, du_r, ds
      real(ar),    parameter        :: two = 2.0_ar, three2nds = 1.5_ar  

      integer                       :: cen, lim, flg, frm

      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flg = 3 )
      parameter( frm = 4 )

      real(ar), allocatable :: scr(:,:)

      slope_order = 2

      if (slope_order .eq. 0) then

         slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0

      else if (slope_order .eq. 2) then
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! X direction
               du_l = u(i,j,k) - u(i-1,j,k)
               du_c = half * ( u(i+1,j,k) - u(i-1,j,k) )
               du_r = u(i+1,j,k) - u(i,j,k)
               
               slopes(i,j,k,1) = mc_limiter ( du_l, du_c, du_r )  

               ! Y direction
               du_l = u(i,j,k) - u(i,j-1,k)
               du_c = half * ( u(i,j+1,k) - u(i,j-1,k) )
               du_r = u(i,j+1,k) - u(i,j,k)

               slopes(i,j,k,2) = mc_limiter ( du_l, du_c, du_r )  

               ! z direction
               du_l = u(i,j,k) - u(i,j,k-1)
               du_c = half * ( u(i,j,k+1) - u(i,j,k-1) )
               du_r = u(i,j,k+1) - u(i,j,k)

               slopes(i,j,k,3) = mc_limiter ( du_l, du_c, du_r )  
               
            end do
         end do
      end do

      else

      allocate (scr(lo(1)-1:hi(1)+1,4))
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)

            do i = lo(1)-1, hi(1)+1
               ! X direction
               du_l = two * (u(i,j,k) - u(i-1,j,k))
               du_c = half * ( u(i+1,j,k) - u(i-1,j,k) )
               du_r = two * (u(i+1,j,k) - u(i,j,k))

               scr(i,cen) = du_c
               scr(i,lim) = min(abs(du_l),abs(du_r))
               scr(i,lim) = merge(scr(i,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(i,flg) = sign(1.d0,scr(i,cen))
               scr(i,frm) = scr(i,flg) * min(abs(du_c),scr(i,lim))
            end do

            do i = lo(1), hi(1)
               ds = (4.d0/3.d0) * scr(i,cen) - &
                    (1.d0/6.d0) * (scr(i+1,frm) + scr(i-1,frm))
               slopes(i,j,k,1) = scr(i,flg)*min(abs(ds),scr(i,lim))
            end do

         end do
      end do

      deallocate (scr)
      allocate (scr(lo(2)-1:hi(2)+1,4))
      
      do k = lo(3), hi(3)
         do i = lo(1), hi(1)

            do j = lo(2)-1, hi(2)+1
               du_l = two * (u(i,j,k) - u(i,j-1,k))
               du_c = half * ( u(i,j+1,k) - u(i,j-1,k) )
               du_r = two * (u(i+1,j,k) - u(i,j,k))

               scr(j,cen) = du_c
               scr(j,lim) = min(abs(du_l),abs(du_r))
               scr(j,lim) = merge(scr(j,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(j,flg) = sign(1.d0,scr(j,cen))
               scr(j,frm) = scr(j,flg) * min(abs(du_c),scr(j,lim))
            end do

            do j = lo(2), hi(2)
               ds = (4.d0/3.d0) * scr(j,cen) - &
                    (1.d0/6.d0) * (scr(j+1,frm) + scr(j-1,frm))
               slopes(i,j,k,2) = scr(j,flg)*min(abs(ds),scr(j,lim))
            end do

         end do
      end do

      deallocate (scr)
      allocate (scr(lo(3)-1:hi(3)+1,4))
      
      do i = lo(1), hi(1)
         do j = lo(2), hi(2)

            do k = lo(3)-1, hi(3)+1
               du_l = two  * (u(i,j,k) - u(i,j,k-1) )
               du_c = half * ( u(i,j,k+1) - u(i,j,k-1) )
               du_r = two  * (u(i,j,k+1) - u(i,j,k) )

               scr(k,cen) = du_c
               scr(k,lim) = min(abs(du_l),abs(du_r))
               scr(k,lim) = merge(scr(k,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(k,flg) = sign(1.d0,scr(k,cen))
               scr(k,frm) = scr(k,flg) * min(abs(du_c),scr(k,lim))
            end do

            do k = lo(3), hi(3)
               ds = (4.d0/3.d0) * scr(k,cen) - &
                    (1.d0/6.d0) * (scr(k+1,frm) + scr(k-1,frm))
               slopes(i,j,k,3) = scr(k,flg)*min(abs(ds),scr(k,lim))
            end do

         end do
      end do

      deallocate (scr)
      end if

      
      ! ! 
      ! ! Compute slopes at boundary where physical BCs are imposed
      ! ! 
      ! if ( lo(1) == domlo(1) ) then

      !    i = lo(1)
         
      !    do k = lo(3), hi(3)
      !       do j = lo(2), hi(2)

               
      !          if ( ( bc_ilo_type(j,k,1) == MINF_ ) .or. &
      !               ( bc_ilo_type(j,k,1) == NSW_ )  .or. &
      !               ( bc_ilo_type(j,k,1) == FSW_ )  ) then

      !             du_l = zero
      !             du_c = - half*u(i+2,j,k) + two*u(i+1,j,k) - three2nds*u(i,j,k) 
      !             du_r = u(i+1,j,k) - u(i,j,k)
               
      !             slopes(i,j,k,1) = mc_limiter ( du_l, du_c, du_r )  
                  
      !          end if
      !       end do
      !    end do
      ! end if

      ! if ( hi(1) == (domhi(1)+1) ) then

      !    i = hi(1) 
         
      !    do k = lo(3), hi(3)
      !       do j = lo(2), hi(2)

               
      !          if ( ( bc_ilo_type(j,k,1) == MINF_ ) .or. &
      !               ( bc_ilo_type(j,k,1) == NSW_ )  .or. &
      !               ( bc_ilo_type(j,k,1) == FSW_ )  ) then

      !             du_l = u(i,j,k) - u(i-1,j,k)
      !             du_c = half*u(i-2,j,k) - two*u(i-1,j,k) + three2nds*u(i,j,k) 
      !             du_r = zero
               
      !             slopes(i,j,k,1) = mc_limiter ( du_l, du_c, du_r )  
                  
      !          end if
      !       end do
      !    end do
      ! end if

      ! block
      !    real(ar)   :: usl_x, usl_y, usl_z
      !    usl_x = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) )
      !    usl_y = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) )
      !    usl_z = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) )
      !    print*, "Max u slopes = ", usl_x, usl_y, usl_z
      ! end block
      
   end subroutine compute_u_slopes

   !
   ! Compute v-velocity slopes
   ! 
   subroutine compute_v_slopes ( lo, hi, v, vlo, vhi, slopes, &
        & domlo, domhi, bc_jlo_type, bc_jhi_type ) bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)

      ! Grid bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! BC types 
      integer(c_int), intent(in   ) ::  &
           & bc_jlo_type(domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2), &
           & bc_jhi_type(domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      
      ! Arrays
      real(ar),       intent(in   ) ::                      &
           &   v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      
      real(ar),       intent(  out) ::                           &
           & slopes(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

      ! Local variables
      integer                       :: i, j, k, slope_order
      real(ar)                      :: du_l, du_c, du_r, ds
      real(ar),    parameter        :: two = 2.0_ar, three2nds = 1.5_ar  
      
      integer                       :: cen, lim, flg, frm

      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flg = 3 )
      parameter( frm = 4 )

      real(ar), allocatable :: scr(:,:)

      slope_order = 2

      if (slope_order .eq. 0) then

         slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0

      else if (slope_order .eq. 2) then
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! X direction
               du_l = v(i,j,k) - v(i-1,j,k)
               du_c = half * ( v(i+1,j,k) - v(i-1,j,k) )
               du_r = v(i+1,j,k) - v(i,j,k)
               
               slopes(i,j,k,1) = mc_limiter ( du_l, du_c, du_r )  

               ! Y direction
               du_l = v(i,j,k) - v(i,j-1,k)
               du_c = half * ( v(i,j+1,k) - v(i,j-1,k) )
               du_r = v(i,j+1,k) - v(i,j,k)

               slopes(i,j,k,2) = mc_limiter ( du_l, du_c, du_r )  

               ! z direction
               du_l = v(i,j,k) - v(i,j,k-1)
               du_c = half * ( v(i,j,k+1) - v(i,j,k-1) )
               du_r = v(i,j,k+1) - v(i,j,k)

               slopes(i,j,k,3) = mc_limiter ( du_l, du_c, du_r )  
               
            end do
         end do
      end do

      else 

      allocate (scr(lo(1)-1:hi(1)+1,4))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)

            do i = lo(1)-1, hi(1)+1
               ! X direction
               du_l = v(i,j,k) - v(i-1,j,k)
               du_c = half * ( v(i+1,j,k) - v(i-1,j,k) )
               du_r = v(i+1,j,k) - v(i,j,k)

               scr(i,cen) = du_c
               scr(i,lim) = 2.0d0*min(abs(du_l),abs(du_r))
               scr(i,lim) = merge(scr(i,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(i,flg) = sign(1.d0,scr(i,cen))
               scr(i,frm) = scr(i,flg) * min(abs(du_c),scr(i,lim))
            end do

            do i = lo(1), hi(1)
               ds = (4.d0/3.d0) * scr(i,cen) - &
                    (1.d0/6.d0) * (scr(i+1,frm) + scr(i-1,frm))
               slopes(i,j,k,1) = scr(i,flg)*min(abs(ds),scr(i,lim))
            end do

         end do
      end do

      deallocate (scr)
      allocate (scr(lo(2)-1:hi(2)+1,4))
      
      do k = lo(3), hi(3)
         do i = lo(1), hi(1)

            do j = lo(2)-1, hi(2)+1
               du_l = two * (v(i,j,k) - v(i,j-1,k))
               du_c = half * ( v(i,j+1,k) - v(i,j-1,k) )
               du_r = two * (v(i+1,j,k) - v(i,j,k))

               scr(j,cen) = du_c
               scr(j,lim) = min(abs(du_l),abs(du_r))
               scr(j,lim) = merge(scr(j,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(j,flg) = sign(1.d0,scr(j,cen))
               scr(j,frm) = scr(j,flg) * min(abs(du_c),scr(j,lim))
            end do

            do j = lo(2), hi(2)
               ds = (4.d0/3.d0) * scr(j,cen) - &
                    (1.d0/6.d0) * (scr(j+1,frm) + scr(j-1,frm))
               slopes(i,j,k,2) = scr(j,flg)*min(abs(ds),scr(j,lim))
            end do

         end do
      end do

      deallocate (scr)
      allocate (scr(lo(3)-1:hi(3)+1,4))
      
      do i = lo(1), hi(1)
         do j = lo(2), hi(2)

            do k = lo(3)-1, hi(3)+1
               du_l = two  * (v(i,j,k) - v(i,j,k-1) )
               du_c = half * ( v(i,j,k+1) - v(i,j,k-1) )
               du_r = two  * (v(i,j,k+1) - v(i,j,k) )

               scr(k,cen) = du_c
               scr(k,lim) = min(abs(du_l),abs(du_r))
               scr(k,lim) = merge(scr(k,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(k,flg) = sign(1.d0,scr(k,cen))
               scr(k,frm) = scr(k,flg) * min(abs(du_c),scr(k,lim))
            end do

            do k = lo(3), hi(3)
               ds = (4.d0/3.d0) * scr(k,cen) - &
                    (1.d0/6.d0) * (scr(k+1,frm) + scr(k-1,frm))
               slopes(i,j,k,3) = scr(k,flg)*min(abs(ds),scr(k,lim))
            end do

         end do
      end do

      deallocate (scr)
      end if
      
      ! ! 
      ! ! Compute slopes at boundary where physical BCs are imposed
      ! ! 
      ! if ( lo(2) == domlo(2) ) then

      !    j = lo(2)
         
      !    do k = lo(3), hi(3)
      !       do i = lo(1), hi(1)
           
      !          if ( ( bc_jlo_type(i,k,1) == MINF_ ) .or. &
      !               ( bc_jlo_type(i,k,1) == NSW_ )  .or. &
      !               ( bc_jlo_type(i,k,1) == FSW_ )  ) then

      !             du_l = zero
      !             du_c = - half*v(i,j+2,k) + two*v(i,j+1,k) - three2nds*v(i,j,k) 
      !             du_r = v(i,j+1,k) - v(i,j,k)
               
      !             slopes(i,j,k,2) = mc_limiter ( du_l, du_c, du_r )  
                  
      !          end if
      !       end do
      !    end do
         
      ! end if

      ! if ( hi(2) == (domhi(2)+1) ) then

      !    j = hi(2) 
         
      !    do k = lo(3), hi(3)
      !       do i = lo(1), hi(1)

      !          if ( ( bc_jlo_type(i,k,1) == MINF_ ) .or. &
      !               ( bc_jlo_type(i,k,1) == NSW_ )  .or. &
      !               ( bc_jlo_type(i,k,1) == FSW_ )  ) then

      !             du_l = v(i,j,k) - v(i,j-1,k)
      !             du_c = half*v(i,j-2,k) - two*v(i,j-1,k) + three2nds*v(i,j,k) 
      !             du_r = zero
               
      !             slopes(i,j,k,2) = mc_limiter ( du_l, du_c, du_r )  
                  
      !          end if
      !       end do
      !    end do
         
      ! end if

      ! block
      !    real(ar)   :: usl_x, usl_y, usl_z
      !    usl_x = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) )
      !    usl_y = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) )
      !    usl_z = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) )
      !    print*, "Max v slopes = ", usl_x, usl_y, usl_z
      ! end block

      
   end subroutine compute_v_slopes


   !
   ! Compute w-velocity slopes
   ! 
   subroutine compute_w_slopes ( lo, hi, w, wlo, whi, slopes, &
        & domlo, domhi, bc_klo_type, bc_khi_type ) bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: wlo(3), whi(3)

      ! Grid bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! BC types 
      integer(c_int), intent(in   ) ::  &
           & bc_klo_type(domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2), &
           & bc_khi_type(domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      
      ! Arrays
      real(ar),       intent(in   ) ::                      &
           &   w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      
      real(ar),       intent(  out) ::                           &
           & slopes(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),3)

      ! Local variables
      integer                       :: i, j, k, slope_order
      real(ar)                      :: du_l, du_c, du_r, ds
      real(ar),    parameter        :: two = 2.0_ar, three2nds = 1.5_ar  
      
      integer                       :: cen, lim, flg, frm

      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flg = 3 )
      parameter( frm = 4 )

      real(ar), allocatable :: scr(:,:)

      slope_order = 2

      if (slope_order .eq. 0) then

         slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0

      else if (slope_order .eq. 2) then
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! X direction
               du_l = w(i,j,k) - w(i-1,j,k)
               du_c = half * ( w(i+1,j,k) - w(i-1,j,k) )
               du_r = w(i+1,j,k) - w(i,j,k)
               
               slopes(i,j,k,1) = mc_limiter ( du_l, du_c, du_r )  

               ! Y direction
               du_l = w(i,j,k) - w(i,j-1,k)
               du_c = half * ( w(i,j+1,k) - w(i,j-1,k) )
               du_r = w(i,j+1,k) - w(i,j,k)

               slopes(i,j,k,2) = mc_limiter ( du_l, du_c, du_r )  

               ! z direction
               du_l = w(i,j,k) - w(i,j,k-1)
               du_c = half * ( w(i,j,k+1) - w(i,j,k-1) )
               du_r = w(i,j,k+1) - w(i,j,k)

               slopes(i,j,k,3) = mc_limiter ( du_l, du_c, du_r )  
               
            end do
         end do
      end do

      else

      allocate (scr(lo(1)-1:hi(1)+1,4))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)

            do i = lo(1)-1, hi(1)+1
               ! X direction
               du_l = two * (w(i,j,k) - w(i-1,j,k))
               du_c = half * ( w(i+1,j,k) - w(i-1,j,k) )
               du_r = two * (w(i+1,j,k) - w(i,j,k))

               scr(i,cen) = du_c
               scr(i,lim) = min(abs(du_l),abs(du_r))
               scr(i,lim) = merge(scr(i,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(i,flg) = sign(1.d0,scr(i,cen))
               scr(i,frm) = scr(i,flg) * min(abs(du_c),scr(i,lim))
            end do

            do i = lo(1), hi(1)
               ds = (4.d0/3.d0) * scr(i,cen) - &
                    (1.d0/6.d0) * (scr(i+1,frm) + scr(i-1,frm))
               slopes(i,j,k,1) = scr(i,flg)*min(abs(ds),scr(i,lim))
            end do

         end do
      end do

      deallocate (scr)
      allocate (scr(lo(2)-1:hi(2)+1,4))
      
      do k = lo(3), hi(3)
         do i = lo(1), hi(1)

            do j = lo(2)-1, hi(2)+1
               du_l = two * (w(i,j,k) - w(i,j-1,k))
               du_c = half * ( w(i,j+1,k) - w(i,j-1,k) )
               du_r = two * (w(i+1,j,k) - w(i,j,k))

               scr(j,cen) = du_c
               scr(j,lim) = min(abs(du_l),abs(du_r))
               scr(j,lim) = merge(scr(j,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(j,flg) = sign(1.d0,scr(j,cen))
               scr(j,frm) = scr(j,flg) * min(abs(du_c),scr(j,lim))
            end do

            do j = lo(2), hi(2)
               ds = (4.d0/3.d0) * scr(j,cen) - &
                    (1.d0/6.d0) * (scr(j+1,frm) + scr(j-1,frm))
               slopes(i,j,k,2) = scr(j,flg)*min(abs(ds),scr(j,lim))
            end do

         end do
      end do

      deallocate (scr)
      allocate (scr(lo(3)-1:hi(3)+1,4))
      
      do i = lo(1), hi(1)
         do j = lo(2), hi(2)

            do k = lo(3)-1, hi(3)+1
               du_l = two  * (w(i,j,k) - w(i,j,k-1) )
               du_c = half * ( w(i,j,k+1) - w(i,j,k-1) )
               du_r = two  * (w(i,j,k+1) - w(i,j,k) )

               scr(k,cen) = du_c
               scr(k,lim) = min(abs(du_l),abs(du_r))
               scr(k,lim) = merge(scr(k,lim),0.d0,du_l*du_r .gt. 0.d0)
               scr(k,flg) = sign(1.d0,scr(k,cen))
               scr(k,frm) = scr(k,flg) * min(abs(du_c),scr(k,lim))
            end do

            do k = lo(3), hi(3)
               ds = (4.d0/3.d0) * scr(k,cen) - &
                    (1.d0/6.d0) * (scr(k+1,frm) + scr(k-1,frm))
               slopes(i,j,k,3) = scr(k,flg)*min(abs(ds),scr(k,lim))
            end do

         end do
      end do

      deallocate (scr)
      end if


      ! ! 
      ! ! Compute slopes at boundary where physical BCs are imposed
      ! ! 
      ! if ( lo(3) == domlo(3) ) then

      !    k = lo(3)
         
      !    do j = lo(2), hi(2)
      !       do i = lo(1), hi(1)
      
      !          if ( ( bc_klo_type(i,j,1) == MINF_ ) .or. &
      !               ( bc_klo_type(i,j,1) == NSW_ )  .or. &
      !               ( bc_klo_type(i,j,1) == FSW_ )  ) then

      !             du_l = zero
      !             du_c = - half*w(i,j,k+2) + two*w(i,j,k+1) - three2nds*w(i,j,k) 
      !             du_r = w(i,j,k+1) - w(i,j,k)
               
      !             slopes(i,j,k,3) = mc_limiter ( du_l, du_c, du_r )  
                  
      !          end if
      !       end do
      !    end do
      ! end if

      ! if ( hi(3) == (domhi(3)+1) ) then

      !    k = hi(3) 
        
      !    do j = lo(2), hi(2)
      !       do i = lo(1), hi(1)
               
      !          if ( ( bc_klo_type(i,j,1) == MINF_ ) .or. &
      !               ( bc_klo_type(i,j,1) == NSW_ )  .or. &
      !               ( bc_klo_type(i,j,1) == FSW_ )  ) then

      !             du_l = w(i,j,k) - w(i,j,k-1)
      !             du_c = half*w(i,j,k-2) - two*w(i,j,k-1) + three2nds*w(i,j,k) 
      !             du_r = zero
               
      !             slopes(i,j,k,3) = mc_limiter ( du_l, du_c, du_r )  
                  
      !          end if
      !       end do
      !    end do
      ! end if

      ! block
      !    real(ar)   :: usl_x, usl_y, usl_z
      !    usl_x = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) )
      !    usl_y = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) )
      !    usl_z = maxval ( slopes(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) )
      !    print*, "Max w slopes = ", usl_x, usl_y, usl_z
      ! end block
      
   end subroutine compute_w_slopes

   
   
   !
   ! Monotonized Central (MC) limiter
   ! 
   function mc_limiter ( dleft, dcenter, dright )  result (slope)

      real(ar), intent(in   ) :: dleft, dcenter, dright
      real(ar)                :: slope
      real(ar), parameter     :: two = 2.0_ar
     
      slope = min ( abs(two * dleft), abs(dcenter), abs(two * dright) )
      slope = merge ( slope, zero, dleft*dright >= zero )
      slope = sign ( one, dcenter ) * slope 
     
   end function mc_limiter

end module
