! 
!              
!  This module contains the subroutines to compute the three components
!  of the convection term (U.grad)U
!
! 
module convection_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   use bc,                only: minf_, nsw_, fsw_, psw_, pinf_, pout_

   implicit none
   private

   ! Public members
   public compute_ugradu

contains

   !
   ! Compute all components of u dot grad u at cell centers
   ! 
   subroutine compute_ugradu ( lo, hi, vel, vlo, vhi, xslopes, yslopes, zslopes, slo, shi, &
                               ugradu, ulo, uhi, domlo, domhi, &
                               bc_ilo_type, bc_ihi_type, &
                               bc_jlo_type, bc_jhi_type, &
                               bc_klo_type, bc_khi_type, ng, dx ) bind(C)


      ! Tile bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: domlo(3), domhi(3), ng

      ! Grid
      real(ar),        intent(in   ) :: dx(3)

      ! Velocity Array
      real(ar),        intent(in   ) ::                            &
           & vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)    , &
           & xslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & yslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & zslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(ar),        intent(  out) ::                           &
           & ugradu(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Local variables
      integer(c_int)                 :: i, j, k
      real(ar)                       :: idx, idy, idz
      real(ar)                       :: udu, vdu, wdu
      real(ar)                       :: udv, vdv, wdv
      real(ar)                       :: udw, vdw, wdw
      real(ar)                       :: upls, umns, vpls, vmns, wpls, wmns 
      real(ar)                       :: u_e, u_w, u_s, u_n, u_b, u_t
      real(ar)                       :: v_e, v_w, v_s, v_n, v_b, v_t
      real(ar)                       :: w_e, w_w, w_s, w_n, w_b, w_t

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! ****************************************************
               ! West face
               ! ****************************************************

               ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
               ! In the case of PINF, POUT          we are using the upwind value
               if (i.eq.domlo(1) .and. &
                     ( ( bc_ilo_type(j,k,1) == MINF_ ) .or. &
                       ( bc_ilo_type(j,k,1) == NSW_ )  .or. &
                       ( bc_ilo_type(j,k,1) == FSW_ )  .or. &
                       ( bc_ilo_type(j,k,1) == PSW_ )  .or. &
                       ( bc_ilo_type(j,k,1) == PINF_ )  .or. &
                       ( bc_ilo_type(j,k,1) == POUT_ )  ) ) then
   
                        u_w =  vel(i-1,j,k,1)
                        v_w =  vel(i-1,j,k,2)
                        w_w =  vel(i-1,j,k,3)

               else

                   upls  = vel(i  ,j,k,1) - half * xslopes(i  ,j,k,1)
                   umns  = vel(i-1,j,k,1) + half * xslopes(i-1,j,k,1)
                   vpls  = vel(i  ,j,k,2) - half * xslopes(i  ,j,k,2)
                   vmns  = vel(i-1,j,k,2) + half * xslopes(i-1,j,k,2)
                   wpls  = vel(i  ,j,k,3) - half * xslopes(i  ,j,k,3)
                   wmns  = vel(i-1,j,k,3) + half * xslopes(i-1,j,k,3)

                   u_w   = upwind_normal ( umns, upls )
                   v_w   = upwind        ( vmns, vpls, u_w)
                   w_w   = upwind        ( wmns, wpls, u_w)
    
               endif 

               ! ****************************************************
               ! East face
               ! ****************************************************

               ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
               ! In the case of PINF, POUT          we are using the upwind value
               if (i.eq.domhi(1) .and. &
                     ( ( bc_ihi_type(j,k,1) == MINF_ ) .or. &
                       ( bc_ihi_type(j,k,1) == NSW_  ) .or. &
                       ( bc_ihi_type(j,k,1) == FSW_  ) .or. &
                       ( bc_ihi_type(j,k,1) == PSW_  ) .or. &
                       ( bc_ihi_type(j,k,1) == PINF_ ) .or. &
                       ( bc_ihi_type(j,k,1) == POUT_ )  ) ) then
   
                        u_e =  vel(i+1,j,k,1)
                        v_e =  vel(i+1,j,k,2)
                        w_e =  vel(i+1,j,k,3)

               else

                   upls  = vel(i+1,j,k,1) - half * xslopes(i+1,j,k,1)
                   umns  = vel(i  ,j,k,1) + half * xslopes(i  ,j,k,1)
                   vpls  = vel(i+1,j,k,2) - half * xslopes(i+1,j,k,2)
                   vmns  = vel(i  ,j,k,2) + half * xslopes(i  ,j,k,2)
                   wpls  = vel(i+1,j,k,3) - half * xslopes(i+1,j,k,3)
                   wmns  = vel(i  ,j,k,3) + half * xslopes(i  ,j,k,3)

                   u_e   = upwind_normal ( umns, upls )
                   v_e   = upwind        ( vmns, vpls, u_e)
                   w_e   = upwind        ( wmns, wpls, u_e)

               endif

               ! ****************************************************
               ! South face
               ! ****************************************************

               ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
               ! In the case of PINF, POUT          we are using the upwind value
               if (j.eq.domlo(2) .and. &
                     ( ( bc_jlo_type(i,k,1) == MINF_ ) .or. &
                       ( bc_jlo_type(i,k,1) == NSW_ )  .or. &
                       ( bc_jlo_type(i,k,1) == FSW_ )  .or. &
                       ( bc_jlo_type(i,k,1) == PSW_ )  .or. &
                       ( bc_jlo_type(i,k,1) == PINF_ )  .or. &
                       ( bc_jlo_type(i,k,1) == POUT_ )  )  ) then
   
                        u_s =  vel(i,j-1,k,1)
                        v_s =  vel(i,j-1,k,2)
                        w_s =  vel(i,j-1,k,3)

               else

                   upls  = vel(i,j  ,k,1) - half * yslopes(i,j  ,k,1)
                   umns  = vel(i,j-1,k,1) + half * yslopes(i,j-1,k,1)
                   vpls  = vel(i,j  ,k,2) - half * yslopes(i,j  ,k,2)
                   vmns  = vel(i,j-1,k,2) + half * yslopes(i,j-1,k,2)
                   wpls  = vel(i,j  ,k,3) - half * yslopes(i,j  ,k,3)
                   wmns  = vel(i,j-1,k,3) + half * yslopes(i,j-1,k,3)

                   v_s   = upwind_normal ( vmns, vpls )
                   u_s   = upwind        ( umns, upls, v_s )
                   w_s   = upwind        ( wmns, wpls, v_s )

               endif

               ! ****************************************************
               ! North face
               ! ****************************************************

               ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
               ! In the case of PINF, POUT          we are using the upwind value
               if (j.eq.domhi(2) .and. &
                     ( ( bc_jhi_type(i,k,1) == MINF_ ) .or. &
                       ( bc_jhi_type(i,k,1) == NSW_ )  .or. &
                       ( bc_jhi_type(i,k,1) == FSW_ )  .or. &
                       ( bc_jhi_type(i,k,1) == PSW_ )  .or. &
                       ( bc_jhi_type(i,k,1) == PINF_ )  .or. &
                       ( bc_jhi_type(i,k,1) == POUT_ )  ) ) then
   
                        u_n =  vel(i,j+1,k,1)
                        v_n =  vel(i,j+1,k,2)
                        w_n =  vel(i,j+1,k,3)

               else

                   upls  = vel(i,j+1,k,1) - half * yslopes(i,j+1,k,1)
                   umns  = vel(i,j  ,k,1) + half * yslopes(i,j  ,k,1)
                   vpls  = vel(i,j+1,k,2) - half * yslopes(i,j+1,k,2)
                   vmns  = vel(i,j  ,k,2) + half * yslopes(i,j  ,k,2)
                   wpls  = vel(i,j+1,k,3) - half * yslopes(i,j+1,k,3)
                   wmns  = vel(i,j  ,k,3) + half * yslopes(i,j  ,k,3)

                   v_n   = upwind_normal ( vmns, vpls )
                   u_n   = upwind        ( umns, upls, v_n)
                   w_n   = upwind        ( wmns, wpls, v_n)

               endif

               ! ****************************************************
               ! Bottom face
               ! ****************************************************

               ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
               ! In the case of PINF, POUT          we are using the upwind value
               if (k.eq.domlo(3) .and. &
                     ( ( bc_klo_type(i,j,1) == MINF_ ) .or. &
                       ( bc_klo_type(i,j,1) == NSW_ )  .or. &
                       ( bc_klo_type(i,j,1) == FSW_ )  .or. &
                       ( bc_klo_type(i,j,1) == PSW_ )  .or. &
                       ( bc_klo_type(i,j,1) == PINF_ )  .or. &
                       ( bc_klo_type(i,j,1) == POUT_ )  ) ) then
   
                        u_b =  vel(i,j,k-1,1)
                        v_b =  vel(i,j,k-1,2)
                        w_b =  vel(i,j,k-1,3)

               else

                   upls  = vel(i,j,k  ,1) - half * zslopes(i,j,k  ,1)
                   umns  = vel(i,j,k-1,1) + half * zslopes(i,j,k-1,1)
                   vpls  = vel(i,j,k  ,2) - half * zslopes(i,j,k  ,2)
                   vmns  = vel(i,j,k-1,2) + half * zslopes(i,j,k-1,2)
                   wpls  = vel(i,j,k  ,3) - half * zslopes(i,j,k  ,3)
                   wmns  = vel(i,j,k-1,3) + half * zslopes(i,j,k-1,3)

                   w_b   = upwind_normal ( wmns, wpls )
                   u_b   = upwind        ( umns, upls, w_b )
                   v_b   = upwind        ( vmns, vpls, w_b )

               endif

               ! ****************************************************
               ! Top face
               ! ****************************************************

               ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
               ! In the case of PINF, POUT          we are using the upwind value
               if (k.eq.domhi(3) .and. &
                     ( ( bc_khi_type(i,j,1) == MINF_ ) .or. &
                       ( bc_khi_type(i,j,1) == NSW_ )  .or. &
                       ( bc_khi_type(i,j,1) == FSW_ )  .or. &
                       ( bc_khi_type(i,j,1) == PSW_ )  .or. &
                       ( bc_khi_type(i,j,1) == PINF_ )  .or. &
                       ( bc_khi_type(i,j,1) == POUT_ )  ) ) then
   
                        u_t =  vel(i,j,k+1,1)
                        v_t =  vel(i,j,k+1,2)
                        w_t =  vel(i,j,k+1,3)

               else

                   upls  = vel(i,j,k+1,1) - half * zslopes(i,j,k+1,1)
                   umns  = vel(i,j,k  ,1) + half * zslopes(i,j,k  ,1)
                   vpls  = vel(i,j,k+1,2) - half * zslopes(i,j,k+1,2)
                   vmns  = vel(i,j,k  ,2) + half * zslopes(i,j,k  ,2)
                   wpls  = vel(i,j,k+1,3) - half * zslopes(i,j,k+1,3)
                   wmns  = vel(i,j,k  ,3) + half * zslopes(i,j,  k,3)

                   w_t   = upwind_normal ( wmns, wpls )
                   u_t   = upwind        ( umns, upls, w_t)
                   v_t   = upwind        ( vmns, vpls, w_t)

               endif

               ! ****************************************************
               ! Define convective terms
               ! ****************************************************

               udu   = vel(i,j,k,1) * (u_e - u_w)
               udv   = vel(i,j,k,1) * (v_e - v_w)
               udw   = vel(i,j,k,1) * (w_e - w_w)

               vdu   = vel(i,j,k,2) * (u_n - u_s)
               vdv   = vel(i,j,k,2) * (v_n - v_s)
               vdw   = vel(i,j,k,2) * (w_n - w_s)

               wdu   = vel(i,j,k,3) * (u_t - u_b)
               wdv   = vel(i,j,k,3) * (v_t - v_b)
               wdw   = vel(i,j,k,3) * (w_t - w_b)

               ! ****************************************************
               ! Assemble terms
               ! ****************************************************
               ugradu(i,j,k,1) = udu*idx + vdu*idy + wdu*idz
               ugradu(i,j,k,2) = udv*idx + vdv*idy + wdv*idz
               ugradu(i,j,k,3) = udw*idx + vdw*idy + wdw*idz

            end do
         end do
      end do

   end subroutine compute_ugradu

   ! Upwind along direction normal to velocity component
   function upwind_normal ( umns, upls ) result (ev)

      real(ar), intent(in) :: umns, upls
      real(ar)             :: ev, avg

      if ( umns < zero .and. upls > zero ) then
         ev = zero
      else 
         avg = half * ( upls + umns )
         ev = merge ( umns, upls, avg >= zero ) 
      end if

   end function upwind_normal

   ! Upwind non-normal velocity
   function upwind ( velmns, velpls, uedge ) result (ev)

      ! Small value to protect against tiny velocities used in upwinding
      real(ar),        parameter     :: small_vel = 1.0d-10

      real(ar), intent(in) :: velmns, velpls, uedge
      real(ar)             :: ev

      if ( abs(uedge) .lt. small_vel) then
         ev = half * ( velpls + velmns )
      else 
         ev = merge ( velmns, velpls, uedge >= zero ) 
      end if

   end function upwind

end module convection_mod
