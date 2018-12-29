subroutine construct_gradp(slo, shi, p_g, p0_g, &
                           gpx, xlo, xhi, gpy, ylo, yhi, gpz, zlo, zhi, &
                           dx, dy, dz, &
                           bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                           bc_klo_type, bc_khi_type, domlo, domhi, ng) &
                           bind(C, name="construct_gradp")

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use bc, only: PINF_,POUT_,MINF_,FSW_,NSW_

   implicit none

   integer(c_int), intent(in   ) :: slo(3),shi(3), &
        xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
        domlo(3), domhi(3), ng

   real(rt), intent(inout) :: &
        p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
       p0_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   real(rt), intent(inout) :: &
        gpx(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)), &
        gpy(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3)), &
        gpz(zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))

   real(rt),     intent(in   ) :: dx, dy, dz

   integer(c_int), intent(in   ) :: &
           bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

   ! Local variables
   !---------------------------------------------------------------------//
   integer ::  i, j, k
   integer ::  imin, imax, jmin, jmax, kmin, kmax
   real(rt) :: odx, ody, odz
   !......................................................................!

   odx = 1.0d0/dx
   ody = 1.0d0/dy
   odz = 1.0d0/dz

   imin = max(xlo(1),slo(1)+1)
   imax = min(xhi(1),shi(1)  )
   do k = xlo(3),xhi(3)
   do j = xlo(2),xhi(2)
      do i = imin,imax
         gpx(i,j,k) = (p_g(i,j,k) + p0_g(i,j,k) - p_g(i-1,j,k) - p0_g(i-1,j,k)) * odx
      end do
   end do
   end do

   jmin = max(ylo(2),slo(2)+1)
   jmax = min(yhi(2),shi(2)  )
   do k = ylo(3),yhi(3)
   do i = ylo(1),yhi(1)
      do j = jmin, jmax
         gpy(i,j,k) = (p_g(i,j,k) + p0_g(i,j,k) - p_g(i,j-1,k) - p0_g(i,j-1,k)) * ody
      end do
   end do
   end do

   kmin = max(zlo(3),slo(3)+1)
   kmax = min(zhi(3),shi(3)  )
   do j = zlo(2),zhi(2)
   do i = zlo(1),zhi(1)
      do k = kmin,kmax
         gpz(i,j,k) = (p_g(i,j,k) + p0_g(i,j,k) - p_g(i,j,k-1) - p0_g(i,j,k-1)) * odz
      end do
   end do
   end do

end subroutine construct_gradp

subroutine set_gradp_bcs ( slo, shi, gpx, ulo, uhi, gpy, vlo, vhi, gpz, wlo, whi, &
     & bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi,               &
     & domlo, domhi, ng ) bind(C)

   use amrex_fort_module,  only: ar => amrex_real
   use iso_c_binding ,     only: c_int
   use bc
   use param,              only: zero

   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)
   integer(c_int), intent(in   ) :: ulo(3), uhi(3)
   integer(c_int), intent(in   ) :: vlo(3), vhi(3)
   integer(c_int), intent(in   ) :: wlo(3), whi(3)

   ! Grid bounds
   integer(c_int), intent(in   ) :: domlo(3), domhi(3)

   ! Number of ghost nodes
   integer(c_int), intent(in   ) :: ng

   ! BCs type
   integer(c_int), intent(in   ) :: &
        bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
        bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

   ! Arrays
   real(ar),      intent(inout) ::  &
        gpx(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),      &
        gpy(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),      &
        gpz(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

   ! Local variables
   integer  :: i, j, k
   integer  :: nlft, nrgt, nbot, ntop, nup, ndwn

   nlft = max(0,domlo(1)-slo(1))
   nbot = max(0,domlo(2)-slo(2))
   ndwn = max(0,domlo(3)-slo(3))

   nrgt = max(0,shi(1)-domhi(1))
   ntop = max(0,shi(2)-domhi(2))
   nup  = max(0,shi(3)-domhi(3))

   if (nlft .gt. 0) then
      do k = slo(3), shi(3)
         do j = slo(2), shi(2)

            select case (bct_ilo(j,k,1))

            case ( pinf_, pout_)

               gpx(ulo(1):domlo(1)-1,j,k) = huge(1.0_rt) ! Make sure we never use this
               gpy(vlo(1):domlo(1)-1,j,k) = zero
               gpz(wlo(1):domlo(1)-1,j,k) = zero

            case ( minf_)

               gpx(       domlo(1)  ,j,k) = gpx(domlo(1)+1,j,k)
               gpy(vlo(1):domlo(1)-1,j,k) = zero
               gpz(wlo(1):domlo(1)-1,j,k) = zero

            case ( nsw_)

               gpx(ulo(1):domlo(1)  ,j,k) =  zero
               gpy(vlo(1):domlo(1)-1,j,k) = -gpy(domlo(1),j,k)
               gpz(wlo(1):domlo(1)-1,j,k) = -gpz(domlo(1),j,k)

            case ( fsw_)

               gpx(ulo(1):domlo(1)  ,j,k) = zero
               gpy(vlo(1):domlo(1)-1,j,k) = gpy(domlo(1),j,k)
               gpz(wlo(1):domlo(1)-1,j,k) = gpz(domlo(1),j,k)

            end select

         end do
      end do
   endif

   if (nrgt .gt. 0) then

      do k = slo(3),shi(3)
         do j = slo(2),shi(2)

            select case ( bct_ihi(j,k,1) )

            case ( pinf_, pout_ )

               gpx(domhi(1)+2:uhi(1),j,k) =  huge(1.0_rt) ! Make sure we never use this
               gpy(domhi(1)+1:vhi(1),j,k) =  gpy(domhi(1)  ,j,k)
               gpz(domhi(1)+1:whi(1),j,k) =  gpz(domhi(1)  ,j,k)

            case ( minf_ )

               gpx(domhi(1)+1       ,j,k) = gpx(domhi(1),j,k)
               gpy(domhi(1)+1:vhi(1),j,k) = zero
               gpz(domhi(1)+1:whi(1),j,k) = zero

            case ( nsw_ )

               gpx(domhi(1)+1:uhi(1),j,k) =  zero
               gpy(domhi(1)+1:vhi(1),j,k) = -gpy(domhi(1),j,k)
               gpz(domhi(1)+1:whi(1),j,k) = -gpz(domhi(1),j,k)

            case ( fsw_ )

               gpx(domhi(1)+1:uhi(1),j,k) = zero
               gpy(domhi(1)+1:vhi(1),j,k) = gpy(domhi(1),j,k)
               gpz(domhi(1)+1:whi(1),j,k) = gpz(domhi(1),j,k)

            end select

         end do
      end do
   endif

   if (nbot .gt. 0) then

      do k = slo(3), shi(3)
         do i = slo(1), shi(1)

            select case ( bct_jlo(i,k,1) )

            case ( pinf_, pout_)

               gpx(i,ulo(2):domlo(2)-1,k) = zero
               gpy(i,vlo(2):domlo(2)-1,k) =huge(1.0_rt) ! Make sure we never use this
               gpz(i,wlo(2):domlo(2)-1,k) = zero

            case ( minf_ )

               gpx(i,vlo(2):domlo(2)-1,k) = zero
               gpy(i,       domlo(2)  ,k) = gpy(i, domlo(2)+1,k)
               gpz(i,wlo(2):domlo(2)-1,k) = zero

            case ( nsw_ )

               gpx(i,ulo(2):domlo(2)-1,k) = -gpx(i,domlo(2),k)
               gpy(i,vlo(2):domlo(2)  ,k) =  zero
               gpz(i,wlo(2):domlo(2)-1,k) = -gpz(i,domlo(2),k)

            case ( fsw_)

               gpx(i,ulo(2):domlo(2)-1,k) = gpx(i,domlo(2),k)
               gpy(i,vlo(2):domlo(2)  ,k) = 0.0d0
               gpz(i,wlo(2):domlo(2)-1,k) = gpz(i,domlo(2),k)

            end select

         end do
      end do
   endif

   if (ntop .gt. 0) then

      do k = slo(3), shi(3)
         do i = slo(1), shi(1)

            select case ( bct_jhi(i,k,1) )

            case ( pinf_, pout_ )

               gpx(i,domhi(2)+1:uhi(2),k) = zero
               gpy(i,domhi(2)+2:vhi(2),k) = huge(1.0_rt) ! Make sure we never use this
               gpz(i,domhi(2)+1:whi(2),k) = zero

            case ( minf_)

               gpx(i,domhi(2)+1:uhi(2),k) = zero
               gpy(i,domhi(2)+1       ,k) = gpy(i,domhi(2),k)
               gpz(i,domhi(2)+1:whi(2),k) = zero

            case ( nsw_)

               gpx(i,domhi(2)+1:uhi(2),k) = -gpx(i,domhi(2),k)
               gpy(i,domhi(2)+1:vhi(2),k) =  zero
               gpz(i,domhi(2)+1:whi(2),k) = -gpz(i,domhi(2),k)

            case ( fsw_)

               gpx(i,domhi(2)+1:uhi(2),k) = gpx(i,domhi(2),k)
               gpy(i,domhi(2)+1:vhi(2),k) = zero
               gpz(i,domhi(2)+1:whi(2),k) = gpz(i,domhi(2),k)

            end select
         end do
      end do
   endif

   if (ndwn .gt. 0) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_ )

               gpx(i,j,ulo(3):domlo(3)-1) = zero
               gpy(i,j,vlo(3):domlo(3)-1) = zero
               gpz(i,j,wlo(3):domlo(3)-1) = huge(1.0_rt) ! Make sure we never use this

            case ( minf_ )

               gpx(i,j,ulo(3):domlo(3)-1) = zero
               gpy(i,j,vlo(3):domlo(3)-1) = zero
               gpz(i,j,       domlo(3)  ) = gpz(i,j,domlo(3)+1)

            case ( nsw_ )

               gpx(i,j,ulo(3):domlo(3)-1) = -gpx(i,j,domlo(3))
               gpy(i,j,vlo(3):domlo(3)-1) = -gpy(i,j,domlo(3))
               gpz(i,j,wlo(3):domlo(3)  ) =  zero

            case ( fsw_ )

               gpx(i,j,ulo(3):domlo(3)-1) = gpx(i,j,domlo(3))
               gpy(i,j,vlo(3):domlo(3)-1) = gpy(i,j,domlo(3))
               gpz(i,j,wlo(3):domlo(3)  ) = zero

            end select
         end do
      end do
   endif

   if (nup .gt. 0) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)

            select case ( bct_khi(i,j,1) )

            case ( pinf_, pout_ )

               gpx(i,j,domhi(3)+1:uhi(3)) = zero
               gpy(i,j,domhi(3)+1:vhi(3)) = zero
               gpz(i,j,domhi(3)+2:whi(3)) = huge(1.0_rt) ! Make sure we never use this

            case ( minf_ )

               gpx(i,j,domhi(3)+1:uhi(3)) = zero
               gpy(i,j,domhi(3)+1:vhi(3)) = zero
               gpz(i,j,domhi(3)+1       ) = gpz(i,j,domhi(3))

            case ( nsw_ )

               gpx(i,j,domhi(3)+1:uhi(3)) = -gpx(i,j,domhi(3))
               gpy(i,j,domhi(3)+1:vhi(3)) = -gpy(i,j,domhi(3))
               gpz(i,j,domhi(3)+1:whi(3)) =  zero

            case ( fsw_ )

               gpx(i,j,domhi(3)+1:uhi(3)) = gpx(i,j,domhi(3))
               gpy(i,j,domhi(3)+1:vhi(3)) = gpy(i,j,domhi(3))
               gpz(i,j,domhi(3)+1:whi(3)) = zero

            end select
         end do
      end do
   endif

end subroutine set_gradp_bcs
