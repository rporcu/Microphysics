subroutine construct_gradp(slo, shi, p_g, p0_g, &
                           gpx, xlo, xhi, gpy, ylo, yhi, gpz, zlo, zhi, &
                           dx, dy, dz, &
                           bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                           bc_klo_type, bc_khi_type, domlo, domhi, ng) &
                           bind(C, name="construct_gradp")

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use bc, only: PINF_,POUT_,MINF_,FSW_,NSW_

   implicit none

   integer(c_int), intent(in   ) :: slo(3),shi(3), &
        xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
        domlo(3), domhi(3), ng

   real(c_real), intent(inout) :: &
        p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
       p0_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   real(c_real), intent(inout) :: &
        gpx(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)), &
        gpy(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3)), &
        gpz(zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))

   real(c_real),     intent(in   ) :: dx, dy, dz

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
   real(c_real) :: odx, ody, odz
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
      if (xlo(1).le.domlo(1)) then
         if(bc_ilo_type(j,k,1) == PINF_ .or. &
            bc_ilo_type(j,k,1) == POUT_) then

            gpy(domlo(1)-1,j,k) = 0.d0
            gpz(domlo(1)-1,j,k) = 0.d0

         else if (bc_ilo_type(j,k,1) == MINF_) then

            gpx(domlo(1)  ,j,k) = gpx(domlo(1)+1,j,k)
            gpy(domlo(1)-1,j,k) = 0.0d0
            gpz(domlo(1)-1,j,k) = 0.0d0

         else if (bc_ilo_type(j,k,1) == NSW_) then

            gpx(domlo(1)  ,j,k) =  0.0d0
            gpy(domlo(1)-1,j,k) = -gpy(domlo(1),j,k)
            gpz(domlo(1)-1,j,k) = -gpz(domlo(1),j,k)

         else if (bc_ilo_type(j,k,1) == FSW_) then

            gpx(domlo(1)  ,j,k) = 0.0d0
            gpy(domlo(1)-1,j,k) = gpy(domlo(1),j,k)
            gpz(domlo(1)-1,j,k) = gpz(domlo(1),j,k)

         end if
      end if
      if (xhi(1).ge.domhi(1)+1) then
         if(bc_ihi_type(j,k,1) == PINF_ .or. &
            bc_ihi_type(j,k,1) == POUT_) then

            gpy(domhi(1)+1,j,k) = 0.d0
            gpz(domhi(1)+1,j,k) = 0.d0

         else if (bc_ihi_type(j,k,1) == MINF_) then

            gpx(domhi(1)+1,j,k) = gpx(domhi(1),j,k)
            gpy(domhi(1)+1,j,k) = 0.0d0
            gpz(domhi(1)+1,j,k) = 0.0d0

         else if (bc_ihi_type(j,k,1) == NSW_) then

            gpx(domhi(1)+1,j,k) =  0.0d0
            gpy(domhi(1)+1,j,k) = -gpy(domhi(1),j,k)
            gpz(domhi(1)+1,j,k) = -gpz(domhi(1),j,k)

         else if (bc_ihi_type(j,k,1) == FSW_) then

            gpx(domhi(1)+1,j,k) = 0.0d0 
            gpy(domhi(1)+1,j,k) = gpy(domhi(1),j,k) 
            gpz(domhi(1)+1,j,k) = gpz(domhi(1),j,k) 
         end if
      end if
   end do
   end do

   jmin = max(ylo(2),slo(2)+1)
   jmax = min(yhi(2),shi(2)  )
   do k = ylo(3),yhi(3)
   do i = ylo(1),yhi(1)
      do j = jmin, jmax
         gpy(i,j,k) = (p_g(i,j,k) + p0_g(i,j,k) - p_g(i,j-1,k) - p0_g(i,j-1,k)) * ody
      end do
      if (ylo(2).le.domlo(2)) then
         if(bc_jlo_type(i,k,1) == PINF_ .or. &
            bc_jlo_type(i,k,1) == POUT_) then

            gpx(i,domlo(2)-1,k) = 0.d0
            gpz(i,domlo(2)-1,k) = 0.d0

         else if (bc_jlo_type(i,k,1) == MINF_) then

            gpx(i,domlo(2)-1,k) = 0.0d0
            gpy(i,domlo(2)  ,k) = gpy(i,domlo(2)+1,k)
            gpz(i,domlo(2)-1,k) = 0.0d0

         else if (bc_jlo_type(i,k,1) == NSW_) then

            gpx(i,domlo(2)-1,k) = -gpx(i,domlo(2),k)
            gpy(i,domlo(2)  ,k) =  0.0d0
            gpz(i,domlo(2)-1,k) = -gpz(i,domlo(2),k)

         else if (bc_jlo_type(i,k,1) == FSW_) then

            gpx(i,domlo(2)-1,k) = gpx(i,domlo(2),k)
            gpy(i,domlo(2)  ,k) = 0.0d0
            gpz(i,domlo(2)-1,k) = gpz(i,domlo(2),k)

         end if
      end if
      if (yhi(2).ge.domhi(2)+1) then
         if(bc_jhi_type(i,k,1) == PINF_ .or. &
            bc_jhi_type(i,k,1) == POUT_) then

            gpx(i,domhi(2)+1,k) = 0.d0
            gpz(i,domhi(2)+1,k) = 0.d0

         else if (bc_jhi_type(i,k,1) == MINF_) then

            gpx(i,domhi(2)+1,k) = 0.0d0
            gpy(i,domhi(2)+1,k) = gpy(i,domhi(2),k)
            gpz(i,domhi(2)+1,k) = 0.0d0

         else if (bc_jhi_type(i,k,1) == NSW_) then

            gpx(i,domhi(2)+1,k) = -gpx(i,domhi(2),k)
            gpy(i,domhi(2)+1,k) =  0.0d0
            gpz(i,domhi(2)+1,k) = -gpz(i,domhi(2),k)

         else if (bc_jhi_type(i,k,1) == FSW_) then

            gpx(i,domhi(2)+1,k) = gpx(i,domhi(2),k) 
            gpy(i,domhi(2)+1,k) = 0.0d0 
            gpz(i,domhi(2)+1,k) = gpz(i,domhi(2),k) 
         end if
      end if
   end do
   end do

   kmin = max(zlo(3),slo(3)+1)
   kmax = min(zhi(3),shi(3)  )
   do j = zlo(2),zhi(2)
   do i = zlo(1),zhi(1)
      do k = kmin,kmax
         gpz(i,j,k) = (p_g(i,j,k) + p0_g(i,j,k) - p_g(i,j,k-1) - p0_g(i,j,k-1)) * odz
      end do
      if (zlo(3).le.domlo(3)) then
         if(bc_klo_type(i,j,1) == PINF_ .or. &
            bc_klo_type(i,j,1) == POUT_) then

            gpx(i,j,domlo(3)-1) = 0.d0
            gpy(i,j,domlo(3)-1) = 0.d0

         else if (bc_klo_type(i,j,1) == MINF_) then

            gpx(i,j,domlo(3)-1) = 0.0d0
            gpy(i,j,domlo(3)-1) = 0.0d0
            gpz(i,j,domlo(3)  ) = gpz(i,j,domlo(3)+1)

         else if (bc_klo_type(i,j,1) == NSW_) then

            gpx(i,j,domlo(3)-1) = -gpx(i,j,domlo(3))
            gpy(i,j,domlo(3)-1) = -gpy(i,j,domlo(3))
            gpz(i,j,domlo(3)  ) =  0.0d0

         else if (bc_klo_type(i,j,1) == FSW_) then

            gpx(i,j,domlo(3)-1) = gpx(i,j,domlo(3))
            gpy(i,j,domlo(3)-1) = gpy(i,j,domlo(3))
            gpz(i,j,domlo(3)  ) = 0.0d0

         end if
      end if
      if (zhi(3).ge.domhi(3)+1) then
         if(bc_khi_type(i,j,1) == PINF_ .or. &
            bc_khi_type(i,j,1) == POUT_) then

            gpx(i,j,domhi(3)+1) = 0.d0
            gpy(i,j,domhi(3)+1) = 0.d0

         else if (bc_khi_type(i,j,1) == MINF_) then

            gpx(i,j,domhi(3)+1) = 0.0d0
            gpy(i,j,domhi(3)+1) = 0.0d0
            gpz(i,j,domhi(3)+1) = gpz(i,j,domhi(3))

         else if (bc_khi_type(i,j,1) == NSW_) then

            gpx(i,j,domhi(3)+1) = -gpx(i,j,domhi(3))
            gpy(i,j,domhi(3)+1) = -gpy(i,j,domhi(3))
            gpz(i,j,domhi(3)+1) =  0.0d0

         else if (bc_khi_type(i,j,1) == FSW_) then

            gpx(i,j,domhi(3)+1) = gpx(i,j,domhi(3))
            gpy(i,j,domhi(3)+1) = gpy(i,j,domhi(3))
            gpz(i,j,domhi(3)+1) = 0.0d0

         end if
      end if
   end do
   end do

end subroutine construct_gradp
