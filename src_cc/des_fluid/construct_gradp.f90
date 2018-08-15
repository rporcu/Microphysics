subroutine construct_gradp(lo, hi, p_g, p0_g, rlo, rhi, &
                           gpx, xlo, xhi, gpy, ylo, yhi, gpz, zlo, zhi, &
                           dx, dy, dz, &
                           bct_ilo, bct_ihi, bct_jlo, bct_jhi, &
                           bct_klo, bct_khi, domlo, domhi, ng, nodal_pressure) &
                           bind(C, name="construct_gradp")

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use bc, only: PINF_,POUT_,MINF_,FSW_,NSW_

   implicit none

   integer(c_int), intent(in   ) :: lo(3), hi(3), rlo(3), rhi(3), &
        xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
        domlo(3), domhi(3), ng, nodal_pressure

   real(rt), intent(inout) :: &
        p_g(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3)), &
       p0_g(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

   real(rt), intent(inout) :: &
        gpx(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)), &
        gpy(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3)), &
        gpz(zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))

   real(rt),     intent(in   ) :: dx, dy, dz

   integer(c_int), intent(in   ) :: &
           bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

   ! Local variables
   !---------------------------------------------------------------------//
   integer ::  i, j, k
   integer ::  imin, imax, jmin, jmax, kmin, kmax
   real(rt) :: odx, ody, odz
   !......................................................................!

   odx = 1.0d0/dx
   ody = 1.0d0/dy
   odz = 1.0d0/dz

   if (nodal_pressure .eq. 1) then

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

           gpx(i,j,k) = 0.25d0 * odx * ( &
                 p_g(i+1,j,k) +  p_g(i+1,j+1,k) +  p_g(i+1,j,k+1) +  p_g(i+1,j+1,k+1) &
              + p0_g(i+1,j,k) + p0_g(i+1,j+1,k) + p0_g(i+1,j,k+1) + p0_g(i+1,j+1,k+1) &
              -  p_g(i  ,j,k) -  p_g(i  ,j+1,k) -  p_g(i  ,j,k+1) -  p_g(i  ,j+1,k+1) &
              - p0_g(i  ,j,k) - p0_g(i  ,j+1,k) - p0_g(i  ,j,k+1) - p0_g(i  ,j+1,k+1)  ) 
           

           gpy(i,j,k) = 0.25d0 * ody * ( &
                 p_g(i,j+1,k) +  p_g(i+1,j+1,k) +  p_g(i,j+1,k+1) +  p_g(i+1,j+1,k+1) & 
              + p0_g(i,j+1,k) + p0_g(i+1,j+1,k) + p0_g(i,j+1,k+1) + p0_g(i+1,j+1,k+1) &
              -  p_g(i,j  ,k) -  p_g(i+1,j  ,k) -  p_g(i,j  ,k+1) -  p_g(i+1,j  ,k+1) &
              - p0_g(i,j  ,k) - p0_g(i+1,j  ,k) - p0_g(i,j  ,k+1) - p0_g(i+1,j  ,k+1) )

           
           gpz(i,j,k) = 0.25d0 * odz * ( &
                 p_g(i,j,k+1) +  p_g(i+1,j,k+1) +  p_g(i,j+1,k+1) +  p_g(i+1,j+1,k+1) &
              + p0_g(i,j,k+1) + p0_g(i+1,j,k+1) + p0_g(i,j+1,k+1) + p0_g(i+1,j+1,k+1) &
              -  p_g(i,j,k  ) -  p_g(i+1,j,k  ) -  p_g(i,j+1,k  ) -  p_g(i+1,j+1,k  ) &
              - p0_g(i,j,k  ) - p0_g(i+1,j,k  ) - p0_g(i,j+1,k  ) - p0_g(i+1,j+1,k  ) )

      end do
      end do
      end do

   else

      imin = max(xlo(1),rlo(1)+1)
      imax = min(xhi(1),rhi(1)  )
      do k = xlo(3),xhi(3)
      do j = xlo(2),xhi(2)
         do i = imin,imax
            gpx(i,j,k) = (p_g(i,j,k) + p0_g(i,j,k) - p_g(i-1,j,k) - p0_g(i-1,j,k)) * odx
         end do
      end do
      end do

      jmin = max(ylo(2),rlo(2)+1)
      jmax = min(yhi(2),rhi(2)  )
      do k = ylo(3),yhi(3)
      do i = ylo(1),yhi(1)
         do j = jmin, jmax
            gpy(i,j,k) = (p_g(i,j,k) + p0_g(i,j,k) - p_g(i,j-1,k) - p0_g(i,j-1,k)) * ody
         end do
      end do
      end do

      kmin = max(zlo(3),rlo(3)+1)
      kmax = min(zhi(3),rhi(3)  )
      do j = zlo(2),zhi(2)
      do i = zlo(1),zhi(1)
         do k = kmin,kmax
            gpz(i,j,k) = (p_g(i,j,k) + p0_g(i,j,k) - p_g(i,j,k-1) - p0_g(i,j,k-1)) * odz
         end do
      end do
      end do

   end if

end subroutine construct_gradp

subroutine set_gradp_bcs ( slo, shi, gpx, xlo, xhi, gpy, ylo, yhi, gpz, zlo, zhi, &
     & bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi,               &
     & domlo, domhi, ng, nodal_pressure ) bind(C) 

   use amrex_fort_module,  only: ar => amrex_real
   use iso_c_binding ,     only: c_int
   use bc

   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)
   integer(c_int), intent(in   ) :: xlo(3), xhi(3)
   integer(c_int), intent(in   ) :: ylo(3), yhi(3)
   integer(c_int), intent(in   ) :: zlo(3), zhi(3)

   ! Grid bounds
   integer(c_int), intent(in   ) :: domlo(3), domhi(3)

   ! Number of ghost nodes
   integer(c_int), intent(in   ) :: ng, nodal_pressure
   
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
        gpx(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)), &
        gpy(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3)), &
        gpz(zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))

   ! Local variables
   integer  :: i, j, k

   if (xlo(1).le.domlo(1)) then
      do k = xlo(3),xhi(3)
      do j = xlo(2),xhi(2)

            select case (bct_ilo(j,k,1))

            case ( pinf_, pout_)

               gpx(domlo(1)-1,j,k) = gpx(domlo(1),j,k)
               gpy(domlo(1)-1,j,k) = gpy(domlo(1),j,k)
               gpz(domlo(1)-1,j,k) = gpz(domlo(1),j,k)

            case ( minf_)

                if (nodal_pressure .eq. 1) then
                   gpx(domlo(1)-1,j,k) = gpx(domlo(1),j,k)
                else
                   gpx(domlo(1)  ,j,k) = gpx(domlo(1)+1,j,k)
                end if

               gpy(domlo(1)-1,j,k) = 0.0d0
               gpz(domlo(1)-1,j,k) = 0.0d0

            case ( nsw_)

               if (nodal_pressure .eq. 1) then
                  gpx(domlo(1)-1,j,k) = -gpx(domlo(1),j,k)
               else
                  gpx(domlo(1)  ,j,k) =  0.0d0
               end if
   
               gpy(domlo(1)-1,j,k) = -gpy(domlo(1),j,k)
               gpz(domlo(1)-1,j,k) = -gpz(domlo(1),j,k)

            case ( fsw_)

               if (nodal_pressure .eq. 1) then
                  gpx(domlo(1)-1,j,k) = -gpx(domlo(1),j,k)
               else
                  gpx(domlo(1)  ,j,k) =  0.0d0
               end if
   
               gpy(domlo(1)-1,j,k) = gpy(domlo(1),j,k)
               gpz(domlo(1)-1,j,k) = gpz(domlo(1),j,k)

            end select
      end do
      end do
   end if

   if (xhi(1).ge.domhi(1)+1) then
      do k = xlo(3),xhi(3)
      do j = xlo(2),xhi(2)

            select case (bct_ihi(j,k,1))

            case ( pinf_, pout_)

               gpx(domhi(1)+1,j,k) = gpx(domhi(1),j,k)
               gpy(domhi(1)+1,j,k) = gpy(domhi(1),j,k)
               gpz(domhi(1)+1,j,k) = gpz(domhi(1),j,k)

            case ( minf_)

               gpx(domhi(1)+1,j,k) = gpx(domhi(1),j,k)
               gpy(domhi(1)+1,j,k) = 0.0d0
               gpz(domhi(1)+1,j,k) = 0.0d0

            case ( nsw_)

               if (nodal_pressure .eq. 1) then
                  gpx(domhi(1)+1,j,k) = -gpx(domhi(1),j,k)
               else
                  gpx(domhi(1)+1,j,k) =  0.0d0
               end if
   
               gpy(domhi(1)+1,j,k) = -gpy(domhi(1),j,k)
               gpz(domhi(1)+1,j,k) = -gpz(domhi(1),j,k)

            case ( fsw_)

               if (nodal_pressure .eq. 1) then
                  gpx(domhi(1)+1,j,k) = -gpx(domhi(1),j,k)
               else
                  gpx(domhi(1)+1,j,k) =  0.0d0
               end if

               gpy(domhi(1)+1,j,k) = gpy(domhi(1),j,k) 
               gpz(domhi(1)+1,j,k) = gpz(domhi(1),j,k) 

            end select
      end do
      end do
   end if

   if (ylo(2).le.domlo(2)) then
      do k = ylo(3),yhi(3)
      do i = ylo(1),yhi(1)

            select case (bct_jlo(i,k,1))

            case ( pinf_, pout_)

               gpx(i,domlo(2)-1,k) = gpx(i,domlo(2),k)
               gpy(i,domlo(2)-1,k) = gpy(i,domlo(2),k)
               gpz(i,domlo(2)-1,k) = gpz(i,domlo(2),k)
   
            case ( minf_)

               if (nodal_pressure .eq. 1) then
                  gpy(i,domlo(2)-1,k) = gpy(i,domlo(2),k)
               else
                  gpy(i,domlo(2)  ,k) = gpy(i,domlo(2)+1,k)
               end if
            
               gpx(i,domlo(2)-1,k) = 0.0d0
               gpz(i,domlo(2)-1,k) = 0.0d0

            case ( nsw_)

               if (nodal_pressure .eq. 1) then
                  gpy(i,domlo(2)-1,k) = -gpy(i,domlo(2),k)
               else
                  gpy(i,domlo(2)  ,k) =  0.0d0
               end if
   
               gpx(i,domlo(2)-1,k) = -gpx(i,domlo(2),k)
               gpz(i,domlo(2)-1,k) = -gpz(i,domlo(2),k)

            case ( fsw_)

               if (nodal_pressure .eq. 1) then
                  gpy(i,domlo(2)-1,k) = -gpy(i,domlo(2),k)
               else
                  gpy(i,domlo(2)  ,k) = 0.0d0
               end if
   
               gpx(i,domlo(2)-1,k) = gpx(i,domlo(2),k)
               gpz(i,domlo(2)-1,k) = gpz(i,domlo(2),k)

            end select
      end do
      end do
   end if

   if (yhi(2).ge.domhi(2)+1) then

      do k = ylo(3),yhi(3)
      do i = ylo(1),yhi(1)

            select case (bct_jhi(i,k,1))

            case ( pinf_, pout_)

               gpx(i,domhi(2)+1,k) = gpx(i,domhi(2),k)
               gpy(i,domhi(2)+1,k) = gpy(i,domhi(2),k)
               gpz(i,domhi(2)+1,k) = gpz(i,domhi(2),k)

            case ( minf_)

               gpx(i,domhi(2)+1,k) = 0.0d0
               gpy(i,domhi(2)+1,k) = gpy(i,domhi(2),k)
               gpz(i,domhi(2)+1,k) = 0.0d0

            case ( nsw_)

               if (nodal_pressure .eq. 1) then
                  gpy(i,domhi(2)+1,k) = -gpy(i,domhi(2),k)
               else
                  gpy(i,domhi(2)+1,k) =  0.0d0
               end if
   
               gpx(i,domhi(2)+1,k) = -gpx(i,domhi(2),k)
               gpz(i,domhi(2)+1,k) = -gpz(i,domhi(2),k)

            case ( fsw_)

               if (nodal_pressure .eq. 1) then
                  gpy(i,domhi(2)+1,k) = -gpy(i,domhi(2),k)
               else
                  gpy(i,domhi(2)+1,k) = 0.0d0 
               end if

               gpx(i,domhi(2)+1,k) = gpx(i,domhi(2),k) 
               gpz(i,domhi(2)+1,k) = gpz(i,domhi(2),k) 

            end select
      end do
      end do
   end if

   if (zlo(3).le.domlo(3)) then
      do j = zlo(2),zhi(2)
      do i = zlo(1),zhi(1)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_)

               gpx(i,j,domlo(3)-1) = gpx(i,j,domlo(3))
               gpy(i,j,domlo(3)-1) = gpy(i,j,domlo(3))
               gpz(i,j,domlo(3)-1) = gpz(i,j,domlo(3))

            case ( minf_)

               if (nodal_pressure .eq. 1) then
                  gpz(i,j,domlo(3)-1) = gpz(i,j,domlo(3))
               else
                  gpz(i,j,domlo(3)  ) = gpz(i,j,domlo(3)+1)
               end if
   
               gpx(i,j,domlo(3)-1) = 0.0d0
               gpy(i,j,domlo(3)-1) = 0.0d0
   
            case ( nsw_)

               if (nodal_pressure .eq. 1) then
                  gpz(i,j,domlo(3)-1) = -gpz(i,j,domlo(3))
               else
                  gpz(i,j,domlo(3)  ) =  0.0d0
               end if
   
               gpx(i,j,domlo(3)-1) = -gpx(i,j,domlo(3))
               gpy(i,j,domlo(3)-1) = -gpy(i,j,domlo(3))

            case ( fsw_)

               if (nodal_pressure .eq. 1) then
                  gpz(i,j,domlo(3)-1) = -gpz(i,j,domlo(3))
               else
                  gpz(i,j,domlo(3)  ) = 0.0d0
               end if
   
               gpx(i,j,domlo(3)-1) = gpx(i,j,domlo(3))
               gpy(i,j,domlo(3)-1) = gpy(i,j,domlo(3))

            end select
      end do
      end do
   end if

   if (zhi(3).ge.domhi(3)+1) then

      do j = zlo(2),zhi(2)
      do i = zlo(1),zhi(1)

            select case (bct_khi(i,j,1))

            case ( pinf_, pout_)

               gpx(i,j,domhi(3)+1) = gpz(i,j,domhi(3))
               gpy(i,j,domhi(3)+1) = gpy(i,j,domhi(3))
               gpz(i,j,domhi(3)+1) = gpz(i,j,domhi(3))
   
            case ( minf_)

               gpz(i,j,domhi(3)+1) = gpz(i,j,domhi(3))
   
               gpx(i,j,domhi(3)+1) = 0.0d0
               gpy(i,j,domhi(3)+1) = 0.0d0
   
            case ( nsw_)

               if (nodal_pressure .eq. 1) then
                  gpz(i,j,domhi(3)+1) = -gpz(i,j,domhi(3))
               else
                  gpz(i,j,domhi(3)+1) =  0.0d0
               end if

               gpx(i,j,domhi(3)+1) = -gpx(i,j,domhi(3))
               gpy(i,j,domhi(3)+1) = -gpy(i,j,domhi(3))
   
            case ( fsw_)
   
               if (nodal_pressure .eq. 1) then
                  gpz(i,j,domhi(3)+1) = -gpz(i,j,domhi(3))
               else
                  gpz(i,j,domhi(3)+1) = 0.0d0
               end if

               gpx(i,j,domhi(3)+1) = gpx(i,j,domhi(3))
               gpy(i,j,domhi(3)+1) = gpy(i,j,domhi(3))

            end select
      end do
      end do
   end if

end subroutine set_gradp_bcs
