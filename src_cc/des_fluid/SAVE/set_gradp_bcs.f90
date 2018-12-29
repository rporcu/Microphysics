subroutine set_gradp_bcs ( slo, shi, gp, glo, ghi, &
     & bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi,               &
     & domlo, domhi, ng) bind(C) 

   use amrex_fort_module,  only: ar => amrex_real
   use param            ,  only: zero
   use iso_c_binding ,     only: c_int
   use bc

   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)
   integer(c_int), intent(in   ) :: glo(3), ghi(3)

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
        gp(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3)

   ! Local variables
   integer  :: i, j, k

   if (glo(1).le.domlo(1)) then
      do k = glo(3),ghi(3)
      do j = glo(2),ghi(2)

            select case (bct_ilo(j,k,1))

            case ( pinf_, pout_)

               gp(domlo(1)-1,j,k,:) = gp(domlo(1),j,k,:)

            case ( minf_)

               gp(domlo(1)-1,j,k,1) = gp(domlo(1),j,k,1  )
               gp(domlo(1)-1,j,k,2) = zero
               gp(domlo(1)-1,j,k,3) = zero

            case ( nsw_)

               gp(domlo(1)-1,j,k,:) = -gp(domlo(1),j,k,:)

            case ( fsw_)

               gp(domlo(1)-1,j,k,1) = -gp(domlo(1),j,k,1)
               gp(domlo(1)-1,j,k,2) =  gp(domlo(1),j,k,2)
               gp(domlo(1)-1,j,k,3) =  gp(domlo(1),j,k,3)

            case ( ignore_)

               gp(domlo(1)-1,j,k,1) =  gp(domlo(1),j,k,1)
               gp(domlo(1)-1,j,k,2) =  gp(domlo(1),j,k,2)
               gp(domlo(1)-1,j,k,3) =  gp(domlo(1),j,k,3)
   
            end select
      end do
      end do
   end if

   if (ghi(1).ge.domhi(1)+1) then
      do k = glo(3),ghi(3)
      do j = glo(2),ghi(2)

            select case (bct_ihi(j,k,1))

            case ( pinf_, pout_)

               gp(domhi(1)+1,j,k,:) = gp(domhi(1),j,k,:)

            case ( minf_)

               gp(domhi(1)+1,j,k,1) = gp(domhi(1),j,k,1  )
               gp(domhi(1)+1,j,k,2) = zero
               gp(domhi(1)+1,j,k,3) = zero

            case ( nsw_)

               gp(domhi(1)+1,j,k,1  ) = -gp(domhi(1),j,k,1  )

            case ( fsw_)

               gp(domhi(1)+1,j,k,1) = -gp(domhi(1),j,k,1)
               gp(domhi(1)+1,j,k,2) =  gp(domhi(1),j,k,2)
               gp(domhi(1)+1,j,k,3) =  gp(domhi(1),j,k,3)

            case ( ignore_)

               gp(domhi(1)+1,j,k,1) =  gp(domhi(1),j,k,1)
               gp(domhi(1)+1,j,k,2) =  gp(domhi(1),j,k,2)
               gp(domhi(1)+1,j,k,3) =  gp(domhi(1),j,k,3)

            end select
      end do
      end do
   end if

   if (glo(2).le.domlo(2)) then
      do k = glo(3),ghi(3)
      do i = glo(1),ghi(1)

            select case (bct_jlo(i,k,1))

            case ( pinf_, pout_)

               gp(i,domlo(2)-1,k,:) = gp(i,domlo(2),k,:)
   
            case ( minf_)

               gp(i,domlo(2)-1,k,2) = gp(i,domlo(2),k,1  )
               gp(i,domlo(2)-1,k,1) = zero
               gp(i,domlo(2)-1,k,3) = zero

            case ( nsw_)

               gp(i,domlo(2)-1,k,:) = -gp(i,domlo(2),k,:)

            case ( fsw_)

               gp(i,domlo(2)-1,k,2) = -gp(i,domlo(2),k,2)
               gp(i,domlo(2)-1,k,1) =  gp(i,domlo(2),k,1)
               gp(i,domlo(2)-1,k,3) =  gp(i,domlo(2),k,3)

            case ( ignore_)

               gp(i,domlo(2)-1,k,2) =  gp(i,domlo(2),k,2)
               gp(i,domlo(2)-1,k,1) =  gp(i,domlo(2),k,1)
               gp(i,domlo(2)-1,k,3) =  gp(i,domlo(2),k,3)

            end select
      end do
      end do
   end if

   if (ghi(2).ge.domhi(2)+1) then

      do k = glo(3),ghi(3)
      do i = glo(1),ghi(1)

            select case (bct_jhi(i,k,1))

            case ( pinf_, pout_)

               gp(i,domhi(2)+1,k,:) = gp(i,domhi(2),k,:)

            case ( minf_)

               gp(i,domhi(2)+1,k,2) = gp(i,domhi(2),k,2)
               gp(i,domhi(2)+1,k,1) = zero
               gp(i,domhi(2)+1,k,3) = zero

            case ( nsw_)

               gp(i,domhi(2)+1,k,:) = -gp(i,domhi(2),k,:)

            case ( fsw_)

               gp(i,domhi(2)+1,k,2) = -gp(i,domhi(2),k,2)
               gp(i,domhi(2)+1,k,1) =  gp(i,domhi(2),k,1)
               gp(i,domhi(2)+1,k,3) =  gp(i,domhi(2),k,3)

            case ( ignore_)

               gp(i,domhi(2)+1,k,2) =  gp(i,domhi(2),k,2)
               gp(i,domhi(2)+1,k,1) =  gp(i,domhi(2),k,1)
               gp(i,domhi(2)+1,k,3) =  gp(i,domhi(2),k,3)

            end select
      end do
      end do
   end if

   if (glo(3).le.domlo(3)) then

      do j = glo(2),ghi(2)
      do i = glo(1),ghi(1)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_)

               gp(i,j,domlo(3)-1,:) = gp(i,j,domlo(3),:)

            case ( minf_)

               gp(i,j,domlo(3)-1,3) = gp(i,j,domlo(3),3)
               gp(i,j,domlo(3)-1,1) = zero
               gp(i,j,domlo(3)-1,2) = zero

            case ( nsw_)

               gp(i,j,domlo(3)-1,:) = -gp(i,j,domlo(3),:)

            case ( fsw_)

               gp(i,j,domlo(3)-1,3) = -gp(i,j,domlo(3),3)
               gp(i,j,domlo(3)-1,1) =  gp(i,j,domlo(3),1)
               gp(i,j,domlo(3)-1,2) =  gp(i,j,domlo(3),2)

            case ( ignore_)

               gp(i,j,domlo(3)-1,3) =  gp(i,j,domlo(3),3)
               gp(i,j,domlo(3)-1,1) =  gp(i,j,domlo(3),1)
               gp(i,j,domlo(3)-1,2) =  gp(i,j,domlo(3),2)

            end select
      end do
      end do
   end if

   if (ghi(3).ge.domhi(3)+1) then

      do j = glo(2),ghi(2)
      do i = glo(1),ghi(1)

            select case (bct_khi(i,j,1))

            case ( pinf_, pout_)

               gp(i,j,domhi(3)+1,:) = gp(i,j,domhi(3),:)
   
            case ( minf_)

               gp(i,j,domhi(3)+1,3) = gp(i,j,domhi(3),3)
               gp(i,j,domhi(3)+1,1) = zero
               gp(i,j,domhi(3)+1,2) = zero
   
            case ( nsw_)

               gp(i,j,domhi(3)+1,:) = -gp(i,j,domhi(3),:)
   
            case ( fsw_)

               gp(i,j,domhi(3)+1,3) = -gp(i,j,domhi(3),3)
               gp(i,j,domhi(3)+1,1) =  gp(i,j,domhi(3),1)
               gp(i,j,domhi(3)+1,2) =  gp(i,j,domhi(3),2)

            case ( ignore_)

               gp(i,j,domhi(3)+1,3) =  gp(i,j,domhi(3),3)
               gp(i,j,domhi(3)+1,1) =  gp(i,j,domhi(3),1)
               gp(i,j,domhi(3)+1,2) =  gp(i,j,domhi(3),2)
   
            end select
      end do
      end do
   end if

end subroutine set_gradp_bcs
