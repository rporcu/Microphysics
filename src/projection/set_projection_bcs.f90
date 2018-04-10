! 
!              
!  This subroutine sets the BCs for all variables involved in
!  the projection, EXCEPT the pressure and velocity components. 
!  To set pressure BCs, use "set_pressure_bcs".
!  To set velocity BCs, use "set_velocity_bcs".
!  
!  Author: Michele Rosso
! 
!  Date: December 20, 2017
!
! 
subroutine set_projection_bcs ( ep_g, slo, shi, ro_g, rop_g, mu_g, lambda_g,&
     & bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi,           &
     & domlo, domhi, ng ) bind(C) 

   use amrex_fort_module,  only: ar => amrex_real
   use iso_c_binding ,     only: c_int
   use bc
   use scales,             only: scale_pressure
   use eos      ,          only: eosg, sutherland
   use fld_const,          only: ro_g0, mw_avg, mu_g0
   use param    ,          only: is_undefined
   
   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)

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
        ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),     &
        ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),     &
        rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),    &
        mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),     &
        lambda_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   ! Local variables
   integer  :: bcv, i, j, k
   integer  :: nlft, nrgt, nbot, ntop, nup, ndwn
   real(ar) :: bc_ro_g, bc_mu_g, bc_lambda_g

   nlft = max(0,domlo(1)-slo(1))
   nbot = max(0,domlo(2)-slo(2))
   ndwn = max(0,domlo(3)-slo(3))

   nrgt = max(0,shi(1)-domhi(1))
   ntop = max(0,shi(2)-domhi(2))
   nup  = max(0,shi(3)-domhi(3))

   bc_ro_g = ro_g0
   
   if (nlft .gt. 0) then
      do k = slo(3), shi(3)
         do j = slo(2), shi(2)

            bcv = bct_ilo(j,k,2)

            select case (bct_ilo(j,k,1))
               
            case ( pinf_, pout_) 

               ep_g(slo(1):domlo(1)-1,j,k) =     ep_g(domlo(1),j,k)
               ro_g(slo(1):domlo(1)-1,j,k) =     ro_g(domlo(1),j,k)
               rop_g(slo(1):domlo(1)-1,j,k) =    rop_g(domlo(1),j,k)
               mu_g(slo(1):domlo(1)-1,j,k) =     mu_g(domlo(1),j,k)
               lambda_g(slo(1):domlo(1)-1,j,k) = lambda_g(domlo(1),j,k)

            case ( minf_)

               if (is_undefined(mu_g0)) then
                  bc_mu_g     = sutherland(bc_t_g(bcv))
                  bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
               else
                  bc_mu_g     = mu_g0
                  bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
               endif

               ! ep_g and rop_g BCs are inconsistent. This is because
               ! projection implementation requires bc_ep_g(bcv) to be defined
               ! at the inlet face. Therefore, ep_g at the first ghost cell will
               ! have a fictitious value. Therefore we cannot set rop_g to be
               ! ro_g * ep_g at the ghost node, otherwise even ro_g would have a
               ! fictitious value. However, this create an inconsistency
               ep_g(domlo(1)-1,j,k)  = 2.0_ar * bc_ep_g(bcv) - ep_g(domlo(1),j,k)
               
               ro_g(slo(1):domlo(1)-1,j,k)     = bc_ro_g
               rop_g(slo(1):domlo(1)-1,j,k)    = bc_ro_g * bc_ep_g(bcv)
               mu_g(slo(1):domlo(1)-1,j,k)     = bc_mu_g
               lambda_g(slo(1):domlo(1)-1,j,k) = bc_lambda_g

            end select
            
         end do
      end do
   endif

   if (nrgt .gt. 0) then
      
      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
            
            bcv = bct_ihi(j,k,2)
            
            select case ( bct_ihi(j,k,1) )

            case ( pinf_, pout_ )

               ep_g(domhi(1)+1:shi(1),j,k) =     ep_g(domhi(1)  ,j,k)
               ro_g(domhi(1)+1:shi(1),j,k) =     ro_g(domhi(1)  ,j,k)
               rop_g(domhi(1)+1:shi(1),j,k) =    rop_g(domhi(1)  ,j,k)
               mu_g(domhi(1)+1:shi(1),j,k) =     mu_g(domhi(1)  ,j,k)
               lambda_g(domhi(1)+1:shi(1),j,k) = lambda_g(domhi(1)  ,j,k)

            case ( minf_ )

               if (is_undefined(mu_g0)) then
                  bc_mu_g     = sutherland(bc_t_g(bcv))
                  bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
               else
                  bc_mu_g     = mu_g0
                  bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
               endif
               
               ep_g(domhi(1)+1,j,k) = 2.0_ar*bc_ep_g(bcv) - ep_g(domhi(1),j,k)   
               ro_g(domhi(1)+1:shi(1),j,k) = bc_ro_g
               rop_g(domhi(1)+1:shi(1),j,k) = bc_ro_g * bc_ep_g(bcv)   
               mu_g(domhi(1)+1:shi(1),j,k) = bc_mu_g
               lambda_g(domhi(1)+1:shi(1),j,k) = bc_lambda_g

            end select

         end do
      end do
   endif

   if (nbot .gt. 0) then
      
      do k = slo(3), shi(3)
         do i = slo(1), shi(1)
            
            bcv = bct_jlo(i,k,2)

            select case ( bct_jlo(i,k,1) )

            case ( pinf_, pout_) 

               ep_g(i,slo(1):domlo(2)-1,k) =     ep_g(i,domlo(2),k)
               ro_g(i,slo(1):domlo(2)-1,k) =     ro_g(i,domlo(2),k)
               rop_g(i,slo(1):domlo(2)-1,k) =    rop_g(i,domlo(2),k)
               mu_g(i,slo(1):domlo(2)-1,k) =     mu_g(i,domlo(2),k)
               lambda_g(i,slo(1):domlo(2)-1,k) = lambda_g(i,domlo(2),k)

            case ( minf_ )

               if (is_undefined(mu_g0)) then
                  bc_mu_g     = sutherland(bc_t_g(bcv))
                  bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
               else
                  bc_mu_g     = mu_g0
                  bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
               endif

               ep_g(i,domlo(2)-1,k) = 2.0_ar*bc_ep_g(bcv) - ep_g(i,domlo(2),k)   
               ro_g(i,slo(1):domlo(2)-1,k) = bc_ro_g
               rop_g(i,slo(1):domlo(2)-1,k) = bc_ro_g * bc_ep_g(bcv)   
               mu_g(i,slo(1):domlo(2)-1,k) = bc_mu_g
               lambda_g(i,slo(1):domlo(2)-1,k) = bc_lambda_g

            end select

         end do
      end do
   endif

   if (ntop .gt. 0) then

      do k = slo(3), shi(3)
         do i = slo(1), shi(1)
            
            bcv = bct_jhi(i,k,2)

            select case ( bct_jhi(i,k,1) )

            case ( pinf_, pout_ )

               ep_g(i,domhi(2)+1:shi(2),k) =     ep_g(i,domhi(2)  ,k)
               ro_g(i,domhi(2)+1:shi(2),k) =     ro_g(i,domhi(2)  ,k)
               rop_g(i,domhi(2)+1:shi(2),k) =    rop_g(i,domhi(2)  ,k)
               mu_g(i,domhi(2)+1:shi(2),k) =     mu_g(i,domhi(2)  ,k)
               lambda_g(i,domhi(2)+1:shi(2),k) = lambda_g(i,domhi(2)  ,k)

            case ( minf_) 

               if (is_undefined(mu_g0)) then
                  bc_mu_g     = sutherland(bc_t_g(bcv))
                  bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
               else
                  bc_mu_g     = mu_g0
                  bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
               endif

               ep_g(i,domhi(2)+1,k) = 2.0_ar*bc_ep_g(bcv) - ep_g(i,domhi(2),k)   
               ro_g(i,domhi(2)+1:shi(2),k) = bc_ro_g
               rop_g(i,domhi(2)+1:shi(2),k) = bc_ro_g * bc_ep_g(bcv)   
               mu_g(i,domhi(2)+1:shi(2),k) = bc_mu_g
               lambda_g(i,domhi(2)+1:shi(2),k) = bc_lambda_g

            end select
         end do
      end do
   endif

   if (ndwn .gt. 0) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)
            
            bcv = bct_klo(i,j,2)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_ ) 

               ep_g(i,j,slo(3):domlo(3)-1) =     ep_g(i,j,domlo(3))
               ro_g(i,j,slo(3):domlo(3)-1) =     ro_g(i,j,domlo(3))
               rop_g(i,j,slo(3):domlo(3)-1) =    rop_g(i,j,domlo(3))
               mu_g(i,j,slo(3):domlo(3)-1) =     mu_g(i,j,domlo(3))
               lambda_g(i,j,slo(3):domlo(3)-1) = lambda_g(i,j,domlo(3))

            case ( minf_ )

               if (is_undefined(mu_g0)) then
                  bc_mu_g     = sutherland(bc_t_g(bcv))
                  bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
               else
                  bc_mu_g     = mu_g0
                  bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
               endif
               
               ep_g(i,j,domlo(3)-1) = 2.0_ar*bc_ep_g(bcv) - ep_g(i,j,domlo(3))   
               ro_g(i,j,slo(3):domlo(3)-1) = bc_ro_g
               rop_g(i,j,slo(3):domlo(3)-1) = bc_ro_g * bc_ep_g(bcv)   
               mu_g(i,j,slo(3):domlo(3)-1) = bc_mu_g
               lambda_g(i,j,slo(3):domlo(3)-1) = bc_lambda_g

            end select
         end do
      end do
   endif

   if (nup .gt. 0) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)
            
            bcv = bct_khi(i,j,2)

            select case ( bct_khi(i,j,1) )

            case ( pinf_, pout_ )
             
               ep_g(i,j,domhi(3)+1:shi(3)) =     ep_g(i,j,domhi(3)  )
               ro_g(i,j,domhi(3)+1:shi(3)) =     ro_g(i,j,domhi(3)  )
               rop_g(i,j,domhi(3)+1:shi(3)) =    rop_g(i,j,domhi(3)  )
               mu_g(i,j,domhi(3)+1:shi(3)) =     mu_g(i,j,domhi(3)  )
               lambda_g(i,j,domhi(3)+1:shi(3)) = lambda_g(i,j,domhi(3)  )

            case ( minf_ ) 

               if (is_undefined(mu_g0)) then
                  bc_mu_g     = sutherland(bc_t_g(bcv))
                  bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
               else
                  bc_mu_g     = mu_g0
                  bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
               endif
               
               ep_g(i,j,domhi(3)+1) = 2.0_ar*bc_ep_g(bcv) - ep_g(i,j,domhi(3)+1)   
               ro_g(i,j,domhi(3)+1:shi(3)) = bc_ro_g
               rop_g(i,j,domhi(3)+1:shi(3)) = bc_ro_g * bc_ep_g(bcv)   
               mu_g(i,j,domhi(3)+1:shi(3)) = bc_mu_g
               lambda_g(i,j,domhi(3)+1:shi(3)) = bc_lambda_g

            end select
         end do
      end do
   endif

end subroutine set_projection_bcs
