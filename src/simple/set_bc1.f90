!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1                                                 C
!  Purpose: This subroutine sets boundary conditions for pressure      C
!  at mass inflow and velocity (normal and tangential at all walls).   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

subroutine set_bc1(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
     u_g, v_g, w_g, p_g, ep_g, ro_g, rop_g, mu_g, lambda_g, &
     bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
     bc_klo_type, bc_khi_type, domlo, domhi, ng ) &
     bind(C, name="set_bc1")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use bc, only: pinf_, pout_, minf_, nsw_, fsw_, psw_
      use bc, only: bc_u_g, bc_v_g, bc_w_g, bc_ep_g, bc_t_g
      use bc, only: bc_uw_g, bc_vw_g, bc_ww_g

      use eos      , only: eosg, sutherland
      use fld_const, only: ro_g0, mu_g0
      use param    , only: is_undefined

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3), ng

      real(c_real), intent(inout) ::  p_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ro_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: mu_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: lambda_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) ::  u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) ::  v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) ::  w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))


      integer(c_int), intent(in   ) :: &
           bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

! Local variables
!--------------------------------------------------------------------//
! local index for boundary condition
      integer :: bcv, i,j,k
      integer :: nlft, nrgt, nbot, ntop, nup, ndwn
      real(c_real) :: bc_ro_g, bc_mu_g, bc_lambda_g
!--------------------------------------------------------------------//

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

      bc_ro_g = ro_g0

      if (nlft .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               bcv = bc_ilo_type(j,k,2)
               if(bc_ilo_type(j,k,1) == PINF_ .or. &
                  bc_ilo_type(j,k,1) == POUT_) then

                       u_g(ulo(1):domlo(1)-1,j,k) =      u_g(domlo(1),j,k)
                       v_g(vlo(1):domlo(1)-1,j,k) =      v_g(domlo(1),j,k)
                       w_g(wlo(1):domlo(1)-1,j,k) =      w_g(domlo(1),j,k)

                      ep_g(slo(1):domlo(1)-1,j,k) =     ep_g(domlo(1),j,k)
                      ro_g(slo(1):domlo(1)-1,j,k) =     ro_g(domlo(1),j,k)
                     rop_g(slo(1):domlo(1)-1,j,k) =    rop_g(domlo(1),j,k)
                      mu_g(slo(1):domlo(1)-1,j,k) =     mu_g(domlo(1),j,k)
                  lambda_g(slo(1):domlo(1)-1,j,k) = lambda_g(domlo(1),j,k)

               else if (bc_ilo_type(j,k,1) == MINF_) then

                  p_g(slo(1):domlo(1)-1,j,k) = &
                       2*p_g(domlo(1),j,k) - p_g(domlo(1)+1,j,k)

                       u_g(ulo(1):domlo(1)  ,j,k) = bc_u_g(bcv)
                       v_g(vlo(1):domlo(1)-1,j,k) = 0.0d0
                       w_g(wlo(1):domlo(1)-1,j,k) = 0.0d0

                       if (is_undefined(mu_g0)) then
                          bc_mu_g     = sutherland(bc_t_g(bcv))
                          bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                       else
                          bc_mu_g     = mu_g0
                          bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                       endif

                      ep_g(slo(1):domlo(1)-1,j,k) = bc_ep_g(bcv)   
                      ro_g(slo(1):domlo(1)-1,j,k) = bc_ro_g
                     rop_g(slo(1):domlo(1)-1,j,k) = bc_ro_g * bc_ep_g(bcv)
                      mu_g(slo(1):domlo(1)-1,j,k) = bc_mu_g
                  lambda_g(slo(1):domlo(1)-1,j,k) = bc_lambda_g

               else if (bc_ilo_type(j,k,1) == NSW_) then

                  u_g(ulo(1):domlo(1)  ,j,k) =  0.0d0
                  v_g(vlo(1):domlo(1)-1,j,k) = -v_g(domlo(1),j,k)
                  w_g(wlo(1):domlo(1)-1,j,k) = -w_g(domlo(1),j,k)

               else if (bc_ilo_type(j,k,1) == FSW_) then

                  u_g(ulo(1):domlo(1)  ,j,k) = 0.0d0
                  v_g(vlo(1):domlo(1)-1,j,k) = v_g(domlo(1),j,k)
                  w_g(wlo(1):domlo(1)-1,j,k) = w_g(domlo(1),j,k)

               end if
            end do
         end do
      endif

      if (nrgt .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               bcv = bc_ihi_type(j,k,2)
               if(bc_ihi_type(j,k,1) == PINF_ .or. &
                  bc_ihi_type(j,k,1) == POUT_) then

                       u_g(domhi(1)+2:uhi(1),j,k) =      u_g(domhi(1)+1,j,k)
                       v_g(domhi(1)+1:vhi(1),j,k) =      v_g(domhi(1)  ,j,k)
                       w_g(domhi(1)+1:whi(1),j,k) =      w_g(domhi(1)  ,j,k)

                      ep_g(domhi(1)+1:shi(1),j,k) =     ep_g(domhi(1)  ,j,k)
                      ro_g(domhi(1)+1:shi(1),j,k) =     ro_g(domhi(1)  ,j,k)
                     rop_g(domhi(1)+1:shi(1),j,k) =    rop_g(domhi(1)  ,j,k)
                      mu_g(domhi(1)+1:shi(1),j,k) =     mu_g(domhi(1)  ,j,k)
                  lambda_g(domhi(1)+1:shi(1),j,k) = lambda_g(domhi(1)  ,j,k)

               else if (bc_ihi_type(j,k,1) == MINF_) then

                  p_g(domhi(1)+1:shi(1),j,k) = &
                       2*p_g(domhi(1),j,k) - p_g(domhi(1)-1,j,k)

                       u_g(domhi(1)+1:uhi(1),j,k) = bc_u_g(bcv)
                       v_g(domhi(1)+1:vhi(1),j,k) = 0.0d0
                       w_g(domhi(1)+1:whi(1),j,k) = 0.0d0

                       if (is_undefined(mu_g0)) then
                          bc_mu_g     = sutherland(bc_t_g(bcv))
                          bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                       else
                          bc_mu_g     = mu_g0
                          bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                       endif

                      ep_g(domhi(1)+1:shi(1),j,k) = bc_ep_g(bcv)   
                      ro_g(domhi(1)+1:shi(1),j,k) = bc_ro_g
                     rop_g(domhi(1)+1:shi(1),j,k) = bc_ro_g * bc_ep_g(bcv)   
                      mu_g(domhi(1)+1:shi(1),j,k) = bc_mu_g
                  lambda_g(domhi(1)+1:shi(1),j,k) = bc_lambda_g

               else if (bc_ihi_type(j,k,1) == NSW_) then

                  u_g(domhi(1)+1:uhi(1),j,k) =  0.0d0
                  v_g(domhi(1)+1:vhi(1),j,k) = -v_g(domhi(1),j,k)
                  w_g(domhi(1)+1:whi(1),j,k) = -w_g(domhi(1),j,k)

               else if (bc_ihi_type(j,k,1) == FSW_) then

                  u_g(domhi(1)+1:uhi(1),j,k) = 0.0d0
                  v_g(domhi(1)+1:vhi(1),j,k) = v_g(domhi(1),j,k)
                  w_g(domhi(1)+1:whi(1),j,k) = w_g(domhi(1),j,k)

               end if

            end do
         end do
      endif

      if (nbot .gt. 0) then
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               bcv = bc_jlo_type(i,k,2)
               if(bc_jlo_type(i,k,1) == PINF_ .or. &
                    bc_jlo_type(i,k,1) == POUT_) then

                       u_g(i,ulo(2):domlo(2)-1,k) =      u_g(i,domlo(2),k)
                       v_g(i,vlo(2):domlo(2)-1,k) =      v_g(i,domlo(2),k)
                       w_g(i,wlo(2):domlo(2)-1,k) =      w_g(i,domlo(2),k)

                      ep_g(i,slo(1):domlo(2)-1,k) =     ep_g(i,domlo(2),k)
                      ro_g(i,slo(1):domlo(2)-1,k) =     ro_g(i,domlo(2),k)
                     rop_g(i,slo(1):domlo(2)-1,k) =    rop_g(i,domlo(2),k)
                      mu_g(i,slo(1):domlo(2)-1,k) =     mu_g(i,domlo(2),k)
                  lambda_g(i,slo(1):domlo(2)-1,k) = lambda_g(i,domlo(2),k)

               else if (bc_jlo_type(i,k,1) == MINF_)then

                  p_g(i,slo(2):domlo(2)-1,k) = &
                       2*p_g(i,domlo(2),k) - p_g(i,domlo(2)+1,k)

                       u_g(i,ulo(2):domlo(2)-1,k) = 0.0d0
                       v_g(i,vlo(2):domlo(2)  ,k) = bc_v_g(bcv)
                       w_g(i,wlo(2):domlo(2)-1,k) = 0.0d0

                       if (is_undefined(mu_g0)) then
                          bc_mu_g     = sutherland(bc_t_g(bcv))
                          bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                       else
                          bc_mu_g     = mu_g0
                          bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                       endif

                      ep_g(i,slo(1):domlo(2)-1,k) = bc_ep_g(bcv)   
                      ro_g(i,slo(1):domlo(2)-1,k) = bc_ro_g
                     rop_g(i,slo(1):domlo(2)-1,k) = bc_ro_g * bc_ep_g(bcv)   
                      mu_g(i,slo(1):domlo(2)-1,k) = bc_mu_g
                  lambda_g(i,slo(1):domlo(2)-1,k) = bc_lambda_g

               else if (bc_jlo_type(i,k,1) == NSW_)then

                  u_g(i,ulo(2):domlo(2)-1,k) = -u_g(i,domlo(2),k)
                  v_g(i,vlo(2):domlo(2)  ,k) =  0.0d0
                  w_g(i,wlo(2):domlo(2)-1,k) = -w_g(i,domlo(2),k)

               else if (bc_jlo_type(i,k,1) == FSW_)then

                  u_g(i,ulo(2):domlo(2)-1,k) = u_g(i,domlo(2),k)
                  v_g(i,vlo(2):domlo(2)  ,k) = 0.0d0
                  w_g(i,wlo(2):domlo(2)-1,k) = w_g(i,domlo(2),k)

               end if

            end do
         end do
      endif

      if (ntop .gt. 0) then
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               bcv = bc_jhi_type(i,k,2)
               if(bc_jhi_type(i,k,1) == PINF_ .or. &
                  bc_jhi_type(i,k,1) == POUT_) then

                       u_g(i,domhi(2)+1:uhi(2),k) =      u_g(i,domhi(2)  ,k)
                       v_g(i,domhi(2)+2:vhi(2),k) =      v_g(i,domhi(2)+1,k)
                       w_g(i,domhi(2)+1:whi(2),k) =      w_g(i,domhi(2)  ,k)

                      ep_g(i,domhi(2)+1:shi(2),k) =     ep_g(i,domhi(2)  ,k)
                      ro_g(i,domhi(2)+1:shi(2),k) =     ro_g(i,domhi(2)  ,k)
                     rop_g(i,domhi(2)+1:shi(2),k) =    rop_g(i,domhi(2)  ,k)
                      mu_g(i,domhi(2)+1:shi(2),k) =     mu_g(i,domhi(2)  ,k)
                  lambda_g(i,domhi(2)+1:shi(2),k) = lambda_g(i,domhi(2)  ,k)

               else if (bc_jhi_type(i,k,1) == MINF_) then

                  p_g(i,domhi(2)+1:shi(2),k) = &
                       2*p_g(i,domhi(2),k) - p_g(i,domhi(2)-1,k)

                       u_g(i,domhi(2)+1:uhi(2),k) = 0.0d0
                       v_g(i,domhi(2)+1:vhi(2),k) = bc_v_g(bcv)
                       w_g(i,domhi(2)+1:whi(2),k) = 0.0d0

                       if (is_undefined(mu_g0)) then
                          bc_mu_g     = sutherland(bc_t_g(bcv))
                          bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                       else
                          bc_mu_g     = mu_g0
                          bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                       endif

                      ep_g(i,domhi(2)+1:shi(2),k) = bc_ep_g(bcv)   
                      ro_g(i,domhi(2)+1:shi(2),k) = bc_ro_g
                     rop_g(i,domhi(2)+1:shi(2),k) = bc_ro_g * bc_ep_g(bcv)   
                      mu_g(i,domhi(2)+1:shi(2),k) = bc_mu_g
                  lambda_g(i,domhi(2)+1:shi(2),k) = bc_lambda_g

               else if (bc_jhi_type(i,k,1) == NSW_) then

                  u_g(i,domhi(2)+1:uhi(2),k) = -u_g(i,domhi(2),k)
                  v_g(i,domhi(2)+1:vhi(2),k) =  0.0d0
                  w_g(i,domhi(2)+1:whi(2),k) = -w_g(i,domhi(2),k)

               else if (bc_jhi_type(i,k,1) == FSW_) then

                  u_g(i,domhi(2)+1:uhi(2),k) = u_g(i,domhi(2),k)
                  v_g(i,domhi(2)+1:vhi(2),k) = 0.0d0
                  w_g(i,domhi(2)+1:whi(2),k) = w_g(i,domhi(2),k)

               end if
            end do
         end do
      endif

      if (ndwn .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               bcv = bc_klo_type(i,j,2)
               if(bc_klo_type(i,j,1) == PINF_ .or. &
                  bc_klo_type(i,j,1) == POUT_) then

                       u_g(i,j,ulo(3):domlo(3)-1) =      u_g(i,j,domlo(3))
                       v_g(i,j,vlo(3):domlo(3)-1) =      v_g(i,j,domlo(3))
                       w_g(i,j,wlo(3):domlo(3)-1) =      w_g(i,j,domlo(3))

                      ep_g(i,j,slo(3):domlo(3)-1) =     ep_g(i,j,domlo(3))
                      ro_g(i,j,slo(3):domlo(3)-1) =     ro_g(i,j,domlo(3))
                     rop_g(i,j,slo(3):domlo(3)-1) =    rop_g(i,j,domlo(3))
                      mu_g(i,j,slo(3):domlo(3)-1) =     mu_g(i,j,domlo(3))
                  lambda_g(i,j,slo(3):domlo(3)-1) = lambda_g(i,j,domlo(3))

               else if (bc_klo_type(i,j,1) == MINF_) then

                  p_g(i,j,slo(3):domlo(3)-1) = &
                       2*p_g(i,j,domlo(3)) - p_g(i,j,domlo(3)+1)

                       u_g(i,j,ulo(3):domlo(3)-1) = 0.0d0
                       v_g(i,j,vlo(3):domlo(3)-1) = 0.0d0
                       w_g(i,j,wlo(3):domlo(3)  ) = bc_w_g(bcv)

                       if (is_undefined(mu_g0)) then
                          bc_mu_g     = sutherland(bc_t_g(bcv))
                          bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                       else
                          bc_mu_g     = mu_g0
                          bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                       endif

                      ep_g(i,j,slo(3):domlo(3)-1) = bc_ep_g(bcv)   
                      ro_g(i,j,slo(3):domlo(3)-1) = bc_ro_g
                     rop_g(i,j,slo(3):domlo(3)-1) = bc_ro_g * bc_ep_g(bcv)   
                      mu_g(i,j,slo(3):domlo(3)-1) = bc_mu_g
                  lambda_g(i,j,slo(3):domlo(3)-1) = bc_lambda_g

               else if (bc_klo_type(i,j,1) == NSW_) then

                  u_g(i,j,ulo(3):domlo(3)-1) = -u_g(i,j,domlo(3))
                  v_g(i,j,vlo(3):domlo(3)-1) = -v_g(i,j,domlo(3))
                  w_g(i,j,wlo(3):domlo(3)  ) =  0.0d0

               else if (bc_klo_type(i,j,1) == FSW_) then

                  u_g(i,j,ulo(3):domlo(3)-1) = u_g(i,j,domlo(3))
                  v_g(i,j,vlo(3):domlo(3)-1) = v_g(i,j,domlo(3))
                  w_g(i,j,wlo(3):domlo(3)  ) = 0.0d0
               end if
            end do
         end do
      endif

      if (nup .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               bcv = bc_khi_type(i,j,2)
               if(bc_khi_type(i,j,1) == PINF_ .or. &
                  bc_khi_type(i,j,1) == POUT_) then

                       u_g(i,j,domhi(3)+1:uhi(3)) =      u_g(i,j,domhi(3)  )
                       v_g(i,j,domhi(3)+1:vhi(3)) =      v_g(i,j,domhi(3)  )
                       w_g(i,j,domhi(3)+2:whi(3)) =      w_g(i,j,domhi(3)+1)

                      ep_g(i,j,domhi(3)+1:shi(3)) =     ep_g(i,j,domhi(3)  )
                      ro_g(i,j,domhi(3)+1:shi(3)) =     ro_g(i,j,domhi(3)  )
                     rop_g(i,j,domhi(3)+1:shi(3)) =    rop_g(i,j,domhi(3)  )
                      mu_g(i,j,domhi(3)+1:shi(3)) =     mu_g(i,j,domhi(3)  )
                  lambda_g(i,j,domhi(3)+1:shi(3)) = lambda_g(i,j,domhi(3)  )

               else if (bc_khi_type(i,j,1) == MINF_) then

                  p_g(i,j,domhi(3)+1:shi(3)) = &
                       2*p_g(i,j,domhi(2)) - p_g(i,j,domhi(3)-1)

                       u_g(i,j,domhi(3)+1:uhi(3)) = 0.0d0
                       v_g(i,j,domhi(3)+1:vhi(3)) = 0.0d0
                       w_g(i,j,domhi(3)+1:whi(3)) = bc_w_g(bcv)

                       if (is_undefined(mu_g0)) then
                          bc_mu_g     = sutherland(bc_t_g(bcv))
                          bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                       else
                          bc_mu_g     = mu_g0
                          bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                       endif

                      ep_g(i,j,domhi(3)+1:shi(3)) = bc_ep_g(bcv)   
                      ro_g(i,j,domhi(3)+1:shi(3)) = bc_ro_g
                     rop_g(i,j,domhi(3)+1:shi(3)) = bc_ro_g * bc_ep_g(bcv)   
                      mu_g(i,j,domhi(3)+1:shi(3)) = bc_mu_g
                  lambda_g(i,j,domhi(3)+1:shi(3)) = bc_lambda_g

               else if (bc_khi_type(i,j,1) == NSW_) then

                  u_g(i,j,domhi(3)+1:uhi(3)) = -u_g(i,j,domhi(3))
                  v_g(i,j,domhi(3)+1:vhi(3)) = -v_g(i,j,domhi(3))
                  w_g(i,j,domhi(3)+1:whi(3)) =  0.0d0

               else if (bc_khi_type(i,j,1) == FSW_) then

                  u_g(i,j,domhi(3)+1:uhi(3)) = u_g(i,j,domhi(3))
                  v_g(i,j,domhi(3)+1:vhi(3)) = v_g(i,j,domhi(3))
                  w_g(i,j,domhi(3)+1:whi(3)) = 0.0d0

               end if
            end do
         end do
      endif

   ! *********************************************************************************
   ! We have to do the PSW bc's last because otherwise non-zero moving wall values
   ! can get over-written
   ! *********************************************************************************

      if (nlft .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               bcv = bc_ilo_type(j,k,2)
               if (bc_ilo_type(j,k,1) == PSW_) then
                  u_g(ulo(1):domlo(1)  ,j,k) = 0.0d0
                  v_g(vlo(1):domlo(1)-1,j,k) = 2.0*bc_vw_g(bcv) - v_g(domlo(1),j,k)
                  w_g(wlo(1):domlo(1)-1,j,k) = 2.0*bc_ww_g(bcv) - w_g(domlo(1),j,k)
               end if
            end do
         end do
      endif

      if (nrgt .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               bcv = bc_ihi_type(j,k,2)
               if (bc_ihi_type(j,k,1) == PSW_) then
                  u_g(domhi(1)+1:uhi(1),j,k) = 0.0d0
                  v_g(domhi(1)+1:vhi(1),j,k) = 2.0*bc_vw_g(bcv) - v_g(domhi(1),j,k)
                  w_g(domhi(1)+1:whi(1),j,k) = 2.0*bc_ww_g(bcv) - w_g(domhi(1),j,k)
               end if
            end do
         end do
      endif

      if (nbot .gt. 0) then
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               bcv = bc_jlo_type(i,k,2)
               if (bc_jlo_type(i,k,1) == PSW_)then
                  u_g(i,ulo(2):domlo(2)-1,k) = 2.0*bc_uw_g(bcv) - u_g(i,domlo(2),k)
                  v_g(i,vlo(2):domlo(2)  ,k) = 0.0d0
                  w_g(i,wlo(2):domlo(2)-1,k) = 2.0*bc_ww_g(bcv) - w_g(i,domlo(2),k)
               end if
            end do
         end do
      endif

      if (ntop .gt. 0) then
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               bcv = bc_jhi_type(i,k,2)
               if (bc_jhi_type(i,k,1) == PSW_)then
                  u_g(i,domhi(2)+1:uhi(2),k) = 2.0*bc_uw_g(bcv) - u_g(i,domhi(2),k)
                  v_g(i,domhi(2)+1:vhi(2),k) = 0.0d0
                  w_g(i,domhi(2)+1:whi(2),k) = 2.0*bc_ww_g(bcv) - w_g(i,domhi(2),k)
               end if
            end do
         end do
      endif

      if (ndwn .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               bcv = bc_klo_type(i,j,2)
               if (bc_klo_type(i,j,1) == PSW_) then
                  u_g(i,j,ulo(3):domlo(3)-1) = 2.0*bc_uw_g(bcv) - u_g(i,j,domlo(3))
                  v_g(i,j,vlo(3):domlo(3)-1) = 2.0*bc_vw_g(bcv) - v_g(i,j,domlo(3))
                  w_g(i,j,wlo(3):domlo(3)  ) = 0.0d0
               end if
            end do
         end do
      endif

      if (nup .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               bcv = bc_khi_type(i,j,2)
               if (bc_khi_type(i,j,1) == PSW_) then
                  u_g(i,j,domhi(3)+1:uhi(3)) = 2.0*bc_uw_g(bcv) - u_g(i,j,domhi(3))
                  v_g(i,j,domhi(3)+1:vhi(3)) = 2.0*bc_vw_g(bcv) - v_g(i,j,domhi(3))
                  w_g(i,j,domhi(3)+1:whi(3)) = 0.0d0
               end if
            end do
         end do
      endif

   end subroutine set_bc1
