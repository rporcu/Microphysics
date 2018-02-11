!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0                                                 C
!  Purpose: This subroutine does the initial setting of all boundary   C
!  conditions. The user specifications of the boundary conditions are  C
!  checked for veracity in various check_data/ routines:               C
!  (e.g., check_boundary_conditions).                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   subroutine set_bc0(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
                       u_g, v_g, w_g, p0_g, ep_g, ro_g, rop_g, mu_g, lambda_g, &
                       bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                       bc_klo_type, bc_khi_type, domlo, domhi) &
      bind(C, name="set_bc0")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use bc, only: bc_p_g, bc_ep_g, bc_t_g
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: pinf_, pout_, minf_

      use eos      , only: eosg, sutherland
      use fld_const, only: ro_g0, mw_avg, mu_g0

      use scales, only: scale_pressure
      use param , only: is_undefined

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(c_real), intent(inout) ::  u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) ::  v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) ::  w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(inout) ::  p0_g&
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

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)

! Local variables
!--------------------------------------------------------------------//
! local index for boundary condition
      integer :: bcv, i,j,k

      integer    nlft, nrgt, nbot, ntop, nup, ndwn

      real(c_real) :: bc_ro_g, bc_mu_g, bc_lambda_g
!--------------------------------------------------------------------//

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

      if (nlft .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)

               bcv = bc_ilo_type(j,k,2)

               if (bc_ilo_type(j,k,1) == PINF_ .or. &
                   bc_ilo_type(j,k,1) == POUT_ .or. &
                   bc_ilo_type(j,k,1) == MINF_) then

                  if (is_undefined(ro_g0)) then
                     bc_ro_g = eosg(mw_avg,bc_p_g(bcv),bc_t_g(bcv))
                  else
                     bc_ro_g = ro_g0
                  endif

                  if (is_undefined(mu_g0)) then
                     bc_mu_g     = sutherland(bc_t_g(bcv))
                     bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                  else
                     bc_mu_g     = mu_g0
                     bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                  endif

                      ep_g(slo(1):domlo(1)-1,j,k) = bc_ep_g(bcv)
                      ro_g(slo(1):domlo(1)-1,j,k) = bc_ro_g
                     rop_g(slo(1):domlo(1)-1,j,k) = bc_ro_g*bc_ep_g(bcv)
                      mu_g(slo(1):domlo(1)-1,j,k) = bc_mu_g
                  lambda_g(slo(1):domlo(1)-1,j,k) = bc_lambda_g

               end if

               if (bc_ilo_type(j,k,1) == PINF_ .or. &
                   bc_ilo_type(j,k,1) == POUT_) then

                   p0_g(slo(1):domlo(1)-1,j,k) = scale_pressure(bc_p_g(bcv))

               else if (bc_ilo_type(j,k,1) == MINF_) then

                   ! Note we index u_g differently to catch the inflow face
                   u_g(ulo(1):domlo(1)  ,j,k) = bc_u_g(bcv)
                   v_g(vlo(1):domlo(1)-1,j,k) = 0.0d0
                   w_g(wlo(1):domlo(1)-1,j,k) = 0.0d0

                   p0_g(slo(1):domlo(1)-1,j,k) = &
                       2.d0 * p0_g(domlo(1),j,k) - p0_g(domlo(1)+1,j,k)

               end if

            end do
         end do
      endif

      if (nrgt .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)

               bcv = bc_ihi_type(j,k,2)

               if (bc_ihi_type(j,k,1) == PINF_ .or. &
                   bc_ihi_type(j,k,1) == POUT_ .or. &
                   bc_ihi_type(j,k,1) == MINF_) then

                   if (is_undefined(ro_g0)) then
                      bc_ro_g = eosg(mw_avg,bc_p_g(bcv),bc_t_g(bcv))
                   else
                      bc_ro_g = ro_g0
                   endif

                   if (is_undefined(mu_g0)) then
                      bc_mu_g     = sutherland(bc_t_g(bcv))
                      bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                   else
                      bc_mu_g     = mu_g0
                      bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                   endif

                        ep_g(domhi(1)+1:shi(1),j,k) = bc_ep_g(bcv)
                        ro_g(domhi(1)+1:shi(1),j,k) = bc_ro_g
                       rop_g(domhi(1)+1:shi(1),j,k) = bc_ro_g*bc_ep_g(bcv)
                        mu_g(domhi(1)+1:shi(1),j,k) = bc_mu_g
                    lambda_g(domhi(1)+1:shi(1),j,k) = bc_lambda_g

               end if

               if (bc_ihi_type(j,k,1) == PINF_ .or. &
                   bc_ihi_type(j,k,1) == POUT_) then

                   p0_g(domhi(1)+1:shi(1),j,k) = scale_pressure(bc_p_g(bcv))

               else if (bc_ihi_type(j,k,1) == MINF_) then

                   ! Note we index the same on the high side
                   u_g(domhi(1)+1:uhi(1),j,k) = bc_u_g(bcv)
                   v_g(domhi(1)+1:vhi(1),j,k) = 0.0d0
                   w_g(domhi(1)+1:whi(1),j,k) = 0.0d0

                   p0_g(domhi(1)+1:shi(1),j,k) = &
                       2.d0 * p0_g(domhi(1),j,k) - p0_g(domhi(1)-1,j,k)
               end if

            end do
         end do
      endif

      if (nbot .gt. 0) then
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)

               bcv = bc_jlo_type(i,k,2)

               if (bc_jlo_type(i,k,1) == PINF_ .or. &
                   bc_jlo_type(i,k,1) == POUT_ .or. &
                   bc_jlo_type(i,k,1) == MINF_) then

                   if (is_undefined(ro_g0)) then
                      bc_ro_g = eosg(mw_avg,bc_p_g(bcv),bc_t_g(bcv))
                   else
                      bc_ro_g = ro_g0
                   endif

                   if (is_undefined(mu_g0)) then
                      bc_mu_g     = sutherland(bc_t_g(bcv))
                      bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                   else
                      bc_mu_g     = mu_g0
                      bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                   endif

                      ep_g(i,slo(2):domlo(2)-1,k) = bc_ep_g(bcv)
                      ro_g(i,slo(2):domlo(2)-1,k) = bc_ro_g
                     rop_g(i,slo(2):domlo(2)-1,k) = bc_ro_g*bc_ep_g(bcv)
                      mu_g(i,slo(2):domlo(2)-1,k) = bc_mu_g
                  lambda_g(i,slo(2):domlo(2)-1,k) = bc_lambda_g

               end if

               if (bc_jlo_type(i,k,1) == PINF_ .or. &
                   bc_jlo_type(i,k,1) == POUT_) then

                   p0_g(i,slo(2):domlo(2)-1,k) = scale_pressure(bc_p_g(bcv))

               else if (bc_jlo_type(i,k,1) == MINF_) then

                   ! Note we index v_g differently to catch the inflow face
                   u_g(i,ulo(2):domlo(2)-1,k) = 0.0d0
                   v_g(i,vlo(2):domlo(2)  ,k) = bc_v_g(bcv)
                   w_g(i,wlo(2):domlo(2)-1,k) = 0.0d0

                   p0_g(i,slo(2):domlo(2)-1,k) = &
                       2.d0 * p0_g(i,domlo(2),k) - p0_g(i,domlo(2)+1,k)

               end if

            end do
         end do
      endif

      if (ntop .gt. 0) then
         do k = slo(3),shi(3)
            do i = slo(1),shi(1)

               bcv = bc_jhi_type(i,k,2)

               if (bc_jhi_type(i,k,1) == PINF_ .or. &
                   bc_jhi_type(i,k,1) == POUT_ .or. &
                   bc_jhi_type(i,k,1) == MINF_) then

                   if (is_undefined(ro_g0)) then
                      bc_ro_g = eosg(mw_avg,bc_p_g(bcv),bc_t_g(bcv))
                   else
                      bc_ro_g = ro_g0
                   endif

                   if (is_undefined(mu_g0)) then
                      bc_mu_g     = sutherland(bc_t_g(bcv))
                      bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                   else
                      bc_mu_g     = mu_g0
                      bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                   endif

                      ep_g(i,domhi(2)+1:shi(2),k) = bc_ep_g(bcv)
                      ro_g(i,domhi(2)+1:shi(2),k) = bc_ro_g
                     rop_g(i,domhi(2)+1:shi(2),k) = bc_ro_g*bc_ep_g(bcv)
                      mu_g(i,domhi(2)+1:shi(2),k) = bc_mu_g
                  lambda_g(i,domhi(2)+1:shi(2),k) = bc_lambda_g

               end if

               if (bc_jhi_type(i,k,1) == PINF_ .or. &
                   bc_jhi_type(i,k,1) == POUT_) then

                   p0_g(i,domhi(2)+1:shi(2),k) = scale_pressure(bc_p_g(bcv))

               else if (bc_jhi_type(i,k,1) == MINF_) then

                   ! Note we index the same on the high side
                   u_g(i,domhi(2)+1:uhi(2),k) = 0.0d0
                   v_g(i,domhi(2)+1:vhi(2),k) = bc_v_g(bcv)
                   w_g(i,domhi(2)+1:whi(2),k) = 0.0d0

                   p0_g(i,domhi(2)+1:shi(2),k) = &
                       2.d0 * p0_g(i,domhi(2),k) - p0_g(i,domhi(2)-1,k)

               end if

            end do
         end do
      endif

      if (ndwn .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)

               bcv = bc_klo_type(i,j,2)

               if (bc_klo_type(i,j,1) == PINF_ .or. &
                   bc_klo_type(i,j,1) == POUT_ .or. &
                   bc_klo_type(i,j,1) == MINF_) then

                   if (is_undefined(ro_g0)) then
                      bc_ro_g = eosg(mw_avg,bc_p_g(bcv),bc_t_g(bcv))
                   else
                      bc_ro_g = ro_g0
                   endif

                   if (is_undefined(mu_g0)) then
                      bc_mu_g     = sutherland(bc_t_g(bcv))
                      bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                   else
                      bc_mu_g     = mu_g0
                      bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                   endif

                       ep_g(i,j,slo(3):domlo(3)-1) = bc_ep_g(bcv)
                       ro_g(i,j,slo(3):domlo(3)-1) = bc_ro_g
                      rop_g(i,j,slo(3):domlo(3)-1) = bc_ro_g*bc_ep_g(bcv)
                       mu_g(i,j,slo(3):domlo(3)-1) = bc_mu_g
                   lambda_g(i,j,slo(3):domlo(3)-1) = bc_lambda_g

               end if

               if (bc_klo_type(i,j,1) == PINF_ .or. &
                   bc_klo_type(i,j,1) == POUT_) then

                   p0_g(i,j,slo(3):domlo(3)-1) = scale_pressure(bc_p_g(bcv))

               else if (bc_klo_type(i,j,1) == MINF_) then

                   ! Note we index w_g differently to catch the inflow face
                   u_g(i,j,ulo(3):domlo(3)-1) = 0.0d0
                   v_g(i,j,vlo(3):domlo(3)-1) = 0.0d0
                   w_g(i,j,wlo(3):domlo(3)  ) = bc_w_g(bcv)

                   p0_g(i,j,slo(3):domlo(3)-1) = &
                       2.d0 * p0_g(i,j,domlo(3)) - p0_g(i,j,domlo(3)+1)

               end if

            end do
         end do
      endif

      if (nup .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)

               bcv = bc_khi_type(i,j,2)

               if (bc_khi_type(i,j,1) == PINF_ .or. &
                   bc_khi_type(i,j,1) == POUT_ .or. &
                   bc_khi_type(i,j,1) == MINF_) then

                   if (is_undefined(ro_g0)) then
                      bc_ro_g = eosg(mw_avg,bc_p_g(bcv),bc_t_g(bcv))
                   else
                      bc_ro_g = ro_g0
                   endif

                   if (is_undefined(mu_g0)) then
                      bc_mu_g     = sutherland(bc_t_g(bcv))
                      bc_lambda_g = -(2.0d0/3.0d0) * bc_mu_g
                   else
                      bc_mu_g     = mu_g0
                      bc_lambda_g = -(2.0d0/3.0d0) * mu_g0
                   endif

                       ep_g(i,j,domhi(3)+1:shi(3)) = bc_ep_g(bcv)
                       ro_g(i,j,domhi(3)+1:shi(3)) = bc_ro_g
                      rop_g(i,j,domhi(3)+1:shi(3)) = bc_ro_g*bc_ep_g(bcv)
                       mu_g(i,j,domhi(3)+1:shi(3)) = bc_mu_g
                   lambda_g(i,j,domhi(3)+1:shi(3)) = bc_lambda_g

               end if

               if (bc_khi_type(i,j,1) == PINF_ .or. &
                   bc_khi_type(i,j,1) == POUT_) then

                   p0_g(i,j,domhi(3)+1:shi(3)) = scale_pressure(bc_p_g(bcv))

               else if (bc_khi_type(i,j,1) == MINF_) then

                   ! Note we index the same on the high side
                   u_g(i,j,domhi(3)+1:uhi(3)) = 0.0d0
                   v_g(i,j,domhi(3)+1:vhi(3)) = 0.0d0
                   w_g(i,j,domhi(3)+1:whi(3)) = bc_w_g(bcv)

                   p0_g(i,j,domhi(3)+1:shi(3)) = &
                       2.d0 * p0_g(i,j,domhi(3)) - p0_g(i,j,domhi(3)-1)

               end if

            end do
         end do
      endif

   end subroutine set_bc0
