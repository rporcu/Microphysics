!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: set_bc0                                                 !
!  Purpose: This subroutine does the initial setting of all boundary   !
!  conditions. The user specifications of the boundary conditions are  !
!  checked for veracity in various check_data/ routines:               !
!  (e.g., check_boundary_conditions).                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine set_bc0(slo, shi, &
                      ep_g, ro_g, rop_g, mu_g, lambda_g, &
                      bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                      bc_klo_type, bc_khi_type, domlo, domhi, ng) &
      bind(C, name="set_bc0")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use bc, only: bc_ep_g, bc_t_g
      use bc, only: pinf_, pout_, minf_

      use eos      , only: eosg, sutherland
      use fld_const, only: ro_g0, mu_g0

      use scales, only: scale_pressure
      use param , only: is_undefined

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)
      integer(c_int), intent(in   ) :: ng

      real(rt), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: ro_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: rop_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: mu_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: lambda_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

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

      integer    nlft, nrgt, nbot, ntop, nup, ndwn

      real(rt) :: bc_ro_g, bc_mu_g, bc_lambda_g
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

                  bc_ro_g = ro_g0

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

                   bc_ro_g = ro_g0

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

                   bc_ro_g = ro_g0

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

                   bc_ro_g = ro_g0

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

                   bc_ro_g = ro_g0

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

                   bc_ro_g = ro_g0

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

            end do
         end do
      endif

   end subroutine set_bc0
