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
                         u_g, v_g, w_g, p_g, ep_g, ro_g, rop_g, mu_g, lambda_g, &
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
      integer    ilo, ihi, jlo, jhi, klo, khi

      real(c_real) :: bc_ro_g, bc_mu_g, bc_lambda_g
!--------------------------------------------------------------------//

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))


      if (nlft .gt. 0) then
         ilo = domlo(1)
         do i = 1, nlft
            do k=slo(3),shi(3)
               do j=slo(2),shi(2)
                  bcv = bc_ilo_type(j,k,2)
                  if(bc_ilo_type(j,k,1) == PINF_ .or. &
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

                     p_g(ilo-i,j,k) = scale_pressure(bc_p_g(bcv))
                     ep_g(ilo-i,j,k) = bc_ep_g(bcv)

                     ro_g(ilo-i,j,k) = bc_ro_g
                     rop_g(ilo-i,j,k) = bc_ro_g*bc_ep_g(bcv)

                     mu_g(ilo-i,j,k) = bc_mu_g

                     if(bc_ilo_type(j,k,1) == PINF_ .or. &
                        bc_ilo_type(j,k,1) == MINF_) then

                        u_g(ilo-i,j,k) = bc_u_g(bcv)
                        v_g(ilo-i,j,k) = 0.0d0
                        w_g(ilo-i,j,k) = 0.0d0

                     end if
                  end if
               end do
            end do
         end do
      endif

      if (nrgt .gt. 0) then
         ihi = domhi(1)
         do i = 1, nrgt
            do k=slo(3),shi(3)
               do j=slo(2),shi(2)
                  bcv = bc_ihi_type(j,k,2)
                  if(bc_ihi_type(j,k,1) == PINF_ .or. &
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

                     p_g(ihi+i,j,k) = scale_pressure(bc_p_g(bcv))
                     ep_g(ihi+i,j,k) = bc_ep_g(bcv)

                     ro_g(ihi+i,j,k) = bc_ro_g
                     rop_g(ihi+i,j,k) = bc_ro_g*bc_ep_g(bcv)

                     mu_g(ihi+i,j,k) = bc_mu_g

                     if(bc_ihi_type(j,k,1) == PINF_ .or. &
                        bc_ihi_type(j,k,1) == MINF_) then

                        u_g(ihi+i-1,j,k) = bc_u_g(bcv)
                        v_g(ihi+i-1,j,k) = 0.0d0
                        w_g(ihi+i-1,j,k) = 0.0d0
                     end if
                  end if
               end do
            end do
         end do
      endif

      if (nbot .gt. 0) then
         jlo = domlo(2)
         do j = 1, nbot
            do k=slo(3),shi(3)
               do i=slo(1),shi(1)
                  bcv = bc_jlo_type(i,k,2)
                  if(bc_jlo_type(i,k,1) == PINF_ .or. &
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

                     p_g(i,jlo-j,k) = scale_pressure(bc_p_g(bcv))
                     ep_g(i,jlo-j,k) = bc_ep_g(bcv)

                     ro_g(i,jlo-j,k) = bc_ro_g
                     rop_g(i,jlo-j,k) = bc_ro_g*bc_ep_g(bcv)

                     mu_g(i,jlo-j,k) = bc_mu_g

                     if(bc_jlo_type(i,k,1) == PINF_ .or. &
                        bc_jlo_type(i,k,1) == MINF_) then

                        u_g(i,jlo-j,k) = 0.0d0
                        v_g(i,jlo-j,k) = bc_v_g(bcv)
                        w_g(i,jlo-j,k) = 0.0d0

                     end if
                  end if
               end do
            end do
         end do
      endif

      if (ntop .gt. 0) then
         jhi = domhi(2)
         do j = 1, ntop
            do k=slo(3),shi(3)
               do i=slo(1),shi(1)
                  bcv = bc_jhi_type(i,k,2)
                  if(bc_jhi_type(i,k,1) == PINF_ .or. &
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

                     p_g(i,jhi+j,k) = scale_pressure(bc_p_g(bcv))
                     ep_g(i,jhi+j,k) = bc_ep_g(bcv)

                     ro_g(i,jhi+j,k) = bc_ro_g
                     rop_g(i,jhi+j,k) = bc_ro_g*bc_ep_g(bcv)

                     mu_g(i,jhi+j,k) = bc_mu_g

                     if(bc_jhi_type(i,k,1) == PINF_ .or. &
                        bc_jhi_type(i,k,1) == MINF_) then

                        u_g(i,jhi+j-1,k) = 0.0d0
                        v_g(i,jhi+j-1,k) = bc_v_g(bcv)
                        w_g(i,jhi+j-1,k) = 0.0d0

                     end if
                  end if
               end do
            end do
         end do
      endif

      if (ndwn .gt. 0) then
         klo = domlo(3)
         do k = 1, ndwn
            do j=slo(2),shi(2)
               do i=slo(1),shi(1)
                  bcv = bc_klo_type(i,j,2)
                  if(bc_klo_type(i,j,1) == PINF_ .or. &
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

                     p_g(i,j,klo-k) = scale_pressure(bc_p_g(bcv))
                     ep_g(i,j,klo-k) = bc_ep_g(bcv)

                     ro_g(i,j,klo-k) = bc_ro_g
                     rop_g(i,j,klo-k) = bc_ro_g*bc_ep_g(bcv)

                     mu_g(i,j,klo-k) = bc_mu_g

                     if(bc_klo_type(i,j,1) == PINF_ .or. &
                          bc_klo_type(i,j,1) == MINF_) then

                        u_g(i,j,klo-k) = 0.0d0
                        v_g(i,j,klo-k) = 0.0d0
                        w_g(i,j,klo-k) = bc_w_g(bcv)
                     end if
                  end if
               end do
            end do
         end do
      endif

      if (nup .gt. 0) then
         khi = domhi(3)
         do k = 1, nup
            do j=slo(2),shi(2)
               do i=slo(1),shi(1)
                  bcv = bc_khi_type(i,j,2)
                  if(bc_khi_type(i,j,1) == PINF_ .or. &
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

                     p_g(i,j,khi+k) = scale_pressure(bc_p_g(bcv))
                     ep_g(i,j,khi+k) = bc_ep_g(bcv)

                     ro_g(i,j,khi+k) = bc_ro_g
                     rop_g(i,j,khi+k) = bc_ro_g*bc_ep_g(bcv)

                     mu_g(i,j,khi+k) = bc_mu_g

                     if(bc_khi_type(i,j,1) == PINF_ .or. &
                          bc_khi_type(i,j,1) == MINF_) then

                        u_g(i,j,khi+k-1) = 0.0d0
                        v_g(i,j,khi+k-1) = 0.0d0
                        w_g(i,j,khi+k-1) = bc_w_g(bcv)
                     end if
                  end if
               end do
            end do
         end do
      endif

   end subroutine set_bc0
