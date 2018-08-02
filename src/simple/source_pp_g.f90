module source_pp_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOURCE_Pp_g                                             !
!  Purpose: Determine source terms for Pressure correction equation.   !
!                                                                      !
!  Notes: The off-diagonal coefficients are positive. The center       !
!         coefficient and the source vector are negative. See          !
!         conv_Pp_g                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine source_pp_g(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, lo, hi, &
      A_m, b_m, b_mmax, dt, &
      u_g, v_g, w_g, p_g, ep_g, rop_g, rop_go, ro_g, d_e, d_n, d_t,&
      dx, dy, dz)

      use eos, only: droodp_g
      use fld_const, only: ro_g0
      use matrix, only: e, w, n, s, t, b
      use param, only: is_defined, is_undefined
      use ur_facs, only: ur_fac

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      integer       , intent(in   ) ::  lo(3), hi(3)

      real(rt), intent(in   ) :: dx, dy, dz, dt

      real(rt), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(rt), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(rt), intent(inout) :: b_mmax&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      real(rt), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(rt), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(rt), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(rt), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in   ) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(rt), intent(in   ) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(rt), intent(in   ) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

!-----------------------------------------------
      integer :: i,j,k

! under relaxation factor for pressure
      real(rt) fac
! terms of bm expression
      real(rt) bma, bme, bmw, bmn, bms, bmt, bmb
! error message
      real(rt) :: odt, vol
!-----------------------------------------------

      odt = 1.0d0/dt
      vol = dx*dy*dz

! Calculate convection-diffusion fluxes through each of the faces

        do k = lo(3),hi(3)
           do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                bma = (rop_g(i,j,k)-rop_go(i,j,k))*vol*odt
                bme = A_m(i,j,k,e)*u_g(i+1,j,k)
                bmw = A_m(i,j,k,w)*u_g(i  ,j,k)
                bmn = A_m(i,j,k,n)*v_g(i,j+1,k)
                bms = A_m(i,j,k,s)*v_g(i,j  ,k)
                bmt = A_m(i,j,k,t)*w_g(i,j,k+1)
                bmb = A_m(i,j,k,b)*w_g(i,j,k  )

                b_m(i,j,k) = -((-(bma + bme - bmw + bmn - bms + bmt - bmb )) )

                b_mmax(i,j,k) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), &
                   abs(bms), abs(bmt), abs(bmb))

                A_m(i,j,k,e) = A_m(i,j,k,e)*d_e(i+1,j,k)
                A_m(i,j,k,w) = A_m(i,j,k,w)*d_e(i  ,j,k)
                A_m(i,j,k,n) = A_m(i,j,k,n)*d_n(i,j+1,k)
                A_m(i,j,k,s) = A_m(i,j,k,s)*d_n(i,j  ,k)
                A_m(i,j,k,t) = A_m(i,j,k,t)*d_t(i,j,k+1)
                A_m(i,j,k,b) = A_m(i,j,k,b)*d_t(i,j,k  )

                A_m(i,j,k,0) = -(A_m(i,j,k,e) + A_m(i,j,k,w) + &
                                 A_m(i,j,k,n) + A_m(i,j,k,s) + &
                                 A_m(i,j,k,t) + A_m(i,j,k,b))

             enddo
          enddo
      enddo

! make correction for compressible flows
      if (is_undefined(ro_g0)) then
         fac = ur_fac(1)  !since p_g = p_g* + ur_fac * pp_g
         do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 A_m(i,j,k,0) = A_m(i,j,k,0) - &
                    fac*droodp_g(ro_g(i,j,k),p_g(i,j,k))*&
                    ep_g(i,j,k)*vol*odt
              end do
           end do
        end do

     endif   ! end if (ro_g0 == undefined); i.e., compressible flow

   end subroutine source_pp_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: source_pp_g_bc                                          !
!  Purpose: Determine source terms for Pressure correction equation.   !
!                                                                      !
!  Notes: The off-diagonal coefficients are positive. The center       !
!         coefficient and the source vector are negative. See          !
!         conv_Pp_g                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   subroutine source_pp_g_bc(alo, ahi, lo, hi, domlo, domhi, ng, A_m, &
                             bc_ilo_type, bc_ihi_type, &
                             bc_jlo_type, bc_jhi_type, &
                             bc_klo_type, bc_khi_type )

      use bc, only: POUT_
      use matrix, only: e, n, t, w, s, b
      use param, only: zero

      implicit none

      integer     , intent(in   ) :: alo(3),ahi(3),lo(3),hi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3), ng

      real(rt), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)


      integer(c_int), intent(in   ) :: &
           bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      integer(c_int) :: i,j,k

!-----------------------------------------------

      ! At west boundary
      if (lo(1) .eq. domlo(1) .and. lo(1).eq.alo(1)) then
         i = alo(1)
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            if (bc_ilo_type(j,k,1) == POUT_) A_m(i,j,k,w) = zero
         end do
         end do
      endif

      ! At east boundary
      if (hi(1) .eq. domhi(1) .and. hi(1).eq.ahi(1)) then
         i = ahi(1)
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            if (bc_ihi_type(j,k,1) == POUT_) A_m(i,j,k,e) = zero
         end do
         end do
      endif

      ! At south boundary
      if (lo(2) .eq. domlo(2) .and. lo(2).eq.alo(2)) then
         j = alo(2)
         do k = lo(3),hi(3)
         do i = lo(1),hi(1)
            if (bc_jlo_type(i,k,1) == POUT_) A_m(i,j,k,s) = zero
         end do
         end do
      endif

      ! At north boundary
      if (hi(2) .eq. domhi(2) .and. hi(2).eq.ahi(2)) then
         j = ahi(2)
         do k = lo(3),hi(3)
         do i = lo(1),hi(1)
            if (bc_jhi_type(i,k,1) == POUT_) A_m(i,j,k,n) = zero
         end do
         end do
      endif

      ! At bottom boundary
      if (lo(3) .eq. domlo(3) .and. lo(3).eq.alo(3)) then
         k = alo(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (bc_klo_type(i,j,1) == POUT_) A_m(i,j,k,b) = zero
         end do
         end do
      endif

      ! At top boundary
      if (hi(3) .eq. domhi(3) .and. hi(3).eq.ahi(3)) then
         k = ahi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (bc_khi_type(i,j,1) == POUT_) A_m(i,j,k,t) = zero
         end do
         end do
      endif

   end subroutine source_pp_g_bc

end module source_pp_module
