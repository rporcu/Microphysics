module source_pp_module

  use bl_fort_module, only : c_real
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
   subroutine source_pp_g(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
      A_m, b_m, b_mmax, dt, &
      u_g, v_g, w_g, p_g, ep_g, rop_g, rop_go, ro_g, d_e, d_n, d_t,&
      dx, dy, dz)

      use bc, only: small_number, one, zero, ijk_p_g
      use eos, only: droodp_g
      use fld_const, only: ro_g0
      use matrix, only: e, w, n, s, t, b
      use param1, only: is_defined, is_undefined
      use run, only: undefined_i
      use ur_facs, only: ur_fac
      use write_error_module, only: write_error

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer       , intent(in   ) ::  lo(3), hi(3)

      real(c_real), intent(in   ) :: dx, dy, dz, dt

      real(c_real), intent(inout) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: b_mmax&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: d_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local Variables
!-----------------------------------------------
! Indices
      integer :: i,j,k
! under relaxation factor for pressure
      real(c_real) fac
! terms of bm expression
      real(c_real) bma, bme, bmw, bmn, bms, bmt, bmb
! error message
      real(c_real) :: oDT, vol
!-----------------------------------------------

      odt = 1.0d0/dt
      vol = dx*dy*dz

! Calculate convection-diffusion fluxes through each of the faces

        do k = lo(3),hi(3)
           do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                bma = (rop_g(i,j,k)-rop_go(i,j,k))*vol*odt
                bme = A_m(i,j,k,e)*u_g(i,j,k)
                bmw = A_m(i,j,k,w)*u_g(i-1,j,k)
                bmn = A_m(i,j,k,n)*v_g(i,j,k)
                bms = A_m(i,j,k,s)*v_g(i,j-1,k)
                bmt = A_m(i,j,k,t)*w_g(i,j,k)
                bmb = A_m(i,j,k,b)*w_g(i,j,k-1)
                b_m(i,j,k) = -((-(bma + bme - bmw + bmn - bms + bmt - bmb )) )
                b_mmax(i,j,k) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), &
                   abs(bms), abs(bmt), abs(bmb))

                A_m(i,j,k,e) = A_m(i,j,k,e)*d_e(i,j,k)
                A_m(i,j,k,w) = A_m(i,j,k,w)*d_e(i-1,j,k)
                A_m(i,j,k,n) = A_m(i,j,k,n)*d_n(i,j,k)
                A_m(i,j,k,s) = A_m(i,j,k,s)*d_n(i,j-1,k)
                A_m(i,j,k,t) = A_m(i,j,k,t)*d_t(i,j,k)
                A_m(i,j,k,b) = A_m(i,j,k,b)*d_t(i,j,k-1)

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

! Specify P' to zero for incompressible flows. Check set_bc0
! for details on selection of IJK_P_g.
      if (ijk_p_g(1) /= undefined_i) then
         b_m(ijk_p_g(1),ijk_p_g(2),ijk_p_g(3)) = zero
         A_m(ijk_p_g(1),ijk_p_g(2),ijk_p_g(3),:) = zero
         A_m(ijk_p_g(1),ijk_p_g(2),ijk_p_g(3),0) = -one
      endif

      return

   end subroutine source_pp_g


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOURCE_Pp_g_bc                                          !
!  Purpose: Determine source terms for Pressure correction equation.   !
!                                                                      !
!  Notes: The off-diagonal coefficients are positive. The center       !
!         coefficient and the source vector are negative. See          !
!         conv_Pp_g                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine source_pp_g_bc(slo, shi, A_m, bc_ilo_type, bc_ihi_type, &
      bc_jlo_type, bc_jhi_type, bc_klo_type, bc_khi_type)

      use ic, only: PINF_, POUT_
      use geometry, only: domlo, domhi
      use matrix, only: e, n, t, w, s, b
      use param1, only: zero

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      real(c_real), intent(inout) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (slo(1):shi(1),slo(2):shi(2),2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (slo(1):shi(1),slo(2):shi(2),2)

! Local Variables
!-----------------------------------------------
      integer :: i,j,k
      integer :: nlft, nrgt, nbot, ntop, nup, ndwn
!-----------------------------------------------

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))


! --- EAST FLUID ---------------------------------------------------------->

      if (nlft .gt. 0) then
         i = domlo(1)
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               if(bc_ilo_type(j,k,1) == PINF_ .or. &
                  bc_ilo_type(j,k,1) == POUT_) then
                  A_m(i,j,k,w) =  zero
               endif
            end do
         end do
      endif

! --- WEST FLUID ---------------------------------------------------------->

      if (nrgt .gt. 0) then
         i = domhi(1)
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               if(bc_ihi_type(j,k,1) == PINF_ .or. &
                  bc_ihi_type(j,k,1) == POUT_) then
                  A_m(i,j,k,e) = zero
               endif
            end do
         end do
      endif

! --- NORTH FLUID --------------------------------------------------------->

      if (nbot .gt. 0) then
         j = domlo(2)
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               if(bc_jlo_type(i,k,1) == PINF_ .or. &
                  bc_jlo_type(i,k,1) == POUT_) then
                  A_m(i,j,k,s) = zero
               endif
            end do
         end do
      endif


! --- SOUTH FLUID --------------------------------------------------------->

      if (ntop .gt. 0) then
         j = domhi(2)
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               if(bc_jhi_type(i,k,1) == PINF_ .or. &
                  bc_jhi_type(i,k,1) == POUT_) then
                  A_m(i,j,k,n) = zero
               endif
            end do
         end do
      endif

! --- TOP FLUID ----------------------------------------------------------->

      if (ndwn .gt. 0) then
         k = domlo(3)
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               if(bc_klo_type(i,j,1) == PINF_ .or. &
                  bc_klo_type(i,j,1) == POUT_) then
                  A_m(i,j,k,b) = zero
               endif
            end do
         end do
      endif

! --- BOTTOM FLUID -------------------------------------------------------->

      if (nup .gt. 0) then
         k = domhi(3)
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               if(bc_khi_type(i,j,1) == PINF_ .or. &
                  bc_khi_type(i,j,1) == POUT_) then
                  A_m(i,j,k,t) = zero
               endif
            end do
         end do
      endif

      return

   end subroutine source_pp_g_bc

end module source_pp_module
