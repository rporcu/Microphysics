module source_pp_module

  use amrex_fort_module, only : c_real => amrex_real
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
      use param1, only: is_defined, is_undefined
      use ur_facs, only: ur_fac

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      integer       , intent(in   ) ::  lo(3), hi(3)

      real(c_real), intent(in   ) :: dx, dy, dz, dt

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(inout) :: b_mmax&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

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
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

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

        do k = alo(3),ahi(3)
           do j = alo(2),ahi(2)
             do i = alo(1),ahi(1)

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
!  Subroutine: SOURCE_Pp_g_bc                                          !
!  Purpose: Determine source terms for Pressure correction equation.   !
!                                                                      !
!  Notes: The off-diagonal coefficients are positive. The center       !
!         coefficient and the source vector are negative. See          !
!         conv_Pp_g                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine source_pp_g_bc(slo, shi, alo, ahi, domlo, domhi, A_m)

      use matrix, only: e, n, t, w, s, b
      use param1, only: zero

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),alo(3),ahi(3),domlo(3),domhi(3)

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

! --- EAST FLUID ---------------------------------------------------------->

      if (slo(1).lt.domlo(1)) &
         A_m(domlo(1),alo(2):ahi(2),alo(3):ahi(3),w) =  zero

! --- WEST FLUID ---------------------------------------------------------->

      if (shi(1).gt.domhi(1)) &
         A_m(domhi(1),alo(2):ahi(2),alo(3):ahi(3),e) =  zero

! --- NORTH FLUID --------------------------------------------------------->

      if (slo(2).lt.domlo(2)) &
         A_m(alo(1):ahi(1),domlo(2),alo(3):ahi(3),s) =  zero

! --- SOUTH FLUID --------------------------------------------------------->

      if (shi(2).gt.domhi(2)) &
         A_m(alo(1):ahi(1),domhi(2),alo(3):ahi(3),n) =  zero

! --- TOP FLUID ----------------------------------------------------------->

      if (slo(3).lt.domlo(3)) &
         A_m(alo(1):ahi(1),alo(2):ahi(2),domlo(3),b) =  zero

! --- BOTTOM FLUID -------------------------------------------------------->

      if (shi(3).gt.domhi(3)) &
         A_m(alo(1):ahi(1),alo(2):ahi(2),domhi(3),t) =  zero

   end subroutine source_pp_g_bc

end module source_pp_module
