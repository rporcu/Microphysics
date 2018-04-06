module source_w_g_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int
   use param        , only: zero, half, one, undefined, is_undefined

   implicit none

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOURCE_W_g                                              !
!  Purpose: Determine source terms for W_g momentum eq. The terms      !
!     appear in the center coefficient and RHS vector. The center      !
!     coefficient and source vector are negative.  The off-diagonal    !
!     coefficients are positive.                                       !
!     The drag terms are excluded from the source at this stage.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine source_w_g(lo, hi, slo, shi, wlo, whi, alo, ahi, dlo, dhi, &
       A_m, b_m, p_g, p0_g, ep_g, ro_g, rop_go, w_go, tau_w_g, f_gds_w, drag_w, &
       dt, dx, dy, dz, domlo, domhi)

! Modules
!---------------------------------------------------------------------//
      use constant, only: gravity

      use matrix, only: e, w, s, n, t, b
      use scales, only: p_scale

      integer     , intent(in   ) :: lo(3),    hi(3)
      integer     , intent(in   ) :: slo(3),   shi(3)
      integer     , intent(in   ) :: wlo(3),   whi(3)
      integer     , intent(in   ) :: alo(3),   ahi(3)
      integer     , intent(in   ) :: dlo(3),   dhi(3)
      integer     , intent(in   ) :: domlo(3), domhi(3)

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: &
           A_m(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3), &
           b_m(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      real(c_real), intent(in   ) :: &
           p_g    (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           p0_g   (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           ep_g   (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           ro_g   (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           rop_go (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &

           w_go   (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           tau_w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &

           f_gds_w(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)), &
           drag_w (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      integer :: i,j,k
! Pressure at top cell
      real(c_real) :: PgB
! Average volume fraction
      real(c_real) :: epga
! Source terms (Surface)
      real(c_real) Sdp
! Source terms (Volumetric)
      real(c_real) V0, Vbf
! jackson terms: local stress tensor quantity
      real(c_real) :: odt
      real(c_real) :: axy, vol
!---------------------------------------------------------------------//

      odt = 1.0d0/dt
      axy = dx*dy
      vol = dx*dy*dz

      do k =  lo(3), hi(3)
         do j =  lo(2), hi(2)
            do i =  lo(1), hi(1)

               epga = half*(ep_g(i,j,k-1) + ep_g(i,j,k))

               ! Pressure gradient term
               sdp = -p_scale*epga*( (p_g(i,j,k  )+p0_g(i,j,k  ))- &
                                     (p_g(i,j,k-1)+p0_g(i,j,k-1)) )*axy

               ! Previous time step
               v0 = half*(rop_go(i,j,k-1) + rop_go(i,j,k))*odt

               ! Body force
               vbf = half*(ep_g(i,j,k-1) * ro_g(i,j,k-1) + &
                    &      ep_g(i,j,k)   * ro_g(i,j,k)  ) *gravity(3)

               ! Collect the terms
               A_m(i,j,k,0) = -(A_m(i,j,k,e) + A_m(i,j,k,w) + &
                                A_m(i,j,k,n) + A_m(i,j,k,s) + &
                                A_m(i,j,k,t) + A_m(i,j,k,b) + &
                                v0*vol + f_gds_w(i,j,k))

               b_m(i,j,k) = -(sdp + tau_w_g(i,j,k) + drag_w(i,j,k) + &
                    ( (v0)*w_go(i,j,k) + vbf)*vol)

            enddo
         enddo
      enddo

      end subroutine source_w_g


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOURCE_W_g_BC                                           !
!  Purpose: Determine source terms for W_g momentum eq. The terms      !
!     appear in the center coefficient and RHS vector. The center      !
!     coefficient and source vector are negative.  The off-diagonal    !
!     coefficients are positive.                                       !
!     The drag terms are excluded from the source at this stage        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine source_w_g_bc(lo,hi,slo,shi,alo,ahi,A_m,b_m, &
         bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
         bc_klo_type, bc_khi_type, domlo, domhi, ng, dx, dy)


      use bc, only: nsw_, fsw_, psw_
      use bc, only: pinf_, pout_
      use bc, only: minf_

      use bc, only: bc_hw_g, bc_ww_g, bc_w_g

      use matrix, only: e, w, s, n, t, b
      use param, only: is_defined

      integer     , intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: slo(3),shi(3),alo(3),ahi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3), ng
      real(c_real), intent(in   ) :: dx, dy

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      integer(c_int), intent(in   ) :: &
           bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

! Local variables
!-----------------------------------------------
      real(c_real) :: odx,ody
      integer      :: bcv,kbc,i,j,k
!---------------------------------------------------------------------//

      odx = 1.0d0/dx
      ody = 1.0d0/dy

! ------------------------------------------------------------->

      ! At left boundary
      if (slo(1) .lt. domlo(1) .and. lo(1).eq.alo(1)) then

         i = alo(1)

         do k = lo(3), hi(3)

           ! bc's on i-faces only defined within (domlo(2):domhi(3),domlo(2):domhi(3))
            kbc = max(min(k,domhi(3)),domlo(3))

            do j = lo(2), hi(2)

               bcv = bc_ilo_type(j,kbc,2)

               if(bc_ilo_type(j,kbc,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,w)
                  A_m(i,j,k,w) = zero

               else if(bc_ilo_type(j,kbc,1) == FSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,w)
                  A_m(i,j,k,w) = zero

               else if(bc_ilo_type(j,kbc,1) == PSW_) then
                  if (is_undefined(bc_hw_g(bcv))) then
                     A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,w)
                     b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,w)*bc_ww_g(bcv)
                  else
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,w)*&
                        (half*bc_hw_g(bcv)-odx)/(half*bc_hw_g(bcv)+odx)
                     b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,w)*&
                        bc_hw_g(bcv)*bc_ww_g(bcv)/(half*bc_hw_g(bcv)+odx)
                  endif
                  A_m(i,j,k,w) = zero

               else if(bc_ilo_type(j,kbc,1) == PINF_ .or. &
                       bc_ilo_type(j,kbc,1) == POUT_ .or. &
                       bc_ilo_type(j,kbc,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif
            end do
         end do
      endif

! ------------------------------------------------------------->

      ! At right boundary
      if (shi(1) .gt. domhi(1) .and. hi(1).eq.ahi(1)) then

         i = ahi(1)
         do k = lo(3), hi(3)

            ! bc's on i-faces only defined within (domlo(2):domhi(3),domlo(2):domhi(3))
            kbc = max(min(k,domhi(3)),domlo(3))

            do j= lo(2), hi(2)

               bcv = bc_ihi_type(j,kbc,2)

               if(bc_ihi_type(j,kbc,1) == NSW_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,e)
                  A_m(i,j,k,e) = zero

               else if(bc_ihi_type(j,kbc,1) == FSW_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,e)
                  A_m(i,j,k,e) = zero

               else if(bc_ihi_type(j,kbc,1) == PSW_) then

                  if (is_undefined(bc_hw_g(bcv))) then
                     A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,e)
                     b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,e)*bc_ww_g(bcv)
                  else
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,e)*&
                        (half*bc_hw_g(bcv)-odx)/(half*bc_hw_g(bcv)+odx)
                     b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,e)*&
                        bc_hw_g(bcv)*bc_ww_g(bcv)/(half*bc_hw_g(bcv)+odx)
                  endif
                  A_m(i,j,k,e) = zero

               else if(bc_ihi_type(j,kbc,1) == PINF_ .or. &
                       bc_ihi_type(j,kbc,1) == POUT_ .or. &
                       bc_ihi_type(j,kbc,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

            end do
         end do
      endif

! ------------------------------------------------------------->

      ! At bottom boundary
      if (slo(2) .lt. domlo(2) .and. lo(2).eq.alo(2)) then

         j = alo(2)

         do k = lo(3), hi(3)

            ! bc's on j-faces only defined within (domlo(1):domhi(1),domlo(3):domhi(3))
            kbc = max(min(k,domhi(3)),domlo(3))

            do i = lo(1), hi(1)

               bcv = bc_jlo_type(i,kbc,2)

               if(bc_jlo_type(i,kbc,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,s)
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(i,kbc,1) == FSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,s)
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(i,kbc,1) == PSW_) then
                  if (is_undefined(bc_hw_g(bcv))) then
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,s)
                     b_m(i,j,k) = b_m(i,j,k) - 2.0*A_m(i,j,k,s)*bc_ww_g(bcv)
                  else
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,s)*&
                        (half*bc_hw_g(bcv)-ody)/(half*bc_hw_g(bcv)+ody)
                     b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,s)*&
                        bc_hw_g(bcv)*bc_ww_g(bcv)/(half*bc_hw_g(bcv)+ody)
                  endif
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(i,kbc,1) == PINF_ .or. &
                       bc_jlo_type(i,kbc,1) == POUT_ .or. &
                       bc_jlo_type(i,kbc,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

            end do
         end do
      endif

! ------------------------------------------------------------->

      ! At top boundary
      if (shi(2) .gt. domhi(2) .and. hi(2).eq.ahi(2)) then

         j = ahi(2)

         do k = lo(3), hi(3)

            ! bc's on j-faces only defined within (domlo(1):domhi(1),domlo(3):domhi(3))
            kbc = max(min(k,domhi(3)),domlo(3))

            do i = lo(1), hi(1)

               bcv = bc_jhi_type(i,kbc,2)

               if(bc_jhi_type(i,kbc,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,n)
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(i,kbc,1) == FSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,n)
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(i,kbc,1) == PSW_) then
                  if (is_undefined(bc_hw_g(bcv))) then
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,n)
                     b_m(i,j,k) = b_m(i,j,k) - 2.0*A_m(i,j,k,n)*bc_ww_g(bcv)
                  else
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,n)*&
                        (half*bc_hw_g(bcv)-ody)/(half*bc_hw_g(bcv)+ody)
                     b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,n)*&
                        bc_hw_g(bcv)*bc_ww_g(bcv)/(half*bc_hw_g(bcv)+ody)
                  endif
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(i,kbc,1) == PINF_ .or. &
                       bc_jhi_type(i,kbc,1) == POUT_ .or. &
                       bc_jhi_type(i,kbc,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

            end do
         end do
      endif

! ------------------------------------------------------------->

      ! At down boundary
      if (slo(3) .lt. domlo(3) .and. lo(3).eq.alo(3)) then

         k = alo(3)

         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               bcv = bc_klo_type(i,j,2)

               if(bc_klo_type(i,j,1) == PINF_ .or. &
                  bc_klo_type(i,j,1) == POUT_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,b)
                  A_m(i,j,k,b) = zero

               else if (bc_klo_type(i,j,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k  ) = -bc_w_g(bcv)

               else if(bc_klo_type(i,j,1) == NSW_ .or. &
                       bc_klo_type(i,j,1) == FSW_ .or. &
                       bc_klo_type(i,j,1) == PSW_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif
            end do
         end do
      endif

! ------------------------------------------------------------->

      ! At up boundary
      if (shi(3) .gt. domhi(3) .and. hi(3).eq.ahi(3)) then

         k = ahi(3)

         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               bcv = bc_khi_type(i,j,2)

               if(bc_khi_type(i,j,1) == PINF_ .or. &
                  bc_khi_type(i,j,1) == POUT_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,t)
                  A_m(i,j,k,t) = zero

               else if(bc_khi_type(i,j,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k  ) = -bc_w_g(bcv)

               else if(bc_khi_type(i,j,1) == NSW_ .or. &
                       bc_khi_type(i,j,1) == FSW_ .or. &
                       bc_khi_type(i,j,1) == PSW_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

            end do
         end do
      endif

      end subroutine source_w_g_bc


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_W_G                                        C
!  Purpose: Adds point sources to the gas phase W-Momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine point_source_w_g(lo, hi, alo,ahi,b_m,vol)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      ! use param  , only: small_number, zero
      ! use ps, only: dim_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      ! use ps, only: ps_w_g

      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      real(c_real)  , intent(inout) :: b_m(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real)  , intent(in   ) :: vol

      ! integer(c_int) :: I, J, K
      ! integer(c_int) :: psv
      ! integer(c_int) :: lKT, lKB
      ! real(c_real)  :: pSource

      ! ! Calculate the mass going into each (i,j,k) cell. This is done for each
      ! ! call in case the point source is time dependent.
      ! ps_lp: do psv = 1, dim_ps
      !    if (.not.ps_defined(psv)) cycle ps_lp
      !    if (abs(PS_W_g(psv)) < small_number) cycle ps_lp

      !    if(PS_W_g(psv) < ZERO) then
      !       lKB = PS_K_B(psv)-1
      !       lKT = PS_K_T(psv)-1
      !    else
      !       lKB = PS_K_B(psv)
      !       lKT = PS_K_T(psv)
      !    endif

      !    do k = lKB, lKT
      !    do j = PS_J_S(psv), PS_J_N(psv)
      !    do i = PS_I_W(psv), PS_I_E(psv)

      !       pSource =  PS_MASSFLOW_G(psv) * (VOL/PS_VOLUME(psv))

      !       b_m(I,J,K) = b_m(I,J,K) - pSource * &
      !          PS_W_g(psv) * PS_VEL_MAG_g(psv)

      !    enddo
      !    enddo
      !    enddo

      ! enddo ps_lp

      end subroutine point_source_w_g
end module source_w_g_module
