module source_u_g_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int
   use param        , only: zero, half, one, undefined, is_undefined, small_number

   implicit none

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOURCE_U_g                                              !
!  Purpose: Determine source terms for U_g momentum eq. The terms      !
!  appear in the center coefficient and RHS vector. The center         !
!  coefficient and source vector are negative.  The off-diagonal       !
!  coefficients are positive.                                          !
!  The drag terms are excluded from the source at this stage.          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine source_u_g(lo, hi, slo, shi, ulo, uhi, alo, ahi, A_m, b_m, &
      dt, p_g, ep_g, ro_g, rop_go, u_go, &
      tau_u_g, dx, dy, dz, domlo, domhi)


! Modules
!---------------------------------------------------------------------//
      use constant, only: gravity
      use bc      , only: delp_x

      use functions, only: avg

      use matrix, only: e, w, s, n, t, b
      use scales, only: p_scale

      integer     , intent(in   ) :: lo(3),hi(3)
      integer     , intent(in   ) :: slo(3),shi(3),ulo(3),uhi(3),alo(3),ahi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      ! Vector b_m
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_go&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: tau_u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

! Local Variables
!---------------------------------------------------------------------//
! Indices
      integer :: i,j,k
! Pressure at east cell
      real(c_real) :: PgW
! Average volume fraction
      real(c_real) :: epga
! Source terms (Surface)
      real(c_real) :: Sdp
! Source terms (Volumetric)
      real(c_real) :: V0, Vbf
! local stress tensor quantity
      real(c_real) :: odt
      real(c_real) :: ayz
      real(c_real) :: vol
!---------------------------------------------------------------------//

      odt = 1.0d0/dt
      ayz = dy*dz
      vol = dx*dy*dz

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               epga = half*(ep_g(i-1,j,k) + ep_g(i,j,k))

               ! Pressure term
               pgw = p_g(i-1,j,k)
               if(i == domlo(1) .or. i==domhi(1)+1) pgw = pgw + delp_x
               sdp = -p_scale*epga*(p_g(i,j,k) - pgw)*ayz

               ! Previous time step
               v0 = half * (rop_go(i-1,j,k) + rop_go(i,j,k))*odt

               ! Body force
               vbf = half*(ro_g(i-1,j,k) + ro_g(i,j,k))*gravity(1)

               ! Collect the terms
               A_m(i,j,k,0) = -(A_m(i,j,k,e) + A_m(i,j,k,w) + &
                                A_m(i,j,k,n) + A_m(i,j,k,s) + &
                                A_m(i,j,k,t) + A_m(i,j,k,b) + &
                                v0*vol)

               b_m(i,j,k) = b_m(i,j,k) -(sdp + tau_u_g(i,j,k) + &
                  ((v0)*u_go(i,j,k) + vbf)*vol)

            enddo
         enddo
      enddo

   end subroutine source_u_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOURCE_U_g_BC                                           !
!                                                                      !
!  Purpose: Determine source terms for U_g momentum eq. The terms      !
!     appear in the center coefficient and RHS vector. The center      !
!     coefficient and source vector are negative. The off-diagonal     !
!     coefficients are positive.                                       !
!     The drag terms are excluded from the source at this stage.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine source_u_g_bc(lo, hi, slo, shi, alo, ahi, A_m, b_m, &
      bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
      bc_klo_type, bc_khi_type, domlo, domhi, dy, dz)

      use bc, only: nsw_, fsw_, psw_
      use bc, only: pinf_, pout_
      use bc, only: minf_

      use bc, only: bc_hw_g, bc_uw_g, bc_u_g

      use matrix, only: e, w, s, n, t, b
      use param, only: is_defined

      integer     , intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3)
      real(c_real), intent(in   ) :: dy, dz

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

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

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      real(c_real) :: ody, odz
      integer      :: bcv,ibc,i,j,k
!-----------------------------------------------

      ody = 1.d0 / dy
      odz = 1.d0 / dz

! ------------------------------------------------------------->

      ! At left boundary
      if (slo(1) .lt. domlo(1) .and. lo(1).eq.alo(1)) then

         i = alo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)

               if (bc_ilo_type(j,k,1) == PINF_ .or. &
                   bc_ilo_type(j,k,1) == POUT_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,w)
                  A_m(i,j,k,w) = zero

               else if (bc_ilo_type(j,k,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one

                  bcv = bc_ilo_type(j,k,2)
                  b_m(i,j,k) = -bc_u_g(bcv)

               else if(bc_ilo_type(j,k,1) == NSW_ .or. &
                       bc_ilo_type(j,k,1) == FSW_ .or. &
                       bc_ilo_type(j,k,1) == PSW_) then

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

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)

               if(bc_ihi_type(j,k,1) == PINF_ .or. &
                  bc_ihi_type(j,k,1) == POUT_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,e)
                  A_m(i,j,k,e) = zero

               else if(bc_ihi_type(j,k,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one

                  bcv = bc_ihi_type(j,k,2)
                  b_m(i,j,k) = -bc_u_g(bcv)

               else if(bc_ihi_type(j,k,1) == NSW_ .or. &
                       bc_ihi_type(j,k,1) == FSW_ .or. &
                       bc_ihi_type(j,k,1) == PSW_) then

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
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)

               ! bc's on j-faces only defined within (domlo(1):domhi(1),domlo(3):domhi(3))
               ibc = max(min(i,domhi(1)),domlo(1))

               if(bc_jlo_type(i,k,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,s)
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(ibc,k,1) == FSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,s)
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(ibc,k,1) == PSW_) then

                  bcv = bc_jlo_type(ibc,k,2)
                  if (is_undefined(bc_hw_g(bcv))) then
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,s)
                     b_m(i,j,k) = b_m(i,j,k) - 2.0*A_m(i,j,k,s)*bc_uw_g(bcv)
                  else
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,s)*&
                        (half*bc_hw_g(bcv)-ody)/(half*bc_hw_g(bcv)+ody)
                     b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,s)*&
                        bc_hw_g(bcv)*bc_uw_g(bcv)/(half*bc_hw_g(bcv)+ody)
                  endif
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(ibc,k,1) == PINF_ .or. &
                       bc_jlo_type(ibc,k,1) == POUT_ .or. &
                       bc_jlo_type(ibc,k,1) == MINF_) then

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
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)

               ! bc's on j-faces only defined within (domlo(1):domhi(1),domlo(3):domhi(3))
               ibc = max(min(i,domhi(1)),domlo(1))

               bcv = bc_jhi_type(ibc,k,2)

               if(bc_jhi_type(ibc,k,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,n)
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(ibc,k,1) == FSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,n)
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(ibc,k,1) == PSW_) then
                  if (is_undefined(bc_hw_g(bcv))) then
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,n)
                     b_m(i,j,k) = b_m(i,j,k) - 2.0*A_m(i,j,k,n)*bc_uw_g(bcv)
                  else
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,n)*&
                        (half*bc_hw_g(bcv)-ody)/(half*bc_hw_g(bcv)+ody)
                     b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,n)*&
                        bc_hw_g(bcv)*bc_uw_g(bcv)/(half*bc_hw_g(bcv)+ody)
                  endif
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(ibc,k,1) == PINF_ .or. &
                       bc_jhi_type(ibc,k,1) == POUT_ .or. &
                       bc_jhi_type(ibc,k,1) == MINF_) then

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
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! bc's on k-faces only defined within (domlo(1):domhi(1),domlo(2):domhi(2))
               ibc = max(min(i,domhi(1)),domlo(1))

               bcv = bc_klo_type(ibc,j,2)

               if(bc_klo_type(ibc,j,1) == NSW_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,b)
                  A_m(i,j,k,b) = zero

               else if(bc_klo_type(ibc,j,1) == FSW_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,b)
                  A_m(i,j,k,b) = zero

               else if(bc_klo_type(ibc,j,1) == PSW_) then

                  if (is_undefined(bc_hw_g(bcv))) then
                     A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,b)
                     b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,b)*bc_uw_g(bcv)
                  else
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,b)*&
                        (half*bc_hw_g(bcv)-odz)/(half*bc_hw_g(bcv)+odz)
                     b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,b)*&
                        bc_hw_g(bcv)*bc_uw_g(bcv)/(half*bc_hw_g(bcv)+odz)
                  endif
                  A_m(i,j,k,b) = zero

               else if(bc_klo_type(ibc,j,1) == PINF_ .or. &
                       bc_klo_type(ibc,j,1) == POUT_ .or. &
                       bc_klo_type(ibc,j,1) == MINF_) then

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
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! bc's on k-faces only defined within (domlo(1):domhi(1),domlo(2):domhi(2))
               ibc = max(min(i,domhi(1)),domlo(1))

               bcv = bc_khi_type(ibc,j,2)

               if(bc_khi_type(ibc,j,1) == NSW_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,t)
                  A_m(i,j,k,t) = zero

               else if(bc_khi_type(ibc,j,1) == FSW_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,t)
                  A_m(i,j,k,t) = zero

               else if(bc_khi_type(ibc,j,1) == PSW_) then

                  if (is_undefined(bc_hw_g(bcv))) then
                     A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,t)
                     b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,t)*bc_uw_g(bcv)
                  else
                     A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,t)*&
                        (half*bc_hw_g(bcv)-odz)/(half*bc_hw_g(bcv)+odz)
                     b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,t)*&
                        bc_hw_g(bcv)*bc_uw_g(bcv)/(half*bc_hw_g(bcv)+odz)
                  endif
                  A_m(i,j,k,t) = zero

               else if(bc_khi_type(ibc,j,1) == PINF_ .or. &
                       bc_khi_type(ibc,j,1) == POUT_ .or. &
                       bc_khi_type(ibc,j,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero

               endif
            end do
         end do
      endif

   end subroutine source_u_g_bc

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_U_G                                        C
!  Purpose: Adds point sources to the gas phase U-momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine point_source_u_g(lo, hi, alo, ahi, b_m, vol)

      ! use ps, only: dim_ps, ps_defined, ps_vel_mag_g, ps_massflow_g
      ! use ps, only: ps_u_g

      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      real(c_real)  , intent(inout) :: b_m(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real)  , intent(in   ) :: vol

      ! integer(c_int) :: I, J, K
      ! integer(c_int) :: psv
      ! integer(c_int) :: lIE, lIW
      ! real(c_real)  :: pSource
!-----------------------------------------------

      ! Calculate the mass going into each (i,j,k) cell. This is done for each
      ! call in case the point source is time dependent.

      ! ps_lp: do psv = 1, dim_ps
      !    if (.not.ps_defined(psv)) cycle PS_LP
      !    if (abs(PS_U_g(psv)) < small_number) cycle PS_LP

      !    if(PS_U_g(psv) < 0.0d0) then
      !       lIW = PS_I_W(psv) - 1
      !       lIE = PS_I_E(psv) - 1
      !    else
      !       lIW = PS_I_W(psv)
      !       lIE = PS_I_E(psv)
      !    endif

      !    do k = PS_K_B(psv), PS_K_T(psv)
      !    do j = PS_J_S(psv), PS_J_N(psv)
      !    do i = lIW, lIE

      !       pSource =  PS_MASSFLOW_G(psv) * (vol/PS_VOLUME(psv))

      !       b_m(I,J,K) = b_m(I,J,K) - pSource *                        &
      !          PS_U_g(psv) * PS_VEL_MAG_g(psv)

      !    enddo
      !    enddo
      !    enddo

      ! enddo ps_lp

      end subroutine point_source_u_g
end module source_u_g_module
