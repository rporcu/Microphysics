module source_w_g_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int
   use param1        , only: zero, half, one, undefined, is_undefined

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
   subroutine source_w_g(slo, shi, wlo, whi, alo, ahi, lo, hi, &
        A_m, b_m, dt, p_g, ep_g, ro_g, rop_go, w_go, &
        tau_w_g, dx, dy, dz, domlo, domhi)

! Modules
!---------------------------------------------------------------------//
      use constant, only: gravity
      use bc, only: delp_z

      use functions, only: avg
      use bc , only: cyclic_z_pd

      use matrix, only: e, w, s, n, t, b
      use scales, only: p_scale

      integer     , intent(in   ) :: slo(3),shi(3),wlo(3),whi(3),alo(3),ahi(3),lo(3),hi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3)

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
      real(c_real), intent(in   ) :: w_go&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: tau_w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: dt, dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      integer :: i,j,k
! Pressure at top cell
      real(c_real) :: PgB
! Average volume fraction
      real(c_real) :: EPGA
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

      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               epga = half*(ep_g(i,j,k-1) + ep_g(i,j,k))

               ! Pressure term
               pgb = p_g(i,j,k-1)
               if(cyclic_z_pd) then
                  if((k==domlo(3)) .or. (k==domhi(3)+1) ) &
                     pgb = pgb + delp_z
               end if
               sdp = -p_scale*epga*(p_g(i,j,k) - pgb)*axy

               ! Previous time step
               v0 = half*(rop_go(i,j,k-1) + rop_go(i,j,k))*odt

               ! Body force
               vbf = half*(ro_g(i,j,k-1) + ro_g(i,j,k))*gravity(3)

               ! Collect the terms
               A_m(i,j,k,0) = -(A_m(i,j,k,e) + A_m(i,j,k,w) + &
                                A_m(i,j,k,n) + A_m(i,j,k,s) + &
                                A_m(i,j,k,t) + A_m(i,j,k,b) + v0*vol)

               b_m(i,j,k) = b_m(i,j,k) - ( sdp + tau_w_g(i,j,k)  + &
                  ( (v0)*w_go(i,j,k) + vbf)*vol)

            enddo
         enddo
      enddo

      return
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
      subroutine source_w_g_bc(slo,shi,alo,ahi,A_m,b_m, &
         bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
         bc_klo_type, bc_khi_type, domlo, domhi, dx, dy)


      use bc, only: nsw_, fsw_, psw_
      use bc, only: pinf_, pout_
      use bc, only: minf_
      use bc, only: cycl_

      use bc, only: bc_hw_g, bc_ww_g, bc_w_g

      use matrix, only: e, w, s, n, t, b
      use param1, only: is_defined

      integer     , intent(in   ) :: slo(3),shi(3),alo(3),ahi(3)
      integer     , intent(in   ) :: domlo(3),domhi(3)
      real(c_real), intent(in   ) :: dx, dy

      real(c_real), intent(INOUT) :: A_m&
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

! Local variables
!-----------------------------------------------
      real(c_real) :: odx, ody

      integer :: bcv, i,j,k

      integer :: nlft, nrgt, nbot, ntop, nup, ndwn
!---------------------------------------------------------------------//

      odx = 1.0d0/dx
      ody = 1.0d0/dy

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

! --- EAST FLUID ---------------------------------------------------------->

      if (nlft .gt. 0) then
         i = alo(1)
         do k=alo(3),ahi(3)
            do j=alo(2),ahi(2)
               bcv = bc_ilo_type(j,k,2)

               if(bc_ilo_type(j,k,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,w)
                  A_m(i,j,k,w) = zero

               else if(bc_ilo_type(j,k,1) == FSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,w)
                  A_m(i,j,k,w) = zero

               else if(bc_ilo_type(j,k,1) == PSW_) then
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

               else if(bc_ilo_type(j,k,1) == PINF_ .or. &
                       bc_ilo_type(j,k,1) == POUT_ .or. &
                       bc_ilo_type(j,k,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif
            end do
         end do
      endif

! --- WEST FLUID ---------------------------------------------------------->

      if (nrgt .gt. 0) then
         i = ahi(1)
         do k=alo(3),ahi(3)
            do j=alo(2),ahi(2)
               bcv = bc_ihi_type(j,k,2)

               if(bc_ihi_type(j,k,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,e)
                  A_m(i,j,k,e) = zero

               else if(bc_ihi_type(j,k,1) == FSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,e)
                  A_m(i,j,k,e) = zero

               else if(bc_ihi_type(j,k,1) == PSW_) then
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

               else if(bc_ihi_type(j,k,1) == PINF_ .or. &
                       bc_ihi_type(j,k,1) == POUT_ .or. &
                       bc_ihi_type(j,k,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

            end do
         end do
      endif

! --- NORTH FLUID --------------------------------------------------------->

      if (nbot .gt. 0) then
         j = alo(2)
         do k=alo(3),ahi(3)
            do i=alo(1),ahi(1)
               bcv = bc_jlo_type(i,k,2)
               if(bc_jlo_type(i,k,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,s)
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(i,k,1) == FSW_ .or. &
                       bc_jlo_type(i,k,1) == cycl_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,s)
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(i,k,1) == PSW_) then
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

               else if(bc_jlo_type(i,k,1) == PINF_ .or. &
                       bc_jlo_type(i,k,1) == POUT_ .or. &
                       bc_jlo_type(i,k,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

            end do
         end do
      endif

! --- SOUTH FLUID --------------------------------------------------------->

      if (ntop .gt. 0) then
         j = ahi(2)
         do k=alo(3),ahi(3)
            do i=alo(1),ahi(1)
               bcv = bc_jhi_type(i,k,2)

               if(bc_jhi_type(i,k,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,n)
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(i,k,1) == FSW_ .or. &
                       bc_jhi_type(i,k,1) == cycl_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,n)
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(i,k,1) == PSW_) then
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

               else if(bc_jhi_type(i,k,1) == PINF_ .or. &
                       bc_jhi_type(i,k,1) == POUT_ .or. &
                       bc_jhi_type(i,k,1) == MINF_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

            end do
         end do
      endif

! --- TOP FLUID ----------------------------------------------------------->

      if (ndwn .gt. 0) then
         k = alo(3)
         do j=alo(2),ahi(2)
            do i=alo(1),ahi(1)
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

! --- BOTTOM FLUID -------------------------------------------------------->

      if (nup .gt. 0) then
         k = ahi(3)
         do j=alo(2),ahi(2)
            do i=alo(1),ahi(1)
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
      subroutine point_source_w_g(alo,ahi,b_m,vol)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      ! use param1  , only: small_number, zero
      ! use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      ! use ps, only: ps_w_g

      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      real(c_real)  , intent(inout) :: b_m(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real)  , intent(in   ) :: vol

      ! integer(c_int) :: I, J, K
      ! integer(c_int) :: psv
      ! integer(c_int) :: lKT, lKB
      ! real(c_real)  :: pSource

      ! ! Calculate the mass going into each (i,j,k) cell. This is done for each
      ! ! call in case the point source is time dependent.
      ! ps_lp: do psv = 1, dimension_ps
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
