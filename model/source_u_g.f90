module source_u_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use param1        , only: zero, half, one, undefined, is_undefined, small_number

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
   subroutine source_u_g(slo, shi, ulo, uhi, alo, ahi, lo, hi, A_m, b_m, &
      dt, p_g, ep_g, ro_g, rop_g, rop_go, u_go, &
      tau_u_g, dx, dy, dz)


! Modules
!---------------------------------------------------------------------//
      USE constant, only: gravity
      USE bc      , only: delp_x

      USE functions, only: avg

      USE geometry, only: domlo, domhi, cyclic_x_pd

      use matrix, only: e, w, s, n, t, b
      USE scales, only: p_scale

      integer     , intent(in   ) :: slo(3),shi(3),ulo(3),uhi(3),alo(3),ahi(3),lo(3),hi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), intent(INOUT) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      ! Vector b_m
      real(c_real), intent(INOUT) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_g&
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
      real(c_real) :: PgE
! Average volume fraction
      real(c_real) :: EPGA
! Average density
      real(c_real) :: ROPGA, ROGA
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
            do i = lo(1)-1, hi(1)

               epga = avg(ep_g(i,j,k),ep_g(i+1,j,k))

! Pressure term
               PGE = P_G(i+1,j,k)
               if(cyclic_x_pd) then
                  if ((i == domlo(1)-1) .or. (i == domhi(1))) &
                     pge = pge - delp_x
               end if
               sdp = -p_scale*epga*(pge - p_g(i,j,k))*ayz

! Volumetric forces
               roga  = half * (ro_g(i,j,k) + ro_g(i+1,j,k))
               ropga = half * (rop_g(i,j,k) + rop_g(i+1,j,k))

! Previous time step
               v0 = half * (rop_go(i,j,k) + rop_go(i+1,j,k))*odt

! Body force
               vbf = roga*gravity(1)

! Collect the terms
               A_m(i,j,k,0) = -(A_m(i,j,k,e) + A_m(i,j,k,w) + &
                                A_m(i,j,k,n) + A_m(i,j,k,s) + &
                                A_m(i,j,k,t) + A_m(i,j,k,b) + &
                                v0*vol)

               b_m(i,j,k) = b_m(i,j,k) -(sdp + tau_u_g(i,j,k) + &
                  ((v0)*u_go(i,j,k) + vbf)*vol)

               ! if (i.eq. 3.and.j.eq.3.and.k.eq.0) print *,'B_M ', &
               ! b_m(i,j,k), sdp , tau_u_g(i,j,k) , v0,u_go(i,j,k) , vbf, vol

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
   subroutine source_u_g_bc(slo, shi, alo, ahi, ulo, uhi, A_m, b_m, u_g, &
      bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
      bc_klo_type, bc_khi_type, dy, dz)

      use ic, only: NSW_, FSW_, PSW_
      use ic, only: PINF_, POUT_
      use ic, only: MINF_, MOUT_
      use ic, only: CYCL_

      use bc, only: bc_hw_g, bc_uw_g, bc_u_g
      use geometry, only: domlo, domhi

      use matrix, only: e, w, s, n, t, b
      use param1, only: is_defined

      integer     , intent(in   ) :: slo(3),shi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      integer     , intent(in   ) :: ulo(3),uhi(3)
      real(c_real), intent(in   ) :: dy, dz

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

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

!-----------------------------------------------
! Local Variables
!-----------------------------------------------

      real(c_real) :: ody, odz

      integer :: bcv,i,j,k

      integer :: nlft, nrgt, nbot, ntop, nup, ndwn
!-----------------------------------------------

      ody = 1.d0 / dy
      odz = 1.d0 / dz

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

               if (bc_ilo_type(j,k,1) == PINF_ .or. &
                   bc_ilo_type(j,k,1) == POUT_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,w)
                  A_m(i,j,k,w) = zero

               else if (bc_ilo_type(j,k,1) == MINF_ .or. &
                        bc_ilo_type(j,k,1) == MOUT_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
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

! --- WEST FLUID ---------------------------------------------------------->

      if (nrgt .gt. 0) then
         i = ahi(1)
         do k=alo(3),ahi(3)
            do j=alo(2),ahi(2)
               bcv = bc_ihi_type(j,k,2)

               if(bc_ihi_type(j,k,1) == PINF_ .or. &
                  bc_ihi_type(j,k,1) == POUT_) then

                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,e)
                  A_m(i,j,k,e) = zero

               else if(bc_ihi_type(j,k,1) == MINF_ .or. &
                       bc_ihi_type(j,k,1) == MOUT_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
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
                       bc_jlo_type(i,k,1) == CYCL_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,s)
                  A_m(i,j,k,s) = zero

               else if(bc_jlo_type(i,k,1) == PSW_) then
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

               else if(bc_jlo_type(i,k,1) == PINF_ .or. &
                       bc_jlo_type(i,k,1) == POUT_ .or. &
                       bc_jlo_type(i,k,1) == MINF_ .or. &
                       bc_jlo_type(i,k,1) == MOUT_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

               ! b_m(i,j-1,k) = zero
               ! A_m(i,j-1,k,:) = zero
               ! A_m(i,j-1,k,0) = -one
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
                       bc_jhi_type(i,k,1) == CYCL_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,n)
                  A_m(i,j,k,n) = zero

               else if(bc_jhi_type(i,k,1) == PSW_) then
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

               else if(bc_jhi_type(i,k,1) == PINF_ .or. &
                       bc_jhi_type(i,k,1) == POUT_ .or. &
                       bc_jhi_type(i,k,1) == MINF_ .or. &
                       bc_jhi_type(i,k,1) == MOUT_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k) = zero
               endif

               ! b_m(i,j+1,k) = zero
               ! A_m(i,j+1,k,:) = zero
               ! A_m(i,j+1,k,0) = -one
            end do
         end do
      endif

! --- TOP FLUID ----------------------------------------------------------->

      if (ndwn .gt. 0) then
         k = alo(3)
         do j=alo(2),ahi(2)
            do i=alo(1),ahi(1)
               bcv = bc_klo_type(i,j,2)
               if(bc_klo_type(i,j,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,b)
                  A_m(i,j,k,b) = zero

               else if(bc_klo_type(i,j,1) == FSW_ .or. &
                       bc_klo_type(i,j,1) == CYCL_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,b)
                  A_m(i,j,k,b) = zero

               else if(bc_klo_type(i,j,1) == PSW_) then
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

               else if(bc_klo_type(i,j,1) == PINF_ .or. &
                       bc_klo_type(i,j,1) == POUT_ .or. &
                       bc_klo_type(i,j,1) == MINF_ .or. &
                       bc_klo_type(i,j,1) == MOUT_) then

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
               if(bc_khi_type(i,j,1) == NSW_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,t)
                  A_m(i,j,k,t) = zero

               else if(bc_khi_type(i,j,1) == FSW_ .or. &
                       bc_khi_type(i,j,1) == CYCL_) then
                  A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,t)
                  A_m(i,j,k,t) = zero

               else if(bc_khi_type(i,j,1) == PSW_) then
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

               else if(bc_khi_type(i,j,1) == PINF_ .or. &
                       bc_khi_type(i,j,1) == POUT_ .or. &
                       bc_khi_type(i,j,1) == MINF_ .or. &
                       bc_khi_type(i,j,1) == MOUT_) then

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
      subroutine point_source_u_g(alo, ahi, b_m, vol)

      use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      use ps, only: ps_u_g, ps_i_e, ps_i_w, ps_j_s, ps_j_n, ps_k_b, ps_k_t

      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      real(c_real)  , intent(inout) :: b_m(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real)  , intent(in   ) :: vol

      integer(c_int) :: I, J, K
      integer(c_int) :: psv
      integer(c_int) :: lIE, lIW
      real(c_real)  :: pSource
!-----------------------------------------------

      ! Calculate the mass going into each (i,j,k) cell. This is done for each
      ! call in case the point source is time dependent.

      ps_lp: do psv = 1, dimension_ps
         if (.not.ps_defined(psv)) cycle PS_LP
         if (abs(PS_U_g(psv)) < small_number) cycle PS_LP

         if(PS_U_g(psv) < 0.0d0) then
            lIW = PS_I_W(psv) - 1
            lIE = PS_I_E(psv) - 1
         else
            lIW = PS_I_W(psv)
            lIE = PS_I_E(psv)
         endif

         do k = PS_K_B(psv), PS_K_T(psv)
         do j = PS_J_S(psv), PS_J_N(psv)
         do i = lIW, lIE

            pSource =  PS_MASSFLOW_G(psv) * (vol/PS_VOLUME(psv))

            b_m(I,J,K) = b_m(I,J,K) - pSource *                        &
               PS_U_g(psv) * PS_VEL_MAG_g(psv)

         enddo
         enddo
         enddo

      enddo ps_lp

      end subroutine point_source_u_g
end module source_u_g_module
