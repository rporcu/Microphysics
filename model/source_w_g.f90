module source_w_g_module

   use bl_fort_module, only : c_real
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
        A_m, b_m, dt, p_g, ep_g, ro_g, rop_g, rop_go, w_go, &
        tau_w_g, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: gravity
      USE bc, only: delp_z

      USE functions, only: avg
      USE geometry , only: domlo, domhi, cyclic_z_pd

      use matrix, only: e, w, s, n, t, b
      USE scales, only: p_scale

      integer     , intent(in   ) :: slo(3),shi(3),wlo(3),whi(3),alo(3),ahi(3),lo(3),hi(3)

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
      real(c_real), intent(in   ) :: rop_g&
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
      INTEGER :: i,j,k
! Pressure at top cell
      real(c_real) :: PgT
! Average volume fraction
      real(c_real) :: EPGA
! Average density
      real(c_real) :: ROPGA, ROGA
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

      do k = lo(3)-1,hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               epga = avg(ep_g(i,j,k),ep_g(i,j,k+1))

! Pressure term
               pgt = p_g(i,j,k+1)
               if(cyclic_z_pd) then
                  if((k==domlo(3)-1) .or. (k==domhi(3)) ) &
                     pgt = p_g(i,j,k+1) - delp_z
               end if
               sdp = -p_scale*epga*(pgt - p_g(i,j,k))*axy

! Volumetric forces
               roga  = avg( ro_g(i,j,k), ro_g(i,j,k+1))
               ropga = avg(rop_g(i,j,k),rop_g(i,j,k+1))

! Previous time step
               v0 = avg(rop_go(i,j,k),rop_go(i,j,k+1))*odt

! Body force
               vbf = ropga*gravity(3)

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
         bc_klo_type, bc_khi_type, dx, dy)


      use ic, only: NSW_, FSW_, PSW_
      use ic, only: PINF_, POUT_
      use ic, only: MINF_, MOUT_
      use ic, only: CYCL_

      use bc, only: bc_hw_g, bc_ww_g, bc_w_g
      use geometry, only: domlo, domhi

      use matrix, only: e, w, s, n, t, b
      use param1, only: is_defined

      integer     , intent(in   ) :: slo(3),shi(3),alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx, dy

      real(c_real), intent(INOUT) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

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
                       bc_ilo_type(j,k,1) == MINF_ .or. &
                       bc_ilo_type(j,k,1) == MOUT_) then

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
                       bc_ihi_type(j,k,1) == MINF_ .or. &
                       bc_ihi_type(j,k,1) == MOUT_) then

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
                       bc_jlo_type(i,k,1) == MINF_ .or. &
                       bc_jlo_type(i,k,1) == MOUT_) then

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
                       bc_jhi_type(i,k,1) == CYCL_) then
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
                       bc_jhi_type(i,k,1) == MINF_ .or. &
                       bc_jhi_type(i,k,1) == MOUT_) then

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

               else if (bc_klo_type(i,j,1) == MINF_ .or. &
                        bc_klo_type(i,j,1) == MOUT_) then

                  A_m(i,j,k,:) =  zero
                  A_m(i,j,k,0) = -one
                  b_m(i,j,k  ) = bc_w_g(bcv)

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

               else if(bc_khi_type(i,j,1) == MINF_ .or. &
                       bc_khi_type(i,j,1) == MOUT_) then

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
      SUBROUTINE POINT_SOURCE_W_G(slo,shi,alo,ahi,b_m,flag,dx,dy,dz)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param1  , only: small_number, zero
      use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      use ps, only: ps_w_g, ps_i_e, ps_i_w, ps_j_s, ps_j_n, ps_k_b, ps_k_t

      integer     , intent(in   ) :: slo(3),shi(3),alo(3),ahi(3)

      ! Vector b_m
      real(c_real), intent(INOUT) :: B_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      integer, intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(IN   ) :: dx,dy,dz
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K
      INTEGER :: PSV
      INTEGER :: lKT, lKB
! terms of bm expression
      real(c_real) :: pSource
      real(c_real) :: vol
!-----------------------------------------------

      vol = dx*dy*dz

! Calculate the mass going into each (i,j,k) cell. This is done for each
! call in case the point source is time dependent.
      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(abs(PS_W_g(PSV)) < small_number) cycle PS_LP

         if(PS_W_g(PSV) < ZERO) then
            lKB = PS_K_B(PSV)-1
            lKT = PS_K_T(PSV)-1
         else
            lKB = PS_K_B(PSV)
            lKT = PS_K_T(PSV)
         endif

         do k = lKB, lKT
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.1.eq.flag(i,j,k,1)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

            b_m(I,J,K) = b_m(I,J,K) - pSource * &
               PS_W_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_W_G
end module source_w_g_module
