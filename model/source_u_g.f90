module source_u_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int
   use param1        , only: zero, half, one, undefined, is_undefined, small_number

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
   subroutine source_u_g(slo, shi, lo, hi, A_m, b_m, &
      dt, p_g, ep_g, ro_g, rop_g, rop_go, u_g, u_go, &
      tau_u_g, flag, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: gravity
      USE bc      , only: delp_x

      USE functions, only: avg, zmax, ieast

      USE geometry, only: domlo, domhi, cyclic_x_pd

      use matrix, only: e, w, s, n, t, b
      USE scales, only: p_scale

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      real(c_real), intent(in   ) :: dt, dx, dy, dz

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

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
      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: tau_u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer, intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: i,j,k
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
            do i = lo(1), hi(1)

               epga = avg(ep_g(i,j,k),ep_g(ieast(i,j,k),j,k))

! Pressure term
               PGE = P_G(ieast(i,j,k),j,k)
               if(cyclic_x_pd) then
                  if ((i == domlo(1)-1) .or. (i == domhi(1))) &
                     pge = pge - delp_x
               end if
               sdp = -p_scale*epga*(pge - p_g(i,j,k))*ayz

! Volumetric forces
               roga  = half * (ro_g(i,j,k) + ro_g(ieast(i,j,k),j,k))
               ropga = half * (rop_g(i,j,k) + rop_g(ieast(i,j,k),j,k))

! Previous time step
               v0 = half * (rop_go(i,j,k) + rop_go(ieast(i,j,k),j,k))*odt

! Body force
               vbf = roga*gravity(1)

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

      return
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
   subroutine source_u_g_bc(slo, shi, lo, hi, A_m, b_m, &
      u_g, flag, dx, dy, dz)

      use ic, only: NSW_, FSW_, PSW_
      use ic, only: PINF_, POUT_
      use ic, only: MINF_, MOUT_

      use bc, only: dimension_bc, bc_type, bc_defined, bc_plane
      use bc, only: bc_i_w, bc_i_e, bc_j_s, bc_j_n, bc_k_b, bc_k_t
      use bc, only: bc_hw_g, bc_uw_g, bc_u_g

      use matrix, only: e, w, s, n, t, b
      use param1, only: is_defined

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Velocity u_g
      real(c_real), INTENT(IN   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER ::  I,  J, K, IM, I1, I2, J1, J2, K1, K2

      real(c_real) :: ody, odz
!-----------------------------------------------

     ody = 1.d0 / dy
     odz = 1.d0 / dz

     do l = 1, dimension_bc
        if (bc_defined(l)) then

            i1 = bc_i_w(l)
            i2 = bc_i_e(l)
            j1 = bc_j_s(l)
            j2 = bc_j_n(l)
            k1 = bc_k_b(l)
            k2 = bc_k_t(l)

            if(i1 == i2) then
               if(i1 == slo(1) ) i1=i1+1
               if(i1 == shi(1) ) i1=i1-1
               i2=i1
            endif

            if(j1 == j2) then
               if(j1 == slo(2) ) j1=j1+1
               if(j1 == shi(2) ) j1=j1-1
               j2=j1
            endif

            if(k1 == k2)then
               if(k1 == slo(3) ) k1=k1+1
               if(k1 == shi(3) ) k1=k1-1
               k2=k1
            endif

            do k = k1, k2
               do j = j1, j2
                  do i = i1, i2

! --- EAST FLUID ---------------------------------------------------------->

! MASS INFLOW
                     if(is_defined(bc_u_g(l)) .and. (&
                        flag(i-1,j,k,1) == MINF_ .or. &
                        flag(i-1,j,k,1) == MOUT_)) then

                        A_m(i-1,j,k,:) =  zero
                        A_m(i-1,j,k,0) = -one
                        b_m(i-1,j,k) = -bc_u_g(l)
                     endif

! --- WEST FLUID ---------------------------------------------------------->

! PRESSURE IN/OUTFLOW
                     if(flag(i+1,j,k,1) == PINF_ .or. &
                        flag(i+1,j,k,1) == POUT_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,e)
                        A_m(i,j,k,e) = zero

                        b_m(i+1,j,k) = zero
                        A_m(i+1,j,k,:) = zero
                        A_m(i+1,j,k,0) = -one

! MASS INFLOW
                     else if(is_defined(bc_u_g(l)) .and. (&
                        flag(i+1,j,k,1) == MINF_ .or. &
                        flag(i+1,j,k,1) == MOUT_)) then
                        A_m(i,j,k,:) =  zero
                        A_m(i,j,k,0) = -one
                        b_m(i,j,k) = -bc_u_g(l)

                        b_m(i+1,j,k) = zero
                        A_m(i+1,j,k,:) = zero
                        A_m(i+1,j,k,0) = -one
                     endif

! --- NORTH FLUID --------------------------------------------------------->

! NO-SLIP WALL
                     if (flag(i,j-1,k,1) == NSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,s)
                        A_m(i,j,k,s) = zero

                        b_m(i,j-1,k) = zero
                        A_m(i,j-1,k,:) = zero
                        A_m(i,j-1,k,0) = -one

! FREE-SLIP WALL
                     else if (flag(i,j-1,k,1) == FSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,s)
                        A_m(i,j,k,s) = zero

                        b_m(i,j-1,k) = zero
                        A_m(i,j-1,k,:) = zero
                        A_m(i,j-1,k,0) = -one

! PARTIAL-SLIP WALL
                     else if (flag(i,j-1,k,1) == PSW_) THEN
                        if (is_undefined(bc_hw_g(l))) then
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,s)
                           b_m(i,j,k) = b_m(i,j,k) - 2.0*A_m(i,j,k,s)*bc_uw_g(l)
                        else
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,s)*&
                              (half*bc_hw_g(l)-ody)/(half*bc_hw_g(l)+ody)
                           b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,s)*&
                              bc_hw_g(l)*bc_uw_g(l)/(half*bc_hw_g(l)+ody)
                        endif
                        A_m(i,j,k,s) = zero

                        A_m(i,j-1,k,:) = zero
                        A_m(i,j-1,k,0) = -one
                        b_m(i,j-1,k) = zero
                     endif

! --- SOUTH FLUID --------------------------------------------------------->

! NO-SLIP WALL
                     if (flag(i,j+1,k,1) == NSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,n)
                        A_m(i,j,k,n) = zero

                        b_m(i,j+1,k) = zero
                        A_m(i,j+1,k,:) = zero
                        A_m(i,j+1,k,0) = -one

! FREE-SLIP WALL
                     else if (flag(i,j+1,k,1) == FSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,n)
                        A_m(i,j,k,n) = zero

                        b_m(i,j+1,k) = zero
                        A_m(i,j+1,k,:) = zero
                        A_m(i,j+1,k,0) = -one

! PARTIAL-SLIP WALL TO NORTH
                     else if (flag(i,j+1,k,1) == PSW_) then
                        if (is_undefined(bc_hw_g(l))) then
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,n)
                           b_m(i,j,k) = b_m(i,j,k) - 2.0*A_m(i,j,k,n)*bc_uw_g(l)
                        else
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,n)*&
                              (half*bc_hw_g(l)-ody)/(half*bc_hw_g(l)+ody)
                           b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,n)*&
                              bc_hw_g(l)*bc_uw_g(l)/(half*bc_hw_g(l)+ody)
                        endif
                        A_m(i,j,k,n) = zero

                        A_m(i,j+1,k,:) = zero
                        A_m(i,j+1,k,0) = -one
                        b_m(i,j+1,k) = zero
                     endif



! --- TOP FLUID ----------------------------------------------------------->

! NO-SLIP WALL
                     if (flag(i,j,k-1,1) == NSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,b)
                        A_m(i,j,k,b) = zero

                        b_m(i,j,k-1) = zero
                        A_m(i,j,k-1,:) = zero
                        A_m(i,j,k-1,0) = -one

! FREE-SLIP WALL
                     else if (flag(i,j,k-1,1) == FSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,b)
                        A_m(i,j,k,b) = zero

                        b_m(i,j,k-1) = zero
                        A_m(i,j,k-1,:) = zero
                        A_m(i,j,k-1,0) = -one

! PARTIAL-SLIP WALL
                     else if (flag(i,j,k-1,1) == PSW_) then
                        if (is_undefined(bc_hw_g(l))) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,b)
                           b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,b)*bc_uw_g(l)
                        else
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,b)*&
                              (half*bc_hw_g(l)-odz)/(half*bc_hw_g(l)+odz)
                           b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,b)*&
                              bc_hw_g(l)*bc_uw_g(l)/(half*bc_hw_g(l)+odz)
                        endif
                        A_m(i,j,k,b) = zero

                        A_m(i,j,k-1,:) = zero
                        A_m(i,j,k-1,0) = -one
                        b_m(i,j,k-1) = zero
                     endif


! --- BOTTOM FLUID -------------------------------------------------------->

! NO-SLIP WALL TO TOP
                     if (flag(i,j,k+1,1) == NSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,t)
                        A_m(i,j,k,t) = zero

                        b_m(i,j,k+1) = zero
                        A_m(i,j,k+1,:) = zero
                        A_m(i,j,k+1,0) = -one

! FREE-SLIP WALL
                     else if (flag(i,j,k+1,1) == FSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,t)
                        A_m(i,j,k,t) = zero

                        b_m(i,j,k+1) = zero
                        A_m(i,j,k+1,:) = zero
                        A_m(i,j,k+1,0) = -one

! PARTIAL-SLIP WALL
                     else if (flag(i,j,k+1,1) == PSW_) then
                        if (is_undefined(bc_hw_g(l))) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,t)
                           b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,t)*bc_uw_g(l)
                        else
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,t)*&
                              (half*bc_hw_g(l)-odz)/(half*bc_hw_g(l)+odz)
                           b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,t)*&
                              bc_hw_g(l)*bc_uw_g(l)/(half*bc_hw_g(l)+odz)
                        endif
                        A_m(i,j,k,t) = zero

                        A_m(i,j,k+1,:) = zero
                        A_m(i,j,k+1,0) = -one
                        b_m(i,j,k+1) = zero
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo

      END SUBROUTINE SOURCE_U_G_BC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_U_G                                        C
!  Purpose: Adds point sources to the gas phase U-momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_U_G(slo, shi, lo, hi, A_m, b_m, flag, dx, dy, dz)

      use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      use ps, only: ps_u_g, ps_i_e, ps_i_w, ps_j_s, ps_j_n, ps_k_b, ps_k_t

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(IN   ) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer, intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), INTENT(IN   ) :: dx,dy,dz

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K
      INTEGER :: PSV, M
      INTEGER :: lIE, lIW
! terms of bm expression
      real(c_real) :: pSource
      real(c_real) :: vol
!-----------------------------------------------

      vol = dx*dy*dz

! Set reference phase to gas
      M = 0

! Calculate the mass going into each (i,j,k) cell. This is done for each
! call in case the point source is time dependent.
      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(abs(PS_U_g(PSV)) < small_number) cycle PS_LP

         if(PS_U_g(PSV) < 0.0d0) then
            lIW = PS_I_W(PSV) - 1
            lIE = PS_I_E(PSV) - 1
         else
            lIW = PS_I_W(PSV)
            lIE = PS_I_E(PSV)
         endif

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = lIW, lIE

            if(.NOT. 1.eq.flag(i,j,k,1)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (vol/PS_VOLUME(PSV))

            b_m(I,J,K) = b_m(I,J,K) - pSource *                        &
               PS_U_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo
      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_U_G
end module source_u_g_module
