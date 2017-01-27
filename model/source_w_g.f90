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
   subroutine source_w_g(slo, shi, lo, hi, A_m, b_m, &
        dt, p_g, ep_g, ro_g, rop_g, rop_go, w_g, w_go, &
        tau_w_g, flag, dx, dy, dz)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: gravity
      USE bc, only: delp_z

      USE functions, only: avg
      USE geometry , only: domlo, domhi, cyclic_z_pd

      use matrix, only: e, w, s, n, t, b
      USE scales, only: p_scale

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), intent(inout) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      ! Vector b_m
      real(c_real), intent(inout) :: b_m&
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
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: tau_w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer, intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
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

      do k = lo(3),hi(3)
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
      subroutine source_w_g_bc(slo, shi, lo, hi, A_m, b_m, &
         w_g, flag, dx, dy, dz)

      use ic, only: NSW_, FSW_, PSW_
      use ic, only: PINF_, POUT_
      use ic, only: MINF_, MOUT_

      use bc, only: dimension_bc, bc_defined, bc_type
      use bc, only: bc_i_w, bc_i_e, bc_j_s, bc_j_n, bc_k_b, bc_k_t
      use bc, only: bc_hw_g, bc_ww_g, bc_w_g

      use matrix, only: e, w, s, n, t, b
      use param1, only: is_defined

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
! Vector b_m
      real(c_real), intent(inout) :: b_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer,      intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in) :: dx, dy, dz

! Local variables
!-----------------------------------------------
! Boundary condition
      integer :: l
! Indices
      integer :: i, j, k, i1, i2, j1, j2, k1, k2, im, jm
      real(c_real) :: odx, ody
!---------------------------------------------------------------------//

      odx = 1.0d0/dx
      ody = 1.0d0/dy

!---------------------------------------------------------------------//

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


            if (bc_type(l) == 'NO_SLIP_WALL') then

               do k = k1, k2
                  do j = j1, j2
                     do i = i1, i2

! East Fluid
                        if (flag(i-1,j,k,1) == NSW_) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,w)
                           A_m(i,j,k,w) = zero

                           b_m(i-1,j,k) = zero
                           A_m(i-1,j,k,:) = zero
                           A_m(i-1,j,k,0) = -one
                        endif
! West fluid
                        if (flag(i+1,j,k,1) == NSW_) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,e)
                           A_m(i,j,k,e) = zero

                           b_m(i+1,j,k) = zero
                           A_m(i+1,j,k,:) = zero
                           A_m(i+1,j,k,0) = -one
                        endif
! North fluid
                        if (flag(i,j-1,k,1) == NSW_) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,s)
                           A_m(i,j,k,s) = zero

                           b_m(i,j-1,k) = zero
                           A_m(i,j-1,k,:) = zero
                           A_m(i,j-1,k,0) = -one
                        endif
! South fluid
                        if (flag(i,j+1,k,1) == NSW_) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,n)
                           A_m(i,j,k,n) = zero

                           b_m(i,j+1,k) = zero
                           A_m(i,j+1,k,:) = zero
                           A_m(i,j+1,k,0) = -one
                        endif
                     enddo
                  enddo
               enddo



            else if (bc_type(l) == 'FREE_SLIP_WALL') then

               do k = k1, k2
                  do j = j1, j2
                     do i = i1, i2

                        if (flag(i-1,j,k,1) == FSW_) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,w)
                           A_m(i,j,k,w) = zero

                           b_m(i-1,j,k) = zero
                           A_m(i-1,j,k,:) = zero
                           A_m(i-1,j,k,0) = -one
                        endif
                        if (flag(i+1,j,k,1) == FSW_) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,e)
                           A_m(i,j,k,e) = zero

                           b_m(i+1,j,k) = zero
                           A_m(i+1,j,k,:) = zero
                           A_m(i+1,j,k,0) = -one
                        endif
                        if (flag(i,j-1,k,1) == FSW_) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,s)
                           A_m(i,j,k,s) = zero

                           b_m(i,j-1,k) = zero
                           A_m(i,j-1,k,:) = zero
                           A_m(i,j-1,k,0) = -one
                        endif
                        if (flag(i,j+1,k,1) == FSW_) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,n)
                           A_m(i,j,k,n) = zero

                           b_m(i,j+1,k) = zero
                           A_m(i,j+1,k,:) = zero
                           A_m(i,j+1,k,0) = -one
                        endif
                     enddo
                  enddo
               enddo

            elseif (bc_type(l) == 'PAR_SLIP_WALL') then

               do k = k1, k2
                  do j = j1, j2
                     do i = i1, i2

                        if (flag(i-1,j,k,1) == PSW_) then
                           if (is_undefined(bc_hw_g(l))) then
                              A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,w)
                              b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,w)*bc_ww_g(l)
                           else
                              A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,w)*&
                                 (half*bc_hw_g(l)-odx)/(half*bc_hw_g(l)+odx)
                              b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,w)*&
                                 bc_hw_g(l)*bc_ww_g(l)/(half*bc_hw_g(l)+odx)
                           endif
                           A_m(i,j,k,w) = zero

                           A_m(i,j,k-1,:) = zero
                           A_m(i,j,k-1,0) = -one
                           b_m(i,j,k-1) = zero
                        endif

                        if (flag(i+1,j,k,1) == PSW_) then
                           if (is_undefined(bc_hw_g(l))) then
                              A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,e)
                              b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,e)*bc_ww_g(l)
                           else
                              A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,e)*&
                                 (half*bc_hw_g(l)-odx)/(half*bc_hw_g(l)+odx)
                              b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,e)*&
                                 bc_hw_g(l)*bc_ww_g(l)/(half*bc_hw_g(l)+odx)
                           endif
                           A_m(i,j,k,e) = zero

                           A_m(i,j,k-1,:) = zero
                           A_m(i,j,k-1,0) = -one
                           b_m(i,j,k-1) = zero
                        endif
                        if (flag(i,j-1,k,1) == PSW_) THEN
                           if (is_undefined(bc_hw_g(l))) then
                              A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,s)
                              b_m(i,j,k) = b_m(i,j,k) - 2.0*A_m(i,j,k,s)*bc_ww_g(l)
                           else
                              A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,s)*&
                                 (half*bc_hw_g(l)-ody)/(half*bc_hw_g(l)+ody)
                              b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,s)*&
                                 bc_hw_g(l)*bc_ww_g(l)/(half*bc_hw_g(l)+ody)
                           endif
                           A_m(i,j,k,s) = zero

                           A_m(i,j-1,k,:) = zero
                           A_m(i,j-1,k,0) = -one
                           b_m(i,j-1,k) = zero
                        endif
                        if (flag(i,j+1,k,1) == PSW_) then
                           if (is_undefined(bc_hw_g(l))) then
                              A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,n)
                              b_m(i,j,k) = b_m(i,j,k) - 2.0*A_m(i,j,k,n)*bc_ww_g(l)
                           else
                              A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,n)*&
                                 (half*bc_hw_g(l)-ody)/(half*bc_hw_g(l)+ody)
                              b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,n)*&
                                 bc_hw_g(l)*bc_ww_g(l)/(half*bc_hw_g(l)+ody)
                           endif
                           A_m(i,j,k,n) = zero

                           A_m(i,j+1,k,:) = zero
                           A_m(i,j+1,k,0) = -one
                           b_m(i,j+1,k) = zero
                        endif

                     enddo
                  enddo
               enddo

            elseif (bc_type(l) == 'P_INFLOW' .or. &
                    bc_type(l) == 'P_OUTFLOW') then

               do k = k1, k2
                  do j = j1, j2
                     do i = i1, i2

                        if (flag(i+1,j,k,1) == PINF_ .or. &
                           flag(i+1,j,k,1) == POUT_ ) then

                           A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,e)
                           A_m(i,j,k,e) = zero

                           b_m(i+1,j,k) = zero
                           A_m(i+1,j,k,:) = zero
                           A_m(i+1,j,k,0) = -one
                        endif

                     enddo
                  enddo
               enddo

            else ! Mass inflow and Mass Outflow

               do k = k1, k2
                  do j = j1, j2
                     do i = i1, i2

                        if(flag(i,j,k-1,1) == MINF_ .or. &
                           flag(i,j,k-1,1) == MOUT_) then

                           A_m(i,j,k-1,:) =  zero
                           A_m(i,j,k-1,0) = -one
                           b_m(i,j,k-1) = -bc_w_g(l)
                        endif

                        if(flag(i,j,k+1,1) == MINF_ .or. &
                           flag(i,j,k+1,1) == MOUT_) then
                           A_m(i,j,k,:) =  zero
                           A_m(i,j,k,0) = -one
                           b_m(i,j,k) = -bc_w_g(l)

                           b_m(i,j,k+1) = zero
                           A_m(i,j,k+1,:) = zero
                           A_m(i,j,k+1,0) = -one
                        endif

                     enddo
                  enddo
               enddo
            endif
         endif
      enddo

      RETURN
      END SUBROUTINE SOURCE_W_G_BC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_W_G                                        C
!  Purpose: Adds point sources to the gas phase W-Momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_W_G(slo, shi, lo, hi, B_M, flag, dx, dy, dz)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param1  , only: small_number, zero
      use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      use ps, only: ps_w_g, ps_i_e, ps_i_w, ps_j_s, ps_j_n, ps_k_b, ps_k_t

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer, intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), INTENT(IN   ) :: dx,dy,dz
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

            B_M(I,J,K) = B_M(I,J,K) - pSource * &
               PS_W_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_W_G
end module source_w_g_module
