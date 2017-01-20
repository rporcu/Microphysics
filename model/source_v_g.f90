module source_v_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   use param1, only: zero, half, one, undefined, is_undefined, small_number

  contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_V_g                                              C
!  Purpose: Determine source terms for V_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector. The center         C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this stage.          C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_V_G(slo, shi, lo, hi, &
                            A_M, B_M, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
                            v_g, v_go, tau_v_g, flag, dx, dy, dz)


! Modules
!---------------------------------------------------------------------//
      USE constant, only: gravity
      USE bc, only: delp_y

      USE functions, only: avg
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus, jnorth
      USE functions, only: jnorth, jsouth
      USE functions, only: zmax
      USE geometry,  only: domlo, domhi, cyclic_y_pd

      use matrix, only: e, w, s, n, t, b

      USE scales, only: p_scale
      USE toleranc, only: dil_ep_s

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      ! Vector b_m
      real(c_real), INTENT(INOUT) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), INTENT(IN   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: rop_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: v_go&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: tau_v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dt, dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: i,j,k
! Pressure at north cell
      real(c_real) :: PgN
! Average volume fraction
      real(c_real) :: EPGA
! Average density
      real(c_real) :: ROPGA, ROGA
! Source terms (Surface)
      real(c_real) :: Sdp
! Source terms (Volumetric)
      real(c_real) :: V0, Vbf
! jackson terms: local stress tensor quantity
      real(c_real) :: ltau_v_g
      real(c_real) :: odt
      real(c_real) :: axz, vol
!---------------------------------------------------------------------//

      odt = 1.0d0/dt
      axz = dx*dz
      vol = dx*dy*dz

      DO K = lo(3), hi(3)
        DO J = lo(2), hi(2)+1
          DO I = lo(1), hi(1)

           EPGA = AVG(EP_G(I,J,K),EP_G(i,jnorth(i,j,k),k))

! Impermeable internal surface
         IF (flag(i,j,k,3)<1000) THEN
            A_M(I,J,K,E) = ZERO
            A_M(I,J,K,W) = ZERO
            A_M(I,J,K,N) = ZERO
            A_M(I,J,K,S) = ZERO
            A_M(I,J,K,T) = ZERO
            A_M(I,J,K,B) = ZERO
            A_M(I,J,K,0) = -ONE
            B_M(I,J,K) = ZERO

! dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_M(I,J,K,E) = ZERO
            A_M(I,J,K,W) = ZERO
            A_M(I,J,K,N) = ZERO
            A_M(I,J,K,S) = ZERO
            A_M(I,J,K,T) = ZERO
            A_M(I,J,K,B) = ZERO
            A_M(I,J,K,0) = -ONE
            B_M(I,J,K) = ZERO
            IF (EP_G(i,jsouth(i,j,k),k) > DIL_EP_S) THEN
               A_M(I,J,K,S) = ONE
            ELSE IF (EP_G(i,jnorth(i,j,k),k) > DIL_EP_S) THEN
               A_M(I,J,K,N) = ONE
            ELSE
               B_M(I,J,K) = -V_G(I,J,K)
            ENDIF

! Normal case
         ELSE

! Surface forces
! Pressure term
            PGN = P_G(i,jnorth(i,j,k),k)
            if ( CYCLIC_Y_PD) then
              if ( (j .eq. domlo(2)-1) .or. (j .eq. domhi(2)) ) &
               PGN = P_G(i,jnorth(i,j,k),k) - DELP_Y
            end if

            SDP = -P_SCALE*EPGA*(PGN - P_G(I,J,K))*AXZ

! Volumetric forces
            ROGA = AVG(RO_G(I,J,K),RO_G(i,jnorth(i,j,k),k))
            ROPGA = AVG(ROP_G(I,J,K),ROP_G(i,jnorth(i,j,k),k))
! Previous time step
            V0 = AVG(ROP_GO(I,J,K),ROP_GO(i,jnorth(i,j,k),k))*ODT

! Body force
            VBF = ROGA*GRAVITY(2)

! if jackson, implement jackson form of governing equations (ep_g dot
! del tau_g): multiply by void fraction otherwise by 1
            ltau_v_g = tau_v_g(i,j,k)


! Collect the terms
            A_M(I,J,K,0) = -(A_M(I,J,K,E)+A_M(I,J,K,W)+&
               A_M(I,J,K,N)+A_M(I,J,K,S)+A_M(I,J,K,T)+A_M(I,J,K,B)+&
               V0*VOL)
            B_M(I,J,K) = B_M(I,J,K) - (SDP + lTAU_V_G +  &
               ((V0)*V_GO(I,J,K) + VBF)*VOL )

         ENDIF
      ENDDO
      ENDDO
      ENDDO

! modifications for bc
      CALL SOURCE_V_G_BC(slo,shi,lo,hi,A_M,B_M,V_G,flag,dx,dy,dz)

      RETURN
      END SUBROUTINE SOURCE_V_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_V_g_BC                                           C
!  Purpose: Determine source terms for V_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_V_G_BC(slo,shi,lo,hi,A_M,B_M,v_g,flag,dx,dy,dz)

      use ic, only: NSW_, FSW_, PSW_
      use ic, only: PINF_, POUT_
      use ic, only: MINF_, MOUT_

      use bc, only: dimension_bc, bc_defined, bc_type, bc_plane
      use bc, only: bc_i_w, bc_i_e, bc_j_s, bc_j_n, bc_k_b, bc_k_t
      use bc, only: bc_hw_g, bc_vw_g, bc_v_g

      use matrix, only: e, w, s, n, t, b
      use param1, only: is_defined

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Velocity v_g
      real(c_real), INTENT(IN   ) :: v_g&
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
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2

      real(c_real) :: odx, odz
!-----------------------------------------------
     odx = 1.d0 / dx
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

! NO-SLIP WALL
                     if (flag(i-1,j,k,1) == NSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,w)
                        A_m(i,j,k,w) = zero

                        b_m(i-1,j,k) = zero
                        A_m(i-1,j,k,:) = zero
                        A_m(i-1,j,k,0) = -one

! FREE-SLIP WALL
                     else if (flag(i-1,j,k,1) == FSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,w)
                        A_m(i,j,k,w) = zero

                        b_m(i-1,j,k) = zero
                        A_m(i-1,j,k,:) = zero
                        A_m(i-1,j,k,0) = -one

! PARTIAL-SLIP WALL
                     else if (flag(i-1,j,k,1) == PSW_) then
                        if (is_undefined(bc_hw_g(l))) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,w)
                           b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,w)*bc_vw_g(l)
                        else
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,w)*&
                              (half*bc_hw_g(l)-odx)/(half*bc_hw_g(l)+odx)
                           b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,w)*&
                              bc_hw_g(l)*bc_vw_g(l)/(half*bc_hw_g(l)+odx)
                        endif
                        A_m(i,j,k,w) = zero

                        A_m(i,j,k-1,:) = zero
                        A_m(i,j,k-1,0) = -one
                        b_m(i,j,k-1) = zero
                     endif

! --- WEST FLUID ---------------------------------------------------------->

! NO-SLIP WALL
                     if (flag(i+1,j,k,1) == NSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,e)
                        A_m(i,j,k,e) = zero

                        b_m(i+1,j,k) = zero
                        A_m(i+1,j,k,:) = zero
                        A_m(i+1,j,k,0) = -one

! FREE-SLIP WALL
                     else if (flag(i+1,j,k,1) == FSW_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,e)
                        A_m(i,j,k,e) = zero

                        b_m(i+1,j,k) = zero
                        A_m(i+1,j,k,:) = zero
                        A_m(i+1,j,k,0) = -one

! PARTIAL-SLIP WALL
                     else if (flag(i+1,j,k,1) == PSW_) then
                        if (is_undefined(bc_hw_g(l))) then
                           A_m(i,j,k,0) = A_m(i,j,k,0)-A_m(i,j,k,e)
                           b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,e)*bc_vw_g(l)
                        else
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,e)*&
                              (half*bc_hw_g(l)-odx)/(half*bc_hw_g(l)+odx)
                           b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,e)*&
                              bc_hw_g(l)*bc_vw_g(l)/(half*bc_hw_g(l)+odx)
                        endif
                        A_m(i,j,k,e) = zero

                        A_m(i,j,k-1,:) = zero
                        A_m(i,j,k-1,0) = -one
                        b_m(i,j,k-1) = zero
                     endif

! --- NORTH FLUID --------------------------------------------------------->

! MASS INFLOW
                     if(is_defined(bc_v_g(l)) .and. (&
                        flag(i,j-1,k,1) == MINF_ .or. &
                        flag(i,j-1,k,1) == MOUT_)) then
                        A_m(i,j-1,k,:) =  zero
                        A_m(i,j-1,k,0) = -one
                        b_m(i,j-1,k) = -bc_v_g(l)
                     endif

! --- SOUTH FLUID --------------------------------------------------------->

! PRESSURE IN/OUTFLOW
                     if(flag(i,j+1,k,1) == PINF_ .or. &
                        flag(i,j+1,k,1) == POUT_) then
                        A_m(i,j,k,0) = A_m(i,j,k,0)+A_m(i,j,k,n)
                        A_m(i,j,k,n) = zero

                        b_m(i,j+1,k) = zero
                        A_m(i,j+1,k,:) = zero
                        A_m(i,j+1,k,0) = -one

! MASS INFLOW
                     else if(is_defined(bc_v_g(l)) .and. (&
                        flag(i,j+1,k,1) == MINF_ .or. &
                        flag(i,j+1,k,1) == MOUT_)) then
                        A_m(i,j,k,:) =  zero
                        A_m(i,j,k,0) = -one
                        b_m(i,j,k) = -bc_v_g(l)

                        b_m(i,j+1,k) = zero
                        A_m(i,j+1,k,:) = zero
                        A_m(i,j+1,k,0) = -one
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
                           b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,b)*bc_vw_g(l)
                        else
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,b)*&
                              (half*bc_hw_g(l)-odz)/(half*bc_hw_g(l)+odz)
                           b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,b)*&
                              bc_hw_g(l)*bc_vw_g(l)/(half*bc_hw_g(l)+odz)
                        endif
                        A_m(i,j,k,b) = zero

                        A_m(i,j,k-1,:) = zero
                        A_m(i,j,k-1,0) = -one
                        b_m(i,j,k-1) = zero
                     endif

! --- BOTTOM FLUID -------------------------------------------------------->

! NO-SLIP WALL
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
                           b_m(i,j,k) = b_m(i,j,k)-2.0*A_m(i,j,k,t)*bc_vw_g(l)
                        else
                           A_m(i,j,k,0) = A_m(i,j,k,0) - A_m(i,j,k,t)*&
                              (half*bc_hw_g(l)-odz)/(half*bc_hw_g(l)+odz)
                           b_m(i,j,k) = b_m(i,j,k) - A_m(i,j,k,t)*&
                              bc_hw_g(l)*bc_vw_g(l)/(half*bc_hw_g(l)+odz)
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

      RETURN
      END SUBROUTINE SOURCE_V_G_BC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_V_G                                        C
!  Purpose: Adds point sources to the gas phase V-Momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_V_G(slo, shi, lo, hi, A_M, B_M, flag, dx, dy, dz)

      use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      use ps, only: ps_v_g, ps_i_e, ps_i_w, ps_j_s, ps_j_n, ps_k_b, ps_k_t

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(IN   ) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: B_m&
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
      INTEGER :: lJN, lJS
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
         if(abs(PS_V_g(PSV)) < small_number) cycle PS_LP

         if(PS_V_g(PSV) < 0.0d0) then
            lJS = PS_J_S(PSV) - 1
            lJN = PS_J_N(PSV) - 1
         else
            lJS = PS_J_S(PSV)
            lJN = PS_J_N(PSV)
         endif

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = lJS, lJN
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.1.eq.flag(i,j,k,1)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

            B_M(I,J,K) = B_M(I,J,K) - pSource * &
               PS_V_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_V_G
end module source_v_g_module
