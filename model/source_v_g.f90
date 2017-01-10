module source_v_g_module

   use bc, only: bc_hw_g, bc_vw_g
   use bc, only: bc_i_w, bc_i_e, bc_j_s, bc_j_n, bc_k_b, bc_k_t
   use bc, only: dimension_bc, bc_defined, bc_type, bc_plane

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
      use compar, only: jmap

      USE constant, only: gravity
      USE bc, only: delp_y

      USE functions, only: avg
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus, jnorth
      USE functions, only: jnorth, jsouth
      USE functions, only: zmax
      USE geometry, only: jmax1, cyclic_y_pd

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
! Phase index
      INTEGER :: m
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

! Set reference phase to gas
      M = 0

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
            IF (CYCLIC_Y_PD) THEN
               IF (JMAP(J).EQ.JMAX1)PGN = P_G(i,jnorth(i,j,k),k) - DELP_Y
            ENDIF
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
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_V_G_BC(slo,shi,lo,hi,A_M,B_M,v_g,flag,dx,dy,dz)

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE functions, only: ieast, iwest, jsouth, jnorth, kbot, ktop
      USE functions, only: jminus, jplus, jm1
      USE geometry  , only: imax2, jmax2, kmax2
      USE geometry  , only: imin2, jmin2, kmin2
      use matrix, only: e, w, s, n, t, b

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
! Phase index
      INTEGER :: M

      real(c_real) :: odx, odz
!-----------------------------------------------

     odx = 1.d0 / dx
     odz = 1.d0 / dz

!-----------------------------------------------

! Set reference phase to gas
      M = 0


! Set the default boundary conditions
! The NS default setting is the where bc_type='dummy' or any default
! (i.e., bc_type=undefined) wall boundary regions are handled. Note that
! the north and south xz planes do not have to be explicitly addressed for
! the v-momentum equation. In this direction the velocities are defined
! at the wall (due staggered grid). They are defined as zero for a
! no penetration condition.
! ---------------------------------------------------------------->>>
! bottom xy plane
      K1 = 1
      if (slo(3) .lt. k1) then
      DO J1 = slo(2),shi(2)
         DO I1 = slo(1),shi(1)
            IF (flag(i1,j1,k1,1) == 100) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and oppostive to the adjacent fluid cell value)
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = -ONE
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (flag(i1,j1,k1,1) == 101) THEN
! Setting the wall velocity equal to the adjacent fluid velocity (set
! the boundary cell value equal to adjacent fluid cell value)
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ONE
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ENDIF
         ENDDO
      ENDDO
      end if

! top xy plane
      K1 = KMAX2
      if (shi(3) .gt. k1) then
      DO J1 = slo(2),shi(2)
         DO I1 = slo(1),shi(1)
            IF (flag(i1,j1,k1,1) == 100) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = -ONE
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (flag(i1,j1,k1,1) == 101) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ONE
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ENDIF
         ENDDO
      ENDDO
      end if

! west zy plane
      I1 = 1
      if (slo(1) .lt. i1) then
      DO K1 = slo(3),shi(3)
         DO J1 = slo(2),shi(2)
            IF (flag(i1,j1,k1,1) == 100) THEN
               A_M(I1,J1,K1,E) = -ONE
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (flag(i1,j1,k1,1) == 101) THEN
               A_M(I1,J1,K1,E) = ONE
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ENDIF
         ENDDO
      ENDDO
      end if

! east zy plane
      I1 = IMAX2
      if (shi(1) .gt. i1) then
      DO K1 = slo(3),shi(3)
         DO J1 = slo(2),shi(2)
            IF (flag(i1,j1,k1,1) == 100) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = -ONE
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (flag(i1,j1,k1,1) == 101) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ONE
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ENDIF
         ENDDO
      ENDDO
      end if
! End setting the default boundary conditions
! ----------------------------------------------------------------<<<

! Setting user specified boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! Setting wall boundary conditions
! ---------------------------------------------------------------->>>
            IF (BC_TYPE(L) == 'NO_SLIP_WALL') THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (flag(i,j,k,1)<100) CYCLE  ! skip redefined cells
                        A_M(I,J,K,E) = ZERO
                        A_M(I,J,K,W) = ZERO
                        A_M(I,J,K,N) = ZERO
                        A_M(I,J,K,S) = ZERO
                        A_M(I,J,K,T) = ZERO
                        A_M(I,J,K,B) = ZERO
                        A_M(I,J,K,0) = -ONE
                        B_M(I,J,K) = ZERO
                        IF (1.eq.flag(ieast(i,j,k),j,k,1)) then
                           A_M(I,J,K,E) = -ONE
                        ELSEIF (1.eq.flag(iwest(i,j,k),j,k,1)) then
                           A_M(I,J,K,W) = -ONE
                        ELSEIF (1.eq.flag(i,j,ktop(i,j,k),1)) then
                           A_M(I,J,K,T) = -ONE
                        ELSEIF (1.eq.flag(i,j,kbot(i,j,k),1)) then
                           A_M(I,J,K,B) = -ONE
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ELSEIF (BC_TYPE(L) == 'FREE_SLIP_WALL') THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (flag(i,j,k,1)<100) CYCLE  ! skip redefined cells
                        A_M(I,J,K,E) = ZERO
                        A_M(I,J,K,W) = ZERO
                        A_M(I,J,K,N) = ZERO
                        A_M(I,J,K,S) = ZERO
                        A_M(I,J,K,T) = ZERO
                        A_M(I,J,K,B) = ZERO
                        A_M(I,J,K,0) = -ONE
                        B_M(I,J,K) = ZERO
                        IF (1.eq.flag(ieast(i,j,k),j,k,1)) then
                           A_M(I,J,K,E) = ONE
                        ELSEIF (1.eq.flag(iwest(i,j,k),j,k,1)) then
                           A_M(I,J,K,W) = ONE
                        ELSEIF (1.eq.flag(i,j,ktop(i,j,k),1)) then
                           A_M(I,J,K,T) = ONE
                        ELSEIF (1.eq.flag(i,j,kbot(i,j,k),1)) then
                           A_M(I,J,K,B) = ONE
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ELSEIF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (flag(i,j,k,1)<100) CYCLE  ! skip redefined cells
                        A_M(I,J,K,E) = ZERO
                        A_M(I,J,K,W) = ZERO
                        A_M(I,J,K,N) = ZERO
                        A_M(I,J,K,S) = ZERO
                        A_M(I,J,K,T) = ZERO
                        A_M(I,J,K,B) = ZERO
                        A_M(I,J,K,0) = -ONE
                        B_M(I,J,K) = ZERO
                        IF (1.eq.flag(ieast(i,j,k),j,k,1)) then
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_M(I,J,K,E) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_VW_G(L)
                           ELSE
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODX)
                              A_M(I,J,K,E) = -(HALF*BC_HW_G(L)-ODX)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (1.eq.flag(iwest(i,j,k),j,k,1)) then
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_M(I,J,K,W) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_VW_G(L)
                           ELSE
                              A_M(I,J,K,W) = -(HALF*BC_HW_G(L)-ODX)
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODX)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (1.eq.flag(i,j,ktop(i,j,k),1)) then
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_M(I,J,K,T) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_VW_G(L)
                           ELSE
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODZ)
                              A_M(I,J,K,T) = -(HALF*BC_HW_G(L)-ODZ)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (1.eq.flag(i,j,kbot(i,j,k),1)) then
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_M(I,J,K,B) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_VW_G(L)
                           ELSE
                              A_M(I,J,K,B) = -(HALF*BC_HW_G(L)-ODZ)
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODZ)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

! end setting of wall boundary conditions
! ----------------------------------------------------------------<<<

! Setting p_inflow or p_outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN
               IF (BC_PLANE(L) == 'S') THEN
! if the fluid cell is on the south side of the outflow/inflow boundary
! then set the velocity in the boundary cell equal to the velocity of
! the adjacent fluid cell
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           A_M(I,J,K,E) = ZERO
                           A_M(I,J,K,W) = ZERO
                           A_M(I,J,K,N) = ZERO
                           A_M(I,J,K,S) = ONE
                           A_M(I,J,K,T) = ZERO
                           A_M(I,J,K,B) = ZERO
                           A_M(I,J,K,0) = -ONE
                           B_M(I,J,K) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
! end setting of p_inflow or p_otuflow flow boundary conditions
! ----------------------------------------------------------------<<<

! Setting outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L) == 'OUTFLOW') THEN
               IF (BC_PLANE(L) == 'S') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           A_M(I,J,K,E) = ZERO
                           A_M(I,J,K,W) = ZERO
                           A_M(I,J,K,N) = ZERO
                           A_M(I,J,K,S) = ONE
                           A_M(I,J,K,T) = ZERO
                           A_M(I,J,K,B) = ZERO
                           A_M(I,J,K,0) = -ONE
                           B_M(I,J,K) = ZERO
                           A_M(i,jminus(i,j,k),k,E) = ZERO
                           A_M(i,jminus(i,j,k),k,W) = ZERO
                           A_M(i,jminus(i,j,k),k,N) = ZERO
                           A_M(i,jminus(i,j,k),k,S) = ONE
                           A_M(i,jminus(i,j,k),k,T) = ZERO
                           A_M(i,jminus(i,j,k),k,B) = ZERO
                           A_M(i,jminus(i,j,k),k,0) = -ONE
                           B_M(i,jminus(i,j,k),k) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ELSEIF (BC_PLANE(L) == 'N') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           A_M(i,jplus(i,j,k),k,E) = ZERO
                           A_M(i,jplus(i,j,k),k,W) = ZERO
                           A_M(i,jplus(i,j,k),k,N) = ONE
                           A_M(i,jplus(i,j,k),k,S) = ZERO
                           A_M(i,jplus(i,j,k),k,T) = ZERO
                           A_M(i,jplus(i,j,k),k,B) = ZERO
                           A_M(i,jplus(i,j,k),k,0) = -ONE
                           B_M(i,jplus(i,j,k),k) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
! end setting of outflow flow boundary conditions
! ----------------------------------------------------------------<<<

! Setting bc that are defined but not nsw, fsw, psw, p_inflow,
! p_outflow, or outflow (at this time, this section addresses
! mass_inflow and mass_outflow type boundaries)
! ---------------------------------------------------------------->>>
            ELSE
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
! setting the velocity in the boundary cell equal to what is known
                        A_M(I,J,K,E) = ZERO
                        A_M(I,J,K,W) = ZERO
                        A_M(I,J,K,N) = ZERO
                        A_M(I,J,K,S) = ZERO
                        A_M(I,J,K,T) = ZERO
                        A_M(I,J,K,B) = ZERO
                        A_M(I,J,K,0) = -ONE
                        B_M(I,J,K) = -V_G(I,J,K)
                        IF (BC_PLANE(L) == 'S') THEN
! if the fluid cell is on the south side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           A_M(i,jsouth(i,j,k),k,E) = ZERO
                           A_M(i,jsouth(i,j,k),k,W) = ZERO
                           A_M(i,jsouth(i,j,k),k,N) = ZERO
                           A_M(i,jsouth(i,j,k),k,S) = ZERO
                           A_M(i,jsouth(i,j,k),k,T) = ZERO
                           A_M(i,jsouth(i,j,k),k,B) = ZERO
                           A_M(i,jsouth(i,j,k),k,0) = -ONE
                           B_M(i,jsouth(i,j,k),k) = -V_G(i,jsouth(i,j,k),k)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF   ! end if/else (bc_type)
                    ! ns, fs, psw; else
                    ! p_inflow, p_outflow, or outflow; else
! end setting of 'else' flow boundary conditions
! (mass_inflow/mass_outflow)
! ----------------------------------------------------------------<<<

         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

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
