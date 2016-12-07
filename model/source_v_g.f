module source_v_g_module

   use bc, only: bc_hw_g, bc_vw_g
   use bc, only: bc_i_w, bc_i_e, bc_j_s, bc_j_n, bc_k_b, bc_k_t
   use bc, only: dimension_bc, bc_defined, bc_type, bc_plane
   use param1, only: zero, half, one, undefined

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
      SUBROUTINE SOURCE_V_G(A_M, B_M, p_g, ep_g, ro_g, rop_g, rop_go, &
                            v_g, v_go, tau_v_g,flag)


! Modules
!---------------------------------------------------------------------//
      use compar, only: jmap
      USE compar  , only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3

      USE constant, only: gravity
      USE bc, only: delp_y

      USE functions, only: avg
      USE functions, only: ip_at_n
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus, jnorth
      USE functions, only: jnorth, jsouth
      USE functions, only: zmax
      USE geometry, only: jmax1, cyclic_y_pd
      USE geometry, only: vol, axz

      use matrix, only: e, w, s, n, t, b

      USE param1, only: zero, one
      USE run, only: momentum_y_eq, odt
      USE scales, only: p_scale
      USE toleranc, only: dil_ep_s

      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, INTENT(IN   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: tau_v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: i,j,k
! Phase index
      INTEGER :: m
! Pressure at north cell
      DOUBLE PRECISION :: PgN
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average density
      DOUBLE PRECISION :: ROPGA, ROGA
! Source terms (Surface)
      DOUBLE PRECISION :: Sdp
! Source terms (Volumetric)
      DOUBLE PRECISION :: V0, Vbf
! jackson terms: local stress tensor quantity
      DOUBLE PRECISION :: ltau_v_g
!---------------------------------------------------------------------//

! Set reference phase to gas
      M = 0

      IF (.NOT.MOMENTUM_Y_EQ(0)) RETURN

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

           EPGA = AVG(EP_G(I,J,K),EP_G(i,jnorth(i,j,k),k))

! Impermeable internal surface
         IF (ip_at_n(i,j,k)) THEN
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
      CALL SOURCE_V_G_BC(A_M, B_M, V_G,flag)

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

      SUBROUTINE SOURCE_V_G_BC(A_M, B_M, v_g,flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE functions, only: ieast, iwest, jsouth, jnorth, kbot, ktop
      USE functions, only: jminus, jplus, jm1
      USE functions, only: fluid_at
      USE geometry  , only: imax2, kmax2
      USE geometry  , only: imin3, imax3, jmin3, jmax3, kmin3, kmax3
      USE geometry  , only: odx, odz
      use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use matrix, only: e, w, s, n, t, b

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Velocity v_g
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2
! Phase index
      INTEGER :: M
!-----------------------------------------------

! Set reference phase to gas
      M = 0


! Set the default boundary conditions
! The NS default setting is the where bc_type='dummy' or any default
! (i.e., bc_type=undefined) wall boundary regions are handled. Note that
! the north and south xz planes do not have to be explicitly addressed for
! the v-momentum equation. In this direction the velocities are defined
! at the wall (due staggered grid). They are defined as zero for a
! no penetration condition (see zero_norm_vel subroutine and code under
! the ip_at_n branch in the above source routine).
! ---------------------------------------------------------------->>>
! bottom xy plane
      K1 = 1
      DO J1 = jmin3,jmax3
         DO I1 = imin3, imax3
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

! top xy plane
      K1 = KMAX2
      DO J1 = jmin3,jmax3
         DO I1 = imin3, imax3
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

! west zy plane
      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
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

! east zy plane
      I1 = IMAX2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
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
                        IF (fluid_at(ieast(i,j,k),j,k)) then
                           A_M(I,J,K,E) = -ONE
                        ELSEIF (fluid_at(iwest(i,j,k),j,k)) then
                           A_M(I,J,K,W) = -ONE
                        ELSEIF (fluid_at(i,j,ktop(i,j,k))) then
                           A_M(I,J,K,T) = -ONE
                        ELSEIF (fluid_at(i,j,kbot(i,j,k))) then
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
                        IF (fluid_at(ieast(i,j,k),j,k)) then
                           A_M(I,J,K,E) = ONE
                        ELSEIF (fluid_at(iwest(i,j,k),j,k)) then
                           A_M(I,J,K,W) = ONE
                        ELSEIF (fluid_at(i,j,ktop(i,j,k))) then
                           A_M(I,J,K,T) = ONE
                        ELSEIF (fluid_at(i,j,kbot(i,j,k))) then
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
                        IF (fluid_at(ieast(i,j,k),j,k)) then
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(I,J,K,E) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_VW_G(L)
                           ELSE
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODX)
                              A_M(I,J,K,E) = -(HALF*BC_HW_G(L)-ODX)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (fluid_at(iwest(i,j,k),j,k)) then
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(I,J,K,W) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_VW_G(L)
                           ELSE
                              A_M(I,J,K,W) = -(HALF*BC_HW_G(L)-ODX)
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODX)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (fluid_at(i,j,ktop(i,j,k))) then
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(I,J,K,T) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_VW_G(L)
                           ELSE
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODZ)
                              A_M(I,J,K,T) = -(HALF*BC_HW_G(L)-ODZ)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (fluid_at(i,j,kbot(i,j,k))) then
                           IF (BC_HW_G(L) == UNDEFINED) THEN
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
      SUBROUTINE POINT_SOURCE_V_G(A_M, B_M)

      use functions, only: fluid_at
      use geometry, only: vol
      use param1  , only: small_number
      use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      use ps, only: ps_v_g, ps_i_e, ps_i_w, ps_j_s, ps_j_n, ps_k_b, ps_k_t

      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN   ) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K
      INTEGER :: PSV, M
      INTEGER :: lJN, lJS
! terms of bm expression
      DOUBLE PRECISION :: pSource
!-----------------------------------------------

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

            if(.NOT.fluid_at(i,j,k)) cycle

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
