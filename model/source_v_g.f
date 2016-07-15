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
      SUBROUTINE SOURCE_V_G(A_M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE bodyforce, only: bfy_g
      USE bc, only: delp_y

      USE compar, only: ijkstart3, ijkend3, jmap

      USE fldvar, only: p_g, ro_g, rop_g, rop_go
      USE fldvar, only: ep_g
      USE fldvar, only: v_g, v_go

      USE fun_avg, only: avg_x, avg_z, avg_y
      USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      USE functions, only: ip_at_n, sip_at_n, is_id_at_n
      USE functions, only: ip_of, jp_of, kp_of, im_of, jm_of, km_of
      USE functions, only: north_of, south_of
      USE functions, only: zmax
      USE geometry, only: jmax1, cyclic_y_pd, do_k, xlength
      USE geometry, only: vol, vol_v
      USE geometry, only: axy, ayz, axz
      USE geometry, only: ox, odx_e


      USE indices, only: i_of, j_of, k_of
      USE indices, only: ip1, im1, kp1
      USE is, only: is_pc
      USE ambm, only: e, w, s, n, t, b

      USE param, only: dimension_3, dimension_m
      USE param1, only: zero, one, half
      USE fldvar, only: mu_g
      USE run, only: momentum_y_eq
      USE run, only: odt
      USE scales, only: p_scale
      USE fldvar, only: tau_v_g
      USE toleranc, only: dil_ep_s
      USE cutcell, only: cartesian_grid, cut_v_treatment_at
      USE cutcell, only: blocked_v_cell_at
      USE cutcell, only: a_vpg_n, a_vpg_s
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK, IJKN, &
                 IMJK, IPJK, IJMK, IJPK, IJKP, IJKM, IMJPK, IJPKM
! Phase index
      INTEGER :: M, L, MM
! Internal surface
      INTEGER :: ISV
! Pressure at north cell
      DOUBLE PRECISION :: PgN
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average density
      DOUBLE PRECISION :: ROPGA, ROGA
! Average viscosity
      DOUBLE PRECISION :: MUGA
! Source terms (Surface)
      DOUBLE PRECISION :: Sdp
! Source terms (Volumetric)
      DOUBLE PRECISION :: V0, Vpm, Vbf
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION :: Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION :: HYS_drag, avgDrag
! virtual (added) mass
      DOUBLE PRECISION :: ROP_MA, Vsn, Vss, U_se, Usw, Vse, Vsw, &
                          Wst, Wsb, Vst, Vsb
! jackson terms: local stress tensor quantity
      DOUBLE PRECISION :: ltau_v_g
!---------------------------------------------------------------------//

! Set reference phase to gas
      M = 0

      IF (.NOT.MOMENTUM_Y_EQ(0)) RETURN


!$omp  parallel do default(shared)                                   &
!$omp  private(I, J, K, IJK, IJKN, IMJK, IPJK, IJMK, IJPK, IMJPK,    &
!$omp          IJKM, IJPKM, IJKP, EPGA, PGN, SDP, ROPGA,             &
!$omp          ROGA, ROP_MA, V0, ISV, MUGA, Vpm, Vbf, L, MM,    &
!$omp          Vsn, Vss, U_se, Usw, Vse, Vsw, Wst, Wsb, Vst,         &
!$omp          Vsb,ghd_drag, avgRop, avgDrag, HYS_drag,      &
!$omp          ltau_v_g)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IJKN = NORTH_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJK = IP_OF(IJK)
         IJMK = JM_OF(IJK)
         IJPK = JP_OF(IJK)
         IMJPK = IM_OF(IJPK)
         IJKM = KM_OF(IJK)
         IJPKM = KM_OF(IJPK)
         IJKP = KP_OF(IJK)

         EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)

! Impermeable internal surface
         IF (IP_AT_N(IJK)) THEN
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO

! dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO
            IF (EP_G(SOUTH_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,S,M) = ONE
            ELSE IF (EP_G(NORTH_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,N,M) = ONE
            ELSE
               B_M(IJK,M) = -V_G(IJK)
            ENDIF

! Cartesian grid implementation
         ELSEIF (BLOCKED_V_CELL_AT(IJK)) THEN
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO

! Normal case
         ELSE

! Surface forces
! Pressure term
            PGN = P_G(IJKN)
            IF (CYCLIC_Y_PD) THEN
               IF (JMAP(J_OF(IJK)).EQ.JMAX1)PGN = P_G(IJKN) - DELP_Y
            ENDIF
            IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
               SDP = -P_SCALE*EPGA*(PGN - P_G(IJK))*AXZ(IJK)
            ELSE
               SDP = -P_SCALE*EPGA*(PGN * A_VPG_N(IJK) - &
                                    P_G(IJK) * A_VPG_S(IJK) )
            ENDIF

            IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
! Volumetric forces
               ROPGA = AVG_Y(ROP_G(IJK),ROP_G(IJKN),J)
               ROGA = AVG_Y(RO_G(IJK),RO_G(IJKN),J)
! Previous time step
               V0 = AVG_Y(ROP_GO(IJK),ROP_GO(IJKN),J)*ODT
            ELSE
! Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + &
                  VOL(IJKN)*ROP_G(IJKN))/(VOL(IJK) + VOL(IJKN))
               ROGA  = (VOL(IJK)*RO_G(IJK)  + &
                  VOL(IJKN)*RO_G(IJKN) )/(VOL(IJK) + VOL(IJKN))
! Previous time step
               V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IJKN)*ROP_GO(IJKN))*&
                  ODT/(VOL(IJK) + VOL(IJKN))
            ENDIF

! pressure drop through porous media
            IF (SIP_AT_N(IJK)) THEN
               ISV = IS_ID_AT_N(IJK)
               MUGA = AVG_Y(MU_G(IJK),MU_G(IJKN),J)
               VPM = MUGA/IS_PC(ISV,1)
               IF (IS_PC(ISV,2) /= ZERO) VPM = VPM + HALF*IS_PC(ISV,2)*&
                   ROPGA*ABS(V_G(IJK))
            ELSE
               VPM = ZERO
            ENDIF

! Body force
            VBF = ROGA*BFY_G(IJK)

! if jackson, implement jackson form of governing equations (ep_g dot
! del tau_g): multiply by void fraction otherwise by 1
            ltau_v_g = tau_v_g(ijk)


! Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+&
               A_M(IJK,N,M)+A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+&
               (V0+VPM)*VOL_V(IJK))
            B_M(IJK,M) = B_M(IJK,M) - (SDP + lTAU_V_G +  &
               ((V0)*V_GO(IJK) + VBF)*VOL_V(IJK) )

         ENDIF
      ENDDO
!$omp end parallel do

! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_V_G(A_M, B_M)
! modifications for bc
      CALL SOURCE_V_G_BC(A_M, B_M)
! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_V_G_BC(A_M, B_M)

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

      SUBROUTINE SOURCE_V_G_BC(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      use ambm, only: e, w, s, n, t, b
      USE scales
      USE constant
      USE fldvar
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK, &
                 IM, KM, IJKS, IJMK, IJPK
! Phase index
      INTEGER :: M
! Turbulent shear stress
      DOUBLE PRECISION :: W_F_Slip
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
      IF (DO_K) THEN
! bottom xy plane
         K1 = 1
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (NS_WALL_AT(IJK)) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and oppostive to the adjacent fluid cell value)
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = -ONE
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ELSEIF (FS_WALL_AT(IJK)) THEN
! Setting the wall velocity equal to the adjacent fluid velocity (set
! the boundary cell value equal to adjacent fluid cell value)
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ONE
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO

! top xy plane
         K1 = KMAX2
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (NS_WALL_AT(IJK)) THEN
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = -ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ELSE IF (FS_WALL_AT(IJK)) THEN
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO
      ENDIF   ! end if (do_k)


! west zy plane
      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = -ONE
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ONE
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! east zy plane
      I1 = IMAX2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = -ONE
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ONE
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
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
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
                        A_M(IJK,E,M) = ZERO
                        A_M(IJK,W,M) = ZERO
                        A_M(IJK,N,M) = ZERO
                        A_M(IJK,S,M) = ZERO
                        A_M(IJK,T,M) = ZERO
                        A_M(IJK,B,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           A_M(IJK,E,M) = -ONE
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           A_M(IJK,W,M) = -ONE
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           A_M(IJK,T,M) = -ONE
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           A_M(IJK,B,M) = -ONE
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
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
                        A_M(IJK,E,M) = ZERO
                        A_M(IJK,W,M) = ZERO
                        A_M(IJK,N,M) = ZERO
                        A_M(IJK,S,M) = ZERO
                        A_M(IJK,T,M) = ZERO
                        A_M(IJK,B,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           A_M(IJK,E,M) = ONE
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           A_M(IJK,W,M) = ONE
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           A_M(IJK,T,M) = ONE
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           A_M(IJK,B,M) = ONE
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
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
                        IM = IM1(I)
                        KM = KM1(K)
                        A_M(IJK,E,M) = ZERO
                        A_M(IJK,W,M) = ZERO
                        A_M(IJK,N,M) = ZERO
                        A_M(IJK,S,M) = ZERO
                        A_M(IJK,T,M) = ZERO
                        A_M(IJK,B,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,E,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_VW_G(L)
                           ELSE
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(I))
                              A_M(IJK,E,M) = -(HALF*BC_HW_G(L)-ODX_E(I))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,W,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_VW_G(L)
                           ELSE
                              A_M(IJK,W,M) = -(HALF*BC_HW_G(L)-ODX_E(IM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(IM))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,T,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_VW_G(L)
                           ELSE
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(K)*OX(I))
                              A_M(IJK,T,M) = -(HALF*BC_HW_G(L)-ODZ_T(K)*OX(I))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,B,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_VW_G(L)
                           ELSE
                              A_M(IJK,B,M) = -(HALF*BC_HW_G(L)-ODZ_T(KM)*OX(I))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(KM)*OX(I))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_VW_G(L)
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
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           A_M(IJK,E,M) = ZERO
                           A_M(IJK,W,M) = ZERO
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ONE
                           A_M(IJK,T,M) = ZERO
                           A_M(IJK,B,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
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
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           A_M(IJK,E,M) = ZERO
                           A_M(IJK,W,M) = ZERO
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ONE
                           A_M(IJK,T,M) = ZERO
                           A_M(IJK,B,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IJMK = JM_OF(IJK)
                           A_M(IJMK,E,M) = ZERO
                           A_M(IJMK,W,M) = ZERO
                           A_M(IJMK,N,M) = ZERO
                           A_M(IJMK,S,M) = ONE
                           A_M(IJMK,T,M) = ZERO
                           A_M(IJMK,B,M) = ZERO
                           A_M(IJMK,0,M) = -ONE
                           B_M(IJMK,M) = ZERO
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
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IJPK = JP_OF(IJK)
                           A_M(IJPK,E,M) = ZERO
                           A_M(IJPK,W,M) = ZERO
                           A_M(IJPK,N,M) = ONE
                           A_M(IJPK,S,M) = ZERO
                           A_M(IJPK,T,M) = ZERO
                           A_M(IJPK,B,M) = ZERO
                           A_M(IJPK,0,M) = -ONE
                           B_M(IJPK,M) = ZERO
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
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
! setting the velocity in the boundary cell equal to what is known
                        A_M(IJK,E,M) = ZERO
                        A_M(IJK,W,M) = ZERO
                        A_M(IJK,N,M) = ZERO
                        A_M(IJK,S,M) = ZERO
                        A_M(IJK,T,M) = ZERO
                        A_M(IJK,B,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = -V_G(IJK)
                        IF (BC_PLANE(L) == 'S') THEN
! if the fluid cell is on the south side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           IJKS = SOUTH_OF(IJK)
                           A_M(IJKS,E,M) = ZERO
                           A_M(IJKS,W,M) = ZERO
                           A_M(IJKS,N,M) = ZERO
                           A_M(IJKS,S,M) = ZERO
                           A_M(IJKS,T,M) = ZERO
                           A_M(IJKS,B,M) = ZERO
                           A_M(IJKS,0,M) = -ONE
                           B_M(IJKS,M) = -V_G(IJKS)
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

      use compar
      use constant
      use geometry
      use indices
      use param1, only: small_number
      use fldvar
      use ps
      use run
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV, M
      INTEGER :: lJN, lJS
! terms of bm expression
      DOUBLE PRECISION :: pSource
!-----------------------------------------------

! Set reference phase to gas
      M = 0

! Calculate the mass going into each IJK cell. This is done for each
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

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(.NOT.fluid_at(ijk)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

            B_M(IJK,M) = B_M(IJK,M) - pSource * &
               PS_V_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_V_G
