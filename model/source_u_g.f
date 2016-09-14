!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_U_g                                              C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector. The center         C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this stage.          C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_U_G(A_M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE bodyforce, only: bfx_g
      USE bc, only: delp_x

      USE compar, only: ijkstart3, ijkend3, imap

      USE fldvar, only: p_g, ro_g, rop_g, rop_go
      USE fldvar, only: ep_g
      USE fldvar, only: u_g, u_go

      USE fun_avg, only: avg_x, avg_z, avg_y
      USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      USE functions, only: ip_at_e, sip_at_e, is_id_at_e
      USE functions, only: ip_of, jp_of, kp_of, im_of, jm_of, km_of
      USE functions, only: east_of, west_of
      USE functions, only: zmax
      USE geometry, only: imax1, cyclic_x_pd
      USE geometry, only: vol, vol_u
      USE geometry, only: ayz


      USE indices, only: i_of, j_of, k_of
      use matrix, only: e, w, s, n, t, b

      USE param, only: dimension_3, dimension_m
      USE param1, only: zero, one, half
      USE run, only: momentum_x_eq
      USE run, only: odt
      USE scales, only: p_scale
      USE fldvar, only: tau_u_g
      USE toleranc, only: dil_ep_s
      USE cutcell, only: cartesian_grid, cut_u_treatment_at
      USE cutcell, only: blocked_u_cell_at
      USE cutcell, only: a_upg_e, a_upg_w
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK, IJKE, IPJK, IJKM, &
                 IPJKM, IMJK, IJMK, IPJMK, IJPK, IJKP
! Phase index
      INTEGER :: M, L, MM
! Internal surface
      INTEGER :: ISV
! Pressure at east cell
      DOUBLE PRECISION :: PgE
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average density
      DOUBLE PRECISION :: ROPGA, ROGA
! Average viscosity
      DOUBLE PRECISION :: MUGA, MUGTA
! Average W_g
      DOUBLE PRECISION :: Wge
! Source terms (Surface)
      DOUBLE PRECISION :: Sdp
! Source terms (Volumetric)
      DOUBLE PRECISION :: V0, Vbf
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION :: Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION :: HYS_drag, avgDrag
! virtual (added) mass
      DOUBLE PRECISION :: ROP_MA, U_se, Usw, Vsw, Vse, Usn,&
                          Uss, Wsb, Wst, Wse, Usb, Ust
! local stress tensor quantity
      DOUBLE PRECISION :: ltau_u_g
!---------------------------------------------------------------------//

! Set reference phase to gas
      M = 0

      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN


!$omp  parallel do default(shared)                                   &
!$omp  private(I, J, K, IJK, IJKE, IJKM, IPJK, IMJK, IPJKM,          &
!$omp          IJMK, IPJMK, IJPK, IJKP, EPGA, PGE, SDP,              &
!$omp           ROPGA, ROGA, ROP_MA, V0, ISV, MUGA, Vbf,   &
!$omp           U_se, Usw, Vsw, Vse, Usn, Uss, Wsb, Wst, Wse,        &
!$omp           Usb, Ust, wGE, MUGTA,              &
!$omp           Ghd_drag, L, MM, avgRop, HYS_drag, avgDrag,          &
!$omp           ltau_u_g)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IJKE = EAST_OF(IJK)
         IJKM = KM_OF(IJK)
         IPJK = IP_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJKM = IP_OF(IJKM)
         IJMK = JM_OF(IJK)
         IPJMK = IP_OF(IJMK)
         IJPK = JP_OF(IJK)
         IJKP = KP_OF(IJK)

         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

! Impermeable internal surface
         IF (IP_AT_E(IJK)) THEN
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO

! Dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO
! set velocity equal to that of west or east cell if solids are present
! in those cells else set velocity equal to known value
            IF (EP_G(WEST_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,W,M) = ONE
            ELSE IF (EP_G(EAST_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,E,M) = ONE
            ELSE
               B_M(IJK,M) = -U_G(IJK)
            ENDIF

! Cartesian grid implementation
         ELSEIF (BLOCKED_U_CELL_AT(IJK)) THEN
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
            PGE = P_G(IJKE)
            IF (CYCLIC_X_PD) THEN
               IF (IMAP(I_OF(IJK)).EQ.IMAX1) PGE = P_G(IJKE) - DELP_X
            ENDIF
            IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                SDP = -P_SCALE*EPGA*(PGE - P_G(IJK))*AYZ(IJK)
            ELSE
                SDP = -P_SCALE*EPGA*(PGE * A_UPG_E(IJK) - P_G(IJK) * A_UPG_W(IJK) )
            ENDIF

            IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
! Volumetric forces
               ROPGA = HALF * (VOL(IJK)*ROP_G(IJK) + &
                               VOL(IPJK)*ROP_G(IJKE))/VOL_U(IJK)
               ROGA  = HALF * (VOL(IJK)*RO_G(IJK) + &
                               VOL(IPJK)*RO_G(IJKE))/VOL_U(IJK)
! Previous time step
               V0 = HALF * (VOL(IJK)*ROP_GO(IJK) + &
                            VOL(IPJK)*ROP_GO(IJKE))*ODT/VOL_U(IJK)
            ELSE
! Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + &
                        VOL(IPJK)*ROP_G(IJKE))/(VOL(IJK) + VOL(IPJK))
               ROGA  = (VOL(IJK)*RO_G(IJK)  + &
                        VOL(IPJK)*RO_G(IJKE) )/(VOL(IJK) + VOL(IPJK))
! Previous time step
               V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IPJK)*ROP_GO(IJKE))*&
                  ODT/(VOL(IJK) + VOL(IPJK))
            ENDIF

! Body force
            VBF = ROGA*BFX_G(IJK)

            ltau_u_g = tau_u_g(ijk)

! Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+&
               A_M(IJK,N,M)+A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+&
               V0*VOL_U(IJK))

            B_M(IJK,M) = B_M(IJK,M) -(SDP + lTAU_U_G + &
               ( (V0)*U_GO(IJK) + VBF)*VOL_U(IJK) )

         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)
      ENDDO   ! end do loop over ijk
!$omp end parallel do

! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_U_G(A_M, B_M)
! modifications for bc
      CALL SOURCE_U_G_BC (A_M, B_M)
! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_U_G_BC(A_M, B_M)

      RETURN
      END SUBROUTINE SOURCE_U_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_U_g_BC                                           C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative. The off-diagonal     C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_U_G_BC(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      use matrix, only: e, w, s, n, t, b
      USE scales
      USE constant
      USE fldvar
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE indices
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
      INTEGER ::  I,  J, K, IM, I1, I2, J1, J2, K1, K2, IJK,&
                  JM, KM, IJKW, IMJK, IP, IPJK
! Phase index
      INTEGER :: M
!-----------------------------------------------

! Set reference phase to gas
      M = 0


! Set the default boundary conditions
! The NS default setting is the where bc_type='dummy' or any default
! (i.e., bc_type=undefined) wall boundary regions are handled. Note that
! the east and west zy planes do not have to be explicitly addressed for
! the u-momentum equation. In this direction the velocities are defined
! at the wall (due staggered grid). They are defined as zero for a
! no penetration condition (see zero_norm_vel subroutine and code under
! the ip_at_e branch in the above source routine).
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
! and opposite to the adjacent fluid cell value)
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
               ELSEIF (FS_WALL_AT(IJK)) THEN
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

! south xz plane
      J1 = 1
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = -ONE
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ONE
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! north xz plane
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = -ONE
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ONE
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
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
                           A_M(IJK,N,M) = -ONE
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           A_M(IJK,S,M) = -ONE
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
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
                           A_M(IJK,N,M) = ONE
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           A_M(IJK,S,M) = ONE
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
                        JM = JM1(J)
                        KM = KM1(K)
                        A_M(IJK,E,M) = ZERO
                        A_M(IJK,W,M) = ZERO
                        A_M(IJK,N,M) = ZERO
                        A_M(IJK,S,M) = ZERO
                        A_M(IJK,T,M) = ZERO
                        A_M(IJK,B,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,N,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(J))
                              A_M(IJK,N,M) = -(HALF*BC_HW_G(L)-ODY_N(J))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,S,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,S,M) = -(HALF*BC_HW_G(L)-ODY_N(JM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(JM))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,T,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,0,M)=-(HALF*BC_HW_G(L)+ODZ_T(K)*OX_E(I))
                              A_M(IJK,T,M)=-(HALF*BC_HW_G(L)-ODZ_T(K)*OX_E(I))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,B,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,B,M) = -(HALF*BC_HW_G(L)-ODZ_T(KM)*OX_E(I&
                                 ))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(KM)*OX_E(I&
                                 ))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO


! Setting p_inflow or p_outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN
               IF (BC_PLANE(L) == 'W') THEN
! if the fluid cell is on the west side of the outflow/inflow boundary
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
                           A_M(IJK,W,M) = ONE
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ZERO
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
               IF (BC_PLANE(L) == 'W') THEN
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
                           A_M(IJK,W,M) = ONE
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ZERO
                           A_M(IJK,T,M) = ZERO
                           A_M(IJK,B,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IM = IM1(I)
                           IMJK = IM_OF(IJK)
                           A_M(IMJK,E,M) = ZERO
                           A_M(IMJK,W,M) = X_E(IM)/X_E(IM1(IM))
                           A_M(IMJK,N,M) = ZERO
                           A_M(IMJK,S,M) = ZERO
                           A_M(IMJK,T,M) = ZERO
                           A_M(IMJK,B,M) = ZERO
                           A_M(IMJK,0,M) = -ONE
                           B_M(IMJK,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ELSEIF (BC_PLANE(L) == 'E') THEN
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
                           IP = IP1(I)
                           IPJK = IP_OF(IJK)
                           A_M(IPJK,E,M) = X_E(IP)/X_E(I)
                           A_M(IPJK,W,M) = ZERO
                           A_M(IPJK,N,M) = ZERO
                           A_M(IPJK,S,M) = ZERO
                           A_M(IPJK,T,M) = ZERO
                           A_M(IPJK,B,M) = ZERO
                           A_M(IPJK,0,M) = -ONE
                           B_M(IPJK,M) = ZERO
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
                        B_M(IJK,M) = -U_G(IJK)
                        IF (BC_PLANE(L) == 'W') THEN
! if the fluid cell is on the west side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           IJKW = WEST_OF(IJK)
                           A_M(IJKW,E,M) = ZERO
                           A_M(IJKW,W,M) = ZERO
                           A_M(IJKW,N,M) = ZERO
                           A_M(IJKW,S,M) = ZERO
                           A_M(IJKW,T,M) = ZERO
                           A_M(IJKW,B,M) = ZERO
                           A_M(IJKW,0,M) = -ONE
                           B_M(IJKW,M) = -U_G(IJKW)
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
      SUBROUTINE POINT_SOURCE_U_G(A_M, B_M)

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
      INTEGER :: lIE, lIW
! terms of bm expression
      DOUBLE PRECISION :: pSource
!-----------------------------------------------

! Set reference phase to gas
      M = 0

! Calculate the mass going into each IJK cell. This is done for each
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

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(.NOT. fluid_at(ijk)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

            B_M(IJK,M) = B_M(IJK,M) - pSource *                        &
               PS_U_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo
      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_U_G
