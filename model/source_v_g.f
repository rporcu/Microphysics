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
      USE constant, only: gravity_y
      USE bc, only: delp_y

      USE fldvar, only: p_g, ro_g, rop_g, rop_go
      USE fldvar, only: ep_g
      USE fldvar, only: v_g, v_go

      USE fun_avg, only: avg
      USE functions, only: ip_at_n, sip_at_n, is_id_at_n
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus, jnorth
      USE functions, only: jnorth, jsouth
      USE functions, only: zmax, funijk, wall_cell
      USE geometry, only: jmax1, cyclic_y_pd
      USE geometry, only: vol
      USE geometry, only: axz

      use matrix, only: e, w, s, n, t, b

      USE param, only: dimension_3
      USE param1, only: zero, one, half
      USE run, only: momentum_y_eq
      USE run, only: odt
      USE scales, only: p_scale
      USE fldvar, only: tau_v_g
      USE toleranc, only: dil_ep_s
      use compar, only: jmap
      use compar, only: istart2, iend2
      use compar, only: jstart2, jend2
      use compar, only: kstart2, kend2
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK, IJKN, &
                 IMJK, IPJK, IJMK, IJPK, IJKP, IJKM, IMJPK, IJPKM
! Phase index
      INTEGER :: M
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

          IJK = FUNIJK(i,j,k)

          IJKN = FUNIJK(i,jnorth(i,j,k),k)
          IMJK = FUNIJK(iminus(i,j,k),j,k)
          IPJK = FUNIJK(iplus(i,j,k),j,k)
          IJMK = FUNIJK(i,jminus(i,j,k),k)
          IJPK = FUNIJK(i,jplus(i,j,k),k)
          IJKM = FUNIJK(i,j,kminus(i,j,k))
          IJKP = FUNIJK(i,j,kplus(i,j,k))

          IMJPK = FUNIJK(IMINUS(I,jplus(i,j,k),K),jplus(i,j,k),k)
          IJPKM = FUNIJK(I,jplus(i,j,k),KMINUS(I,jplus(i,j,k),k))

         EPGA = AVG(EP_G(IJK),EP_G(IJKN))

! Impermeable internal surface
         IF (ip_at_n(i,j,k)) THEN
            A_M(IJK,E) = ZERO
            A_M(IJK,W) = ZERO
            A_M(IJK,N) = ZERO
            A_M(IJK,S) = ZERO
            A_M(IJK,T) = ZERO
            A_M(IJK,B) = ZERO
            A_M(IJK,0) = -ONE
            B_M(IJK) = ZERO

! dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_M(IJK,E) = ZERO
            A_M(IJK,W) = ZERO
            A_M(IJK,N) = ZERO
            A_M(IJK,S) = ZERO
            A_M(IJK,T) = ZERO
            A_M(IJK,B) = ZERO
            A_M(IJK,0) = -ONE
            B_M(IJK) = ZERO
            IF (EP_G(FUNIJK(i,jsouth(i,j,k),k)) > DIL_EP_S) THEN
               A_M(IJK,S) = ONE
            ELSE IF (EP_G(FUNIJK(i,jnorth(i,j,k),k)) > DIL_EP_S) THEN
               A_M(IJK,N) = ONE
            ELSE
               B_M(IJK) = -V_G(IJK)
            ENDIF

! Normal case
         ELSE

! Surface forces
! Pressure term
            PGN = P_G(IJKN)
            IF (CYCLIC_Y_PD) THEN
               IF (JMAP(J).EQ.JMAX1)PGN = P_G(IJKN) - DELP_Y
            ENDIF
            SDP = -P_SCALE*EPGA*(PGN - P_G(IJK))*AXZ

! Volumetric forces
            ROPGA = AVG(ROP_G(IJK),ROP_G(IJKN))
            ROGA = AVG(RO_G(IJK),RO_G(IJKN))
! Previous time step
            V0 = AVG(ROP_GO(I,J,K),ROP_GO(i,jnorth(i,j,k),k))*ODT

! Body force
            VBF = ROGA*GRAVITY_Y

! if jackson, implement jackson form of governing equations (ep_g dot
! del tau_g): multiply by void fraction otherwise by 1
            ltau_v_g = tau_v_g(ijk)


! Collect the terms
            A_M(IJK,0) = -(A_M(IJK,E)+A_M(IJK,W)+&
               A_M(IJK,N)+A_M(IJK,S)+A_M(IJK,T)+A_M(IJK,B)+&
               V0*VOL)
            B_M(IJK) = B_M(IJK) - (SDP + lTAU_V_G +  &
               ((V0)*V_GO(I,J,K) + VBF)*VOL )

         ENDIF
      ENDDO
      ENDDO
      ENDDO

! modifications for bc
      CALL SOURCE_V_G_BC(A_M, B_M)

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
      use matrix, only: e, w, s, n, t, b
      USE scales
      USE constant
      USE fldvar
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions, only: is_on_mype_plus2layers
      USE functions, only: funijk, wall_cell, fluid_cell
      USE functions, only: ieast,iwest,jnorth,jsouth,kbot,ktop
      USE functions, only: jminus,jplus
      USE functions, only: fs_wall_cell, ns_wall_cell
      USE functions, only: im1, jm1, km1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
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
               IF (ns_wall_cell(i1,j1,k1)) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and oppostive to the adjacent fluid cell value)
                  A_M(IJK,E) = ZERO
                  A_M(IJK,W) = ZERO
                  A_M(IJK,N) = ZERO
                  A_M(IJK,S) = ZERO
                  A_M(IJK,T) = -ONE
                  A_M(IJK,B) = ZERO
                  A_M(IJK,0) = -ONE
                  B_M(IJK) = ZERO
               ELSEIF (fs_wall_cell(i1,j1,k1)) THEN
! Setting the wall velocity equal to the adjacent fluid velocity (set
! the boundary cell value equal to adjacent fluid cell value)
                  A_M(IJK,E) = ZERO
                  A_M(IJK,W) = ZERO
                  A_M(IJK,N) = ZERO
                  A_M(IJK,S) = ZERO
                  A_M(IJK,T) = ONE
                  A_M(IJK,B) = ZERO
                  A_M(IJK,0) = -ONE
                  B_M(IJK) = ZERO
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
               IF (ns_wall_cell(i1,j1,k1)) THEN
                  A_M(IJK,E) = ZERO
                  A_M(IJK,W) = ZERO
                  A_M(IJK,N) = ZERO
                  A_M(IJK,S) = ZERO
                  A_M(IJK,T) = ZERO
                  A_M(IJK,B) = -ONE
                  A_M(IJK,0) = -ONE
                  B_M(IJK) = ZERO
               ELSEIF (fs_wall_cell(i1,j1,k1)) THEN
                  A_M(IJK,E) = ZERO
                  A_M(IJK,W) = ZERO
                  A_M(IJK,N) = ZERO
                  A_M(IJK,S) = ZERO
                  A_M(IJK,T) = ZERO
                  A_M(IJK,B) = ONE
                  A_M(IJK,0) = -ONE
                  B_M(IJK) = ZERO
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
            IF (ns_wall_cell(i1,j1,k1)) THEN
               A_M(IJK,E) = -ONE
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
            ELSEIF (fs_wall_cell(i1,j1,k1)) THEN
               A_M(IJK,E) = ONE
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
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
            IF (ns_wall_cell(i1,j1,k1)) THEN
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = -ONE
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
            ELSEIF (fs_wall_cell(i1,j1,k1)) THEN
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = ONE
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
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
                        IF (.NOT.wall_cell(i,j,k)) CYCLE  ! skip redefined cells
                        A_M(IJK,E) = ZERO
                        A_M(IJK,W) = ZERO
                        A_M(IJK,N) = ZERO
                        A_M(IJK,S) = ZERO
                        A_M(IJK,T) = ZERO
                        A_M(IJK,B) = ZERO
                        A_M(IJK,0) = -ONE
                        B_M(IJK) = ZERO
                        IF (fluid_cell(ieast(i,j,k),j,k)) then
                           A_M(IJK,E) = -ONE
                        ELSEIF (fluid_cell(iwest(i,j,k),j,k)) then
                           A_M(IJK,W) = -ONE
                        ELSEIF (fluid_cell(i,j,ktop(i,j,k))) then
                           A_M(IJK,T) = -ONE
                        ELSEIF (fluid_cell(i,j,kbot(i,j,k))) then
                           A_M(IJK,B) = -ONE
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
                        IF (.NOT.wall_cell(i,j,k)) CYCLE  ! skip redefined cells
                        A_M(IJK,E) = ZERO
                        A_M(IJK,W) = ZERO
                        A_M(IJK,N) = ZERO
                        A_M(IJK,S) = ZERO
                        A_M(IJK,T) = ZERO
                        A_M(IJK,B) = ZERO
                        A_M(IJK,0) = -ONE
                        B_M(IJK) = ZERO
                        IF (fluid_cell(ieast(i,j,k),j,k)) then
                           A_M(IJK,E) = ONE
                        ELSEIF (fluid_cell(iwest(i,j,k),j,k)) then
                           A_M(IJK,W) = ONE
                        ELSEIF (fluid_cell(i,j,ktop(i,j,k))) then
                           A_M(IJK,T) = ONE
                        ELSEIF (fluid_cell(i,j,kbot(i,j,k))) then
                           A_M(IJK,B) = ONE
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
                        IF (.NOT.wall_cell(i,j,k)) CYCLE  ! skip redefined cells
                        IM = IM1(I)
                        KM = KM1(K)
                        A_M(IJK,E) = ZERO
                        A_M(IJK,W) = ZERO
                        A_M(IJK,N) = ZERO
                        A_M(IJK,S) = ZERO
                        A_M(IJK,T) = ZERO
                        A_M(IJK,B) = ZERO
                        A_M(IJK,0) = -ONE
                        B_M(IJK) = ZERO
                        IF (fluid_cell(ieast(i,j,k),j,k)) then
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,E) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_VW_G(L)
                           ELSE
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODX_E(I))
                              A_M(IJK,E) = -(HALF*BC_HW_G(L)-ODX_E(I))
                              B_M(IJK) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (fluid_cell(iwest(i,j,k),j,k)) then
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,W) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_VW_G(L)
                           ELSE
                              A_M(IJK,W) = -(HALF*BC_HW_G(L)-ODX_E(IM))
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODX_E(IM))
                              B_M(IJK) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (fluid_cell(i,j,ktop(i,j,k))) then
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,T) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_VW_G(L)
                           ELSE
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODZ_T(K))
                              A_M(IJK,T) = -(HALF*BC_HW_G(L)-ODZ_T(K))
                              B_M(IJK) = -BC_HW_G(L)*BC_VW_G(L)
                           ENDIF
                        ELSEIF (fluid_cell(i,j,kbot(i,j,k))) then
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,B) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_VW_G(L)
                           ELSE
                              A_M(IJK,B) = -(HALF*BC_HW_G(L)-ODZ_T(KM))
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODZ_T(KM))
                              B_M(IJK) = -BC_HW_G(L)*BC_VW_G(L)
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
                           A_M(IJK,E) = ZERO
                           A_M(IJK,W) = ZERO
                           A_M(IJK,N) = ZERO
                           A_M(IJK,S) = ONE
                           A_M(IJK,T) = ZERO
                           A_M(IJK,B) = ZERO
                           A_M(IJK,0) = -ONE
                           B_M(IJK) = ZERO
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
                           A_M(IJK,E) = ZERO
                           A_M(IJK,W) = ZERO
                           A_M(IJK,N) = ZERO
                           A_M(IJK,S) = ONE
                           A_M(IJK,T) = ZERO
                           A_M(IJK,B) = ZERO
                           A_M(IJK,0) = -ONE
                           B_M(IJK) = ZERO
                           IJMK = FUNIJK(i,jminus(i,j,k),k)
                           A_M(IJMK,E) = ZERO
                           A_M(IJMK,W) = ZERO
                           A_M(IJMK,N) = ZERO
                           A_M(IJMK,S) = ONE
                           A_M(IJMK,T) = ZERO
                           A_M(IJMK,B) = ZERO
                           A_M(IJMK,0) = -ONE
                           B_M(IJMK) = ZERO
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
                           IJMK = FUNIJK(i,jplus(i,j,k),k)
                           A_M(IJPK,E) = ZERO
                           A_M(IJPK,W) = ZERO
                           A_M(IJPK,N) = ONE
                           A_M(IJPK,S) = ZERO
                           A_M(IJPK,T) = ZERO
                           A_M(IJPK,B) = ZERO
                           A_M(IJPK,0) = -ONE
                           B_M(IJPK) = ZERO
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
                        A_M(IJK,E) = ZERO
                        A_M(IJK,W) = ZERO
                        A_M(IJK,N) = ZERO
                        A_M(IJK,S) = ZERO
                        A_M(IJK,T) = ZERO
                        A_M(IJK,B) = ZERO
                        A_M(IJK,0) = -ONE
                        B_M(IJK) = -V_G(IJK)
                        IF (BC_PLANE(L) == 'S') THEN
! if the fluid cell is on the south side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           IJKS = FUNIJK(i,jsouth(i,j,k),k)
                           A_M(IJKS,E) = ZERO
                           A_M(IJKS,W) = ZERO
                           A_M(IJKS,N) = ZERO
                           A_M(IJKS,S) = ZERO
                           A_M(IJKS,T) = ZERO
                           A_M(IJKS,B) = ZERO
                           A_M(IJKS,0) = -ONE
                           B_M(IJKS) = -V_G(IJKS)
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
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
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
            if(.NOT.fluid_cell(i,j,k)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

            B_M(IJK) = B_M(IJK) - pSource * &
               PS_V_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_V_G
