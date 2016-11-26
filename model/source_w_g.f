!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_W_g                                              C
!  Purpose: Determine source terms for W_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
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
      SUBROUTINE SOURCE_W_G(A_M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: gravity_z
      USE bc, only: delp_z

      USE fldvar, only: p_g, ro_g, rop_g, rop_go
      USE fldvar, only: ep_g
      USE fldvar, only: w_g, w_go

      USE fun_avg, only: avg_h
      USE fun_avg, only: avg
      USE functions, only: ip_at_t, sip_at_t, is_id_at_t
      USE functions, only: ieast, iwest, jnorth, jsouth, kbot, ktop
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus,ktop
      USE functions, only: zmax, funijk, wall_at
      USE geometry, only: kmax1, cyclic_z_pd
      USE geometry, only: vol
      USE geometry, only: axy

      use matrix, only: e, w, s, n, t, b

      USE param, only: dimension_3
      USE param1, only: zero, one, half
      USE run, only: momentum_z_eq
      USE run, only: odt
      USE scales, only: p_scale
      USE fldvar, only: tau_w_g
      USE toleranc, only: dil_ep_s
      use compar, only: kmap
      use compar, only: istart2, iend2, istart3, iend3
      use compar, only: jstart2, jend2, jstart3, jend3
      use compar, only: kstart2, kend2, kstart3, kend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK, IJKT, IMJK, IJKP, IMJKP,&
                 IPJK, IJKM, IJMK, IJMKP, IJPK
! Phase index
      INTEGER :: M
! Pressure at top cell
      DOUBLE PRECISION :: PgT
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average density
      DOUBLE PRECISION :: ROPGA, ROGA
! Source terms (Surface)
      DOUBLE PRECISION Sdp
! Source terms (Volumetric)
      DOUBLE PRECISION V0, Vbf
! jackson terms: local stress tensor quantity
      DOUBLE PRECISION :: ltau_w_g
!---------------------------------------------------------------------//

! Set reference phase to gas
      M = 0

      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         ! Original

          IJK = FUNIJK(i,j,k)

          IJKT = FUNIJK(i,j,ktop(i,j,k))
          IMJK = FUNIJK(iminus(i,j,k),j,k)
          IPJK = FUNIJK(iplus(i,j,k),j,k)
          IJMK = FUNIJK(i,jminus(i,j,k),k)
          IJPK = FUNIJK(i,jplus(i,j,k),k)
          IJKM = FUNIJK(i,j,kminus(i,j,k))
          IJKP = FUNIJK(i,j,kplus(i,j,k))

          IMJKP = FUNIJK(iminus(i,j,k),j,KPLUS(iminus(i,j,k),j,k))
          IJMKP = FUNIJK(i,jminus(i,j,k),KPLUS(i,jminus(i,j,k),k))

         EPGA = AVG(EP_G(IJK),EP_G(IJKT))

! Impermeable internal surface
         IF (ip_at_t(i,j,k)) THEN
            A_M(I,J,K,E) = ZERO
            A_M(I,J,K,W) = ZERO
            A_M(I,J,K,N) = ZERO
            A_M(I,J,K,S) = ZERO
            A_M(I,J,K,T) = ZERO
            A_M(I,J,K,B) = ZERO
            A_M(I,J,K,0) = -ONE
            B_M(I,J,K) = ZERO

! Dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_M(I,J,K,E) = ZERO
            A_M(I,J,K,W) = ZERO
            A_M(I,J,K,N) = ZERO
            A_M(I,J,K,S) = ZERO
            A_M(I,J,K,T) = ZERO
            A_M(I,J,K,B) = ZERO
            A_M(I,J,K,0) = -ONE
            B_M(I,J,K) = ZERO
! set velocity equal to that of bottom or top cell if solids are
! present in those cells else set velocity equal to known value
            IF (EP_G(FUNIJK(i,j,kbot(i,j,k))) > DIL_EP_S) THEN
               A_M(I,J,K,B) = ONE
            ELSE IF (EP_G(FUNIJK(i,j,ktop(i,j,k))) > DIL_EP_S) THEN
               A_M(I,J,K,T) = ONE
            ELSE
               B_M(I,J,K) = -W_G(IJK)
            ENDIF

! Normal case
         ELSE

! Surface forces

! Pressure term
            PGT = P_G(IJKT)
            IF (CYCLIC_Z_PD) THEN
               IF (KMAP(K).EQ.KMAX1) PGT = P_G(IJKT) - DELP_Z
            ENDIF
            SDP = -P_SCALE*EPGA*(PGT - P_G(IJK))*AXY

! Volumetric forces
            ROPGA = AVG(ROP_G(IJK),ROP_G(IJKT))
            ROGA = AVG(RO_G(IJK),RO_G(IJKT))

! Previous time step
            V0 = AVG(ROP_GO(I,J,K),ROP_GO(i,j,ktop(i,j,k)))*ODT

! Body force
            VBF = ROPGA*GRAVITY_Z


            ltau_w_g = tau_w_g(ijk)

! Collect the terms

            A_M(I,J,K,0) = -(A_M(I,J,K,E)+A_M(I,J,K,W)+&
               A_M(I,J,K,N)+A_M(I,J,K,S)+A_M(I,J,K,T)+A_M(I,J,K,B)+&
               V0*VOL)

            B_M(I,J,K) = B_M(I,J,K) - ( SDP + lTAU_W_G  + &
               ( (V0)*W_GO(I,J,K) + VBF)*VOL)

         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)
      ENDDO   ! end do loop over ijk
      ENDDO   ! end do loop over ijk
      ENDDO   ! end do loop over ijk

! modifications for bc
      CALL SOURCE_W_G_BC (A_M, B_M)

      RETURN
      END SUBROUTINE SOURCE_W_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_W_g_BC                                           C
!  Purpose: Determine source terms for W_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_W_G_BC(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      use matrix, only: e, w, s, n, t, b
      USE scales
      USE constant
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions, only: ieast, iwest, jnorth, jsouth, kbot, ktop
      USE functions, only: kminus, kplus
      USE functions, only: funijk, fluid_at, wall_at, ns_wall_at, fs_wall_at
      USE functions, only: is_on_mype_plus2layers
      USE functions, only: im1, jm1

      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! C_mu is constant in turbulent viscosity
      DOUBLE PRECISION, PARAMETER :: C_mu = 0.09D0
! Kappa is Von Karmen constant
      DOUBLE PRECISION, PARAMETER :: Kappa = 0.42D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK,&
                 IM, JM, IJKB, IJKM, IJKP
! Phase index
      INTEGER :: M

!-----------------------------------------------

! Set reference phase to gas
      M = 0


! Set the default boundary conditions
! The NS default setting is the where bc_type='dummy' or any default
! (i.e., bc_type=undefined) wall boundary regions are handled. Note that
! the top and bottom xy planes do not have to be explicitly addressed for
! the w-momentum equation. In this direction the velocities are defined
! at the wall (due staggered grid). They are defined as zero for a
! no penetration condition (see zero_norm_vel subroutine and code under
! the ip_at_t branch in the above source routine).
! ---------------------------------------------------------------->>>

! south xz plane
      J1 = 1
      DO K1 = kmin3,kmax3
         DO I1 = imin3,imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (ns_wall_at(i1,j1,k1)) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = -ONE
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (fs_wall_at(i1,j1,k1)) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ONE
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
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
            IF (ns_wall_at(i1,j1,k1)) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and oppostive to the adjacent fluid cell value)
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = -ONE
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (fs_wall_at(i1,j1,k1)) THEN
! Setting the wall velocity equal to the adjacent fluid velocity (set
! the boundary cell value equal to adjacent fluid cell value)
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ONE
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ENDIF
         ENDDO
      ENDDO

! west yz plane
      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (ns_wall_at(i1,j1,k1)) THEN
               A_M(I1,J1,K1,E) = -ONE
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (fs_wall_at(i1,j1,k1)) THEN
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

! east yz plane
      I1 = IMAX2
      DO K1 = kmin3,kmax3
         DO J1 = jmin3,jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (ns_wall_at(i1,j1,k1)) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = -ONE
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (fs_wall_at(i1,j1,k1)) THEN
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
                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (.NOT.wall_at(i,j,k)) CYCLE  !skip redefined cells
                        A_M(I,J,K,E) = ZERO
                        A_M(I,J,K,W) = ZERO
                        A_M(I,J,K,N) = ZERO
                        A_M(I,J,K,S) = ZERO
                        A_M(I,J,K,T) = ZERO
                        A_M(I,J,K,B) = ZERO
                        A_M(I,J,K,0) = -ONE
                        B_M(I,J,K) = ZERO
                        if (fluid_at(ieast(i,j,k),j,k)) THEN
                           A_M(I,J,K,E) = -ONE
                        else if (fluid_at(iwest(i,j,k),j,k)) THEN
                           A_M(I,J,K,W) = -ONE
                        else if (fluid_at(i,jnorth(i,j,k),k)) THEN
                           A_M(I,J,K,N) = -ONE
                        else if (fluid_at(i,jsouth(i,j,k),k)) THEN
                           A_M(I,J,K,S) = -ONE
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
                        IF (.NOT.wall_at(i,j,k)) CYCLE  !skip redefined cells
                        A_M(I,J,K,E) = ZERO
                        A_M(I,J,K,W) = ZERO
                        A_M(I,J,K,N) = ZERO
                        A_M(I,J,K,S) = ZERO
                        A_M(I,J,K,T) = ZERO
                        A_M(I,J,K,B) = ZERO
                        A_M(I,J,K,0) = -ONE
                        B_M(I,J,K) = ZERO
                        if (fluid_at(ieast(i,j,k),j,k)) THEN
                           A_M(I,J,K,E) = ONE
                        else if (fluid_at(iwest(i,j,k),j,k)) THEN
                           A_M(I,J,K,W) = ONE
                        else if (fluid_at(i,jnorth(i,j,k),k)) THEN
                           A_M(I,J,K,N) = ONE
                        else if (fluid_at(i,jsouth(i,j,k),k)) THEN
                           A_M(I,J,K,S) = ONE
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
                        IF (.NOT.wall_at(i,j,k)) CYCLE  ! skip redefined cells
                        IM = IM1(I)
                        JM = JM1(J)
                        A_M(I,J,K,E) = ZERO
                        A_M(I,J,K,W) = ZERO
                        A_M(I,J,K,N) = ZERO
                        A_M(I,J,K,S) = ZERO
                        A_M(I,J,K,T) = ZERO
                        A_M(I,J,K,B) = ZERO
                        A_M(I,J,K,0) = -ONE
                        B_M(I,J,K) = ZERO
                        if (fluid_at(ieast(i,j,k),j,k)) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(I,J,K,E) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_WW_G(L)
                           ELSE
                                 A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODX)
                                 A_M(I,J,K,E) = -(HALF*BC_HW_G(L)-ODX)
                                 B_M(I,J,K) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        else if (fluid_at(iwest(i,j,k),j,k)) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(I,J,K,W) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_WW_G(L)
                           ELSE
                                 A_M(I,J,K,W) = -(HALF*BC_HW_G(L)-ODX)
                                 A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODX)
                                 B_M(I,J,K) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        else if (fluid_at(i,jnorth(i,j,k),k)) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(I,J,K,N) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_WW_G(L)
                           ELSE
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODY)
                              A_M(I,J,K,N) = -(HALF*BC_HW_G(L)-ODY)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        else if (fluid_at(i,jsouth(i,j,k),k)) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(I,J,K,S) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_WW_G(L)
                           ELSE
                              A_M(I,J,K,S) = -(HALF*BC_HW_G(L)-ODY)
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODY)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO


! end setting of wall boundary conditions
! ----------------------------------------------------------------<<<

! Setting p_inflow or p_outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN
               IF (BC_PLANE(L) == 'B') THEN
! if the fluid cell is on the bottom side of the outflow/inflow boundary
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
                           A_M(I,J,K,E) = ZERO
                           A_M(I,J,K,W) = ZERO
                           A_M(I,J,K,N) = ZERO
                           A_M(I,J,K,S) = ZERO
                           A_M(I,J,K,T) = ZERO
                           A_M(I,J,K,B) = ONE
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
               IF (BC_PLANE(L) == 'B') THEN
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
                           A_M(I,J,K,E) = ZERO
                           A_M(I,J,K,W) = ZERO
                           A_M(I,J,K,N) = ZERO
                           A_M(I,J,K,S) = ZERO
                           A_M(I,J,K,T) = ZERO
                           A_M(I,J,K,B) = ONE
                           A_M(I,J,K,0) = -ONE
                           B_M(I,J,K) = ZERO

                           A_M(I,J,kminus(i,j,k),E) = ZERO
                           A_M(I,J,kminus(i,j,k),W) = ZERO
                           A_M(I,J,kminus(i,j,k),N) = ZERO
                           A_M(I,J,kminus(i,j,k),S) = ZERO
                           A_M(I,J,kminus(i,j,k),T) = ZERO
                           A_M(I,J,kminus(i,j,k),B) = ONE
                           A_M(I,J,kminus(i,j,k),0) = -ONE
                           B_M(I,J,kminus(i,j,k)) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ELSEIF (BC_PLANE(L) == 'T') THEN
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

                           A_M(i,j,kplus(i,j,k),E) = ZERO
                           A_M(i,j,kplus(i,j,k),W) = ZERO
                           A_M(i,j,kplus(i,j,k),N) = ZERO
                           A_M(i,j,kplus(i,j,k),S) = ZERO
                           A_M(i,j,kplus(i,j,k),T) = ONE
                           A_M(i,j,kplus(i,j,k),B) = ZERO
                           A_M(i,j,kplus(i,j,k),0) = -ONE
                           B_M(i,j,kplus(i,j,k)) = ZERO
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
                        A_M(I,J,K,E) = ZERO
                        A_M(I,J,K,W) = ZERO
                        A_M(I,J,K,N) = ZERO
                        A_M(I,J,K,S) = ZERO
                        A_M(I,J,K,T) = ZERO
                        A_M(I,J,K,B) = ZERO
                        A_M(I,J,K,0) = -ONE
                        B_M(I,J,K) = -W_G(IJK)
                        IF (BC_PLANE(L) == 'B') THEN
! if the fluid cell is on the bottom side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           A_M(i,j,kbot(i,j,k),E) = ZERO
                           A_M(i,j,kbot(i,j,k),W) = ZERO
                           A_M(i,j,kbot(i,j,k),N) = ZERO
                           A_M(i,j,kbot(i,j,k),S) = ZERO
                           A_M(i,j,kbot(i,j,k),T) = ZERO
                           A_M(i,j,kbot(i,j,k),B) = ZERO
                           A_M(i,j,kbot(i,j,k),0) = -ONE
                           B_M(i,j,kbot(i,j,k)) = -W_G(FUNIJK(i,j,kbot(i,j,k)))
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
      SUBROUTINE POINT_SOURCE_W_G(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar
      use constant
      use geometry
      use param1, only: small_number, zero
      use fldvar
      use ps
      use run
      use functions
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV, M
      INTEGER :: lKT, lKB
! terms of bm expression
      DOUBLE PRECISION :: pSource
!-----------------------------------------------

! Set reference phase to gas
      M = 0

! Calculate the mass going into each IJK cell. This is done for each
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

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(.NOT.fluid_at(i,j,k)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

            B_M(I,J,K) = B_M(I,J,K) - pSource * &
               PS_W_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_W_G
