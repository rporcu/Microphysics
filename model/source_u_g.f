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
      USE constant, only: gravity_x
      USE bc, only: delp_x

      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3, imap
      USE compar, only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar, only: istart1, iend1, jstart1, jend1, kstart1, kend1

      USE fldvar, only: p_g, ro_g, rop_g, rop_go
      USE fldvar, only: ep_g
      USE fldvar, only: u_g, u_go

      USE fun_avg, only: avg_x, avg_z, avg_y
      USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      USE functions, only: funijk
      USE functions, only: ip_at_e, sip_at_e, is_id_at_e
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus,ieast,iwest
      USE functions, only: zmax

      USE geometry, only: imax1, cyclic_x_pd, imax
      USE geometry, only: vol
      USE geometry, only: ayz

      use matrix, only: e, w, s, n, t, b

      USE param, only: dimension_3
      USE param1, only: zero, one, half
      USE run, only: momentum_x_eq
      USE run, only: odt
      USE scales, only: p_scale
      USE fldvar, only: tau_u_g
      USE toleranc, only: dil_ep_s
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)

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

      integer err

! Set reference phase to gas
      M = 0

      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

          IJK = FUNIJK(i,j,k)
          IJKE = FUNIJK(ieast(i,j,k),j,k)
          IMJK = FUNIJK(iminus(i,j,k),j,k)
          IPJK = FUNIJK(iplus(i,j,k),j,k)
          IJMK = FUNIJK(i,jminus(i,j,k),k)
          IJPK = FUNIJK(i,jplus(i,j,k),k)
          IJKM = FUNIJK(i,j,kminus(i,j,k))
          IJKP = FUNIJK(i,j,kplus(i,j,k))

          IPJKM = FUNIJK(IPLUS(I,J,kminus(i,j,k)),J,kminus(i,j,k))
          IPJMK = FUNIJK(IPLUS(I,jminus(i,j,k),k),jminus(i,j,k),k)

           EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

! Impermeable internal surface
         IF (ip_at_e(i,j,k)) THEN
            A_M(IJK,E) = ZERO
            A_M(IJK,W) = ZERO
            A_M(IJK,N) = ZERO
            A_M(IJK,S) = ZERO
            A_M(IJK,T) = ZERO
            A_M(IJK,B) = ZERO
            A_M(IJK,0) = -ONE
            B_M(IJK) = ZERO

! Dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_M(IJK,E) = ZERO
            A_M(IJK,W) = ZERO
            A_M(IJK,N) = ZERO
            A_M(IJK,S) = ZERO
            A_M(IJK,T) = ZERO
            A_M(IJK,B) = ZERO
            A_M(IJK,0) = -ONE
            B_M(IJK) = ZERO
! set velocity equal to that of west or east cell if solids are present
! in those cells else set velocity equal to known value
            IF (EP_G(FUNIJK(iwest(i,j,k),j,k)) > DIL_EP_S) THEN
               A_M(IJK,W) = ONE
            ELSE IF (EP_G(FUNIJK(ieast(i,j,k),j,k)) > DIL_EP_S) THEN
               A_M(IJK,E) = ONE
            ELSE
               B_M(IJK) = -U_G(IJK)
            ENDIF

! Normal case
         ELSE

! Surface forces
! Pressure term
            PGE = P_G(IJKE)
            IF (CYCLIC_X_PD) THEN
               IF (IMAP(I).EQ.IMAX1) PGE = P_G(IJKE) - DELP_X
            ENDIF
            SDP = -P_SCALE*EPGA*(PGE - P_G(IJK))*AYZ

! Volumetric forces
            ROPGA = HALF * (ROP_G(IJK) + ROP_G(IJKE))
            ROGA  = HALF * (RO_G(IJK) + RO_G(IJKE))
! Previous time step
            V0 = HALF * (ROP_GO(I,J,K) + ROP_GO(ieast(i,j,k),j,k))*ODT

! Body force
            VBF = ROGA*GRAVITY_X

            ltau_u_g = tau_u_g(ijk)

! Collect the terms
            A_M(IJK,0) = -(A_M(IJK,E)+A_M(IJK,W)+&
               A_M(IJK,N)+A_M(IJK,S)+A_M(IJK,T)+A_M(IJK,B)+&
               V0*VOL)

            B_M(IJK) = B_M(IJK) -(SDP + lTAU_U_G + &
               ( (V0)*U_GO(I,J,K) + VBF)*VOL )

         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)

      ENDDO   ! end do loop over ijk
      ENDDO   ! end do loop over ijk
      ENDDO   ! end do loop over ijk

! modifications for bc
      CALL SOURCE_U_G_BC (A_M, B_M)

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
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions, only: funijk, fs_wall_cell, ns_wall_cell
      USE functions, only: is_on_mype_plus2layers
      USE functions, only: wall_cell, fluid_cell
      USE functions, only: ieast, iwest, jsouth, jnorth, kbot, ktop
      USE functions, only: iminus, iplus
      USE functions, only: im1, ip1, jm1, jp1, km1, kp1
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
               IF (ns_wall_cell(i1,j1,k1)) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and opposite to the adjacent fluid cell value)
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

! south xz plane
      J1 = 1
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (ns_wall_cell(i1,j1,k1)) THEN
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = -ONE
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
            ELSEIF (fs_wall_cell(i1,j1,k1)) THEN
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = ONE
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
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
            IF (ns_wall_cell(i1,j1,k1)) THEN
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = -ONE
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
            ELSEIF (fs_wall_cell(i1,j1,k1)) THEN
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = ONE
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
                        if (fluid_cell(i,jnorth(i,j,k),k)) THEN
                           A_M(IJK,N) = -ONE
                        else iF (fluid_cell(i,jsouth(i,j,k),k)) THEN
                           A_M(IJK,S) = -ONE
                        else iF (fluid_cell(i,j,ktop(i,j,k))) THEN
                           A_M(IJK,T) = -ONE
                        else iF (fluid_cell(i,j,kbot(i,j,k))) THEN
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
                        if (fluid_cell(i,jnorth(i,j,k),k)) THEN
                           A_M(IJK,N) = ONE
                        else if (fluid_cell(i,jsouth(i,j,k),k)) THEN
                           A_M(IJK,S) = ONE
                        else if (fluid_cell(i,j,ktop(i,j,k))) THEN
                           A_M(IJK,T) = ONE
                        else if (fluid_cell(i,j,kbot(i,j,k))) THEN
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
                        JM = JM1(J)
                        KM = KM1(K)
                        A_M(IJK,E) = ZERO
                        A_M(IJK,W) = ZERO
                        A_M(IJK,N) = ZERO
                        A_M(IJK,S) = ZERO
                        A_M(IJK,T) = ZERO
                        A_M(IJK,B) = ZERO
                        A_M(IJK,0) = -ONE
                        B_M(IJK) = ZERO
                        if (fluid_cell(i,jnorth(i,j,k),k)) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,N) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODY_N(J))
                              A_M(IJK,N) = -(HALF*BC_HW_G(L)-ODY_N(J))
                              B_M(IJK) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        else if (fluid_cell(i,jsouth(i,j,k),k)) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,S) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,S) = -(HALF*BC_HW_G(L)-ODY_N(JM))
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODY_N(JM))
                              B_M(IJK) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        else if (fluid_cell(i,j,ktop(i,j,k))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,T) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,0)=-(HALF*BC_HW_G(L)+ODZ_T(K)*OX_E(I))
                              A_M(IJK,T)=-(HALF*BC_HW_G(L)-ODZ_T(K)*OX_E(I))
                              B_M(IJK) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        else if (fluid_cell(i,j,kbot(i,j,k))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,B) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,B) = -(HALF*BC_HW_G(L)-ODZ_T(KM)*OX_E(I&
                                 ))
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODZ_T(KM)*OX_E(I&
                                 ))
                              B_M(IJK) = -BC_HW_G(L)*BC_UW_G(L)
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
                           A_M(IJK,E) = ZERO
                           A_M(IJK,W) = ONE
                           A_M(IJK,N) = ZERO
                           A_M(IJK,S) = ZERO
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
                           A_M(IJK,E) = ZERO
                           A_M(IJK,W) = ONE
                           A_M(IJK,N) = ZERO
                           A_M(IJK,S) = ZERO
                           A_M(IJK,T) = ZERO
                           A_M(IJK,B) = ZERO
                           A_M(IJK,0) = -ONE
                           B_M(IJK) = ZERO
                           IM = IM1(I)
                           IMJK = FUNIJK(iminus(i,j,k),j,k)
                           A_M(IMJK,E) = ZERO
                           A_M(IMJK,W) = X_E(IM)/X_E(IM1(IM))
                           A_M(IMJK,N) = ZERO
                           A_M(IMJK,S) = ZERO
                           A_M(IMJK,T) = ZERO
                           A_M(IMJK,B) = ZERO
                           A_M(IMJK,0) = -ONE
                           B_M(IMJK) = ZERO
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
                           IPJK = FUNIJK(iplus(i,j,k),j,k)
                           A_M(IPJK,E) = X_E(IP)/X_E(I)
                           A_M(IPJK,W) = ZERO
                           A_M(IPJK,N) = ZERO
                           A_M(IPJK,S) = ZERO
                           A_M(IPJK,T) = ZERO
                           A_M(IPJK,B) = ZERO
                           A_M(IPJK,0) = -ONE
                           B_M(IPJK) = ZERO
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
                        B_M(IJK) = -U_G(IJK)
                        IF (BC_PLANE(L) == 'W') THEN
! if the fluid cell is on the west side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           IJKW = FUNIJK(iwest(i,j,k),j,k)
                           A_M(IJKW,E) = ZERO
                           A_M(IJKW,W) = ZERO
                           A_M(IJKW,N) = ZERO
                           A_M(IJKW,S) = ZERO
                           A_M(IJKW,T) = ZERO
                           A_M(IJKW,B) = ZERO
                           A_M(IJKW,0) = -ONE
                           B_M(IJKW) = -U_G(IJKW)
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
            if(.NOT. fluid_cell(i,j,k)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

            B_M(IJK) = B_M(IJK) - pSource *                        &
               PS_U_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo
      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_U_G
