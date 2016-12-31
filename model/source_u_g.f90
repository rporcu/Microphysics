module source_u_g_module

   use param1, only: zero, half, one, undefined, is_undefined, small_number

  contains
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
      SUBROUTINE SOURCE_U_G(A_M, B_M, dt, p_g, ep_g, ro_g, rop_g, rop_go, &
         u_g, u_go, tau_u_g, flag)

! Modules
!---------------------------------------------------------------------//
      use compar  , only: imap
      USE compar  , only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE constant, only: gravity
      USE bc      , only: delp_x

      USE functions, only: avg
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus,ieast,iwest
      USE functions, only: zmax

      USE geometry, only: imax1, cyclic_x_pd
      USE geometry, only: vol
      USE geometry, only: ayz

      use matrix, only: e, w, s, n, t, b

      USE scales, only: p_scale
      USE toleranc, only: dil_ep_s

      IMPLICIT NONE

      double precision, intent(in   ) :: dt

      double precision, intent(in   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: rop_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: u_go&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) :: tau_u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer, intent(in   ) :: flag &
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

! Dummy arguments
!--------o------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: i,j,k
! Phase index
      INTEGER :: m
! Pressure at east cell
      DOUBLE PRECISION :: PgE
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average density
      DOUBLE PRECISION :: ROPGA, ROGA
! Source terms (Surface)
      DOUBLE PRECISION :: Sdp
! Source terms (Volumetric)
      DOUBLE PRECISION :: V0, Vbf
! local stress tensor quantity
      DOUBLE PRECISION :: ltau_u_g
      double precision :: odt
!---------------------------------------------------------------------//

      odt = 1.0d0/dt

! Set reference phase to gas
      M = 0

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

          EPGA = AVG(EP_G(I,J,K),EP_G(ieast(i,j,k),j,k))

! Impermeable internal surface
         IF (flag(i,j,k,2) < 1000) THEN
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
! set velocity equal to that of west or east cell if solids are present
! in those cells else set velocity equal to known value
            IF (EP_G(iwest(i,j,k),j,k) > DIL_EP_S) THEN
               A_M(I,J,K,W) = ONE
            ELSE IF (EP_G(ieast(i,j,k),j,k) > DIL_EP_S) THEN
               A_M(I,J,K,E) = ONE
            ELSE
               B_M(I,J,K) = -U_G(I,J,K)
            ENDIF

! Normal case
         ELSE

! Surface forces
! Pressure term
            PGE = P_G(ieast(i,j,k),j,k)
            IF (CYCLIC_X_PD) THEN
               IF (IMAP(I).EQ.IMAX1) PGE = P_G(ieast(i,j,k),j,k) - DELP_X
            ENDIF
            SDP = -P_SCALE*EPGA*(PGE - P_G(I,J,K))*AYZ

! Volumetric forces
            ROGA  = HALF * (RO_G(I,J,K) + RO_G(ieast(i,j,k),j,k))
            ROPGA = HALF * (ROP_G(I,J,K) + ROP_G(ieast(i,j,k),j,k))
! Previous time step
            V0 = HALF * (ROP_GO(I,J,K) + ROP_GO(ieast(i,j,k),j,k))*ODT

! Body force
            VBF = ROGA*GRAVITY(1)

            ltau_u_g = tau_u_g(i,j,k)

! Collect the terms
            A_M(I,J,K,0) = -(A_M(I,J,K,E)+A_M(I,J,K,W)+&
               A_M(I,J,K,N)+A_M(I,J,K,S)+A_M(I,J,K,T)+A_M(I,J,K,B)+&
               V0*VOL)

            B_M(I,J,K) = B_M(I,J,K) -(SDP + lTAU_U_G + &
               ( (V0)*U_GO(I,J,K) + VBF)*VOL )

         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)

      ENDDO   ! end do loop over ijk
      ENDDO   ! end do loop over ijk
      ENDDO   ! end do loop over ijk

! modifications for bc
      CALL SOURCE_U_G_BC (A_M, B_M, U_G, flag)

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

      SUBROUTINE SOURCE_U_G_BC(A_M, B_M, U_G, flag)

      USE bc, only: bc_hw_g, bc_uw_g
      USE bc, only: bc_i_w, bc_i_e, bc_j_s, bc_j_n, bc_k_b, bc_k_t
      USE bc, only: dimension_bc, bc_type, bc_defined, bc_plane
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE functions, only: ieast, iwest, jsouth, jnorth, kbot, ktop
      USE functions, only: iminus, iplus, im1
      USE geometry  , only: imin3, imax3, jmin3, jmax3, kmin3, kmax3
      USE geometry  , only: jmax2, kmax2
      USE geometry  , only: ody, odz
      USE matrix, only: e, w, s, n, t, b

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Velocity u_g
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER ::  I,  J, K, IM, I1, I2, J1, J2, K1, K2
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
! no penetration condition.
! ---------------------------------------------------------------->>>
! bottom xy plane
      K1 = 1
      DO J1 = jmin3,jmax3
         DO I1 = imin3, imax3
            IF (flag(i1,j1,k1,1) == 100) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and opposite to the adjacent fluid cell value)
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

! south xz plane
      J1 = 1
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (flag(i1,j1,k1,1) == 100) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = -ONE
               A_M(I1,J1,K1,S) = ZERO
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (flag(i1,j1,k1,1) == 101) THEN
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
            IF (flag(i1,j1,k1,1) == 100) THEN
               A_M(I1,J1,K1,E) = ZERO
               A_M(I1,J1,K1,W) = ZERO
               A_M(I1,J1,K1,N) = ZERO
               A_M(I1,J1,K1,S) = -ONE
               A_M(I1,J1,K1,T) = ZERO
               A_M(I1,J1,K1,B) = ZERO
               A_M(I1,J1,K1,0) = -ONE
               B_M(I1,J1,K1) = ZERO
            ELSEIF (flag(i1,j1,k1,1) == 101) THEN
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
                        if (1.eq.flag(i,jnorth(i,j,k),k,1)) THEN
                           A_M(I,J,K,N) = -ONE
                        else iF (1.eq.flag(i,jsouth(i,j,k),k,1)) THEN
                           A_M(I,J,K,S) = -ONE
                        else iF (1.eq.flag(i,j,ktop(i,j,k),1)) THEN
                           A_M(I,J,K,T) = -ONE
                        else iF (1.eq.flag(i,j,kbot(i,j,k),1)) THEN
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
                        if (1.eq.flag(i,jnorth(i,j,k),k,1)) THEN
                           A_M(I,J,K,N) = ONE
                        else if (1.eq.flag(i,jsouth(i,j,k),k,1)) THEN
                           A_M(I,J,K,S) = ONE
                        else if (1.eq.flag(i,j,ktop(i,j,k),1)) THEN
                           A_M(I,J,K,T) = ONE
                        else if (1.eq.flag(i,j,kbot(i,j,k),1)) THEN
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
                        if (1.eq.flag(i,jnorth(i,j,k),k,1)) THEN
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_M(I,J,K,N) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_UW_G(L)
                           ELSE
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODY)
                              A_M(I,J,K,N) = -(HALF*BC_HW_G(L)-ODY)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        else if (1.eq.flag(i,jsouth(i,j,k),k,1)) THEN
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_M(I,J,K,S) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_UW_G(L)
                           ELSE
                              A_M(I,J,K,S) = -(HALF*BC_HW_G(L)-ODY)
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODY)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        else if (1.eq.flag(i,j,ktop(i,j,k),1)) THEN
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_M(I,J,K,T) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_UW_G(L)
                           ELSE
                              A_M(I,J,K,0)=-(HALF*BC_HW_G(L)+ODZ)
                              A_M(I,J,K,T)=-(HALF*BC_HW_G(L)-ODZ)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        else if (1.eq.flag(i,j,kbot(i,j,k),1)) THEN
                           IF (IS_UNDEFINED(BC_HW_G(L))) THEN
                              A_M(I,J,K,B) = -HALF
                              A_M(I,J,K,0) = -HALF
                              B_M(I,J,K) = -BC_UW_G(L)
                           ELSE
                              A_M(I,J,K,B) = -(HALF*BC_HW_G(L)-ODZ)
                              A_M(I,J,K,0) = -(HALF*BC_HW_G(L)+ODZ)
                              B_M(I,J,K) = -BC_HW_G(L)*BC_UW_G(L)
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
                           A_M(I,J,K,E) = ZERO
                           A_M(I,J,K,W) = ONE
                           A_M(I,J,K,N) = ZERO
                           A_M(I,J,K,S) = ZERO
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
                           A_M(I,J,K,E) = ZERO
                           A_M(I,J,K,W) = ONE
                           A_M(I,J,K,N) = ZERO
                           A_M(I,J,K,S) = ZERO
                           A_M(I,J,K,T) = ZERO
                           A_M(I,J,K,B) = ZERO
                           A_M(I,J,K,0) = -ONE
                           B_M(I,J,K) = ZERO
                           IM = IM1(I)
                           A_M(iminus(i,j,k),j,k,E) = ZERO
                           A_M(iminus(i,j,k),j,k,W) = ONE
                           A_M(iminus(i,j,k),j,k,N) = ZERO
                           A_M(iminus(i,j,k),j,k,S) = ZERO
                           A_M(iminus(i,j,k),j,k,T) = ZERO
                           A_M(iminus(i,j,k),j,k,B) = ZERO
                           A_M(iminus(i,j,k),j,k,0) = -ONE
                           B_M(iminus(i,j,k),j,k) = ZERO
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
                           A_M(iplus(i,j,k),j,k,E) = ONE
                           A_M(iplus(i,j,k),j,k,W) = ZERO
                           A_M(iplus(i,j,k),j,k,N) = ZERO
                           A_M(iplus(i,j,k),j,k,S) = ZERO
                           A_M(iplus(i,j,k),j,k,T) = ZERO
                           A_M(iplus(i,j,k),j,k,B) = ZERO
                           A_M(iplus(i,j,k),j,k,0) = -ONE
                           B_M(iplus(i,j,k),j,k) = ZERO
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
                        B_M(I,J,K) = -U_G(I,J,K)
                        IF (BC_PLANE(L) == 'W') THEN
! if the fluid cell is on the west side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           A_M(iwest(i,j,k),j,k,E) = ZERO
                           A_M(iwest(i,j,k),j,k,W) = ZERO
                           A_M(iwest(i,j,k),j,k,N) = ZERO
                           A_M(iwest(i,j,k),j,k,S) = ZERO
                           A_M(iwest(i,j,k),j,k,T) = ZERO
                           A_M(iwest(i,j,k),j,k,B) = ZERO
                           A_M(iwest(i,j,k),j,k,0) = -ONE
                           B_M(iwest(i,j,k),j,k) = -U_G(iwest(i,j,k),j,k)
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
      SUBROUTINE POINT_SOURCE_U_G(A_M, B_M, flag)

      use compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use geometry , only: vol
      use ps, only: dimension_ps, ps_defined, ps_volume, ps_vel_mag_g, ps_massflow_g
      use ps, only: ps_u_g, ps_i_e, ps_i_w, ps_j_s, ps_j_n, ps_k_b, ps_k_t

      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN   ) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      integer, intent(in   ) :: flag &
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K
      INTEGER :: PSV, M
      INTEGER :: lIE, lIW
! terms of bm expression
      DOUBLE PRECISION :: pSource
!-----------------------------------------------

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

            pSource =  PS_MASSFLOW_G(PSV) * (VOL/PS_VOLUME(PSV))

            B_M(I,J,K) = B_M(I,J,K) - pSource *                        &
               PS_U_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo
      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_U_G
end module source_u_g_module
