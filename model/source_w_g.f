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

!     USE compar, only: ijkstart3, ijkend3, kmap
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3, kmap
      USE compar, only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar, only: istart1, iend1, jstart1, jend1, kstart1, kend1

      USE fldvar, only: p_g, ro_g, rop_g, rop_go
      USE fldvar, only: ep_g
      USE fldvar, only: w_g, w_go

      USE fun_avg, only: avg_x, avg_z, avg_y
      USE fun_avg, only: avg_x_h, avg_z_h
      USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      USE functions, only: ip_at_t, sip_at_t, is_id_at_t
      USE functions, only: ip_of, jp_of, kp_of, im_of, jm_of, km_of
      USE functions, only: east_of, west_of, top_of, bottom_of
      USE functions, only: iminus,iplus,jminus,jplus,kminus,kplus,new_top_of
      USE functions, only: zmax, funijk, wall_at
      USE geometry, only: kmax1, cyclic_z_pd
      USE geometry, only: vol, vol_w
      USE geometry, only: axy

      USE indices, only: i_of, j_of, k_of
      use matrix, only: e, w, s, n, t, b

      USE param, only: dimension_3
      USE param1, only: zero, one, half
      USE run, only: momentum_z_eq
      USE run, only: odt
      USE scales, only: p_scale
      USE fldvar, only: tau_w_g
      USE toleranc, only: dil_ep_s
      USE cutcell, only: cartesian_grid, cut_w_treatment_at
      USE cutcell, only: blocked_w_cell_at
      USE cutcell, only: a_wpg_t, a_wpg_b
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK, IJKT, IMJK, IJKP, IMJKP,&
                 IJKE, IJKW, IJKTE, IJKTW, IM, IPJK,   &
                 IJKM, IJMK, IJMKP, IJPK
! Phase index
      INTEGER :: M, L, MM
! Internal surface
      INTEGER :: ISV
! Pressure at top cell
      DOUBLE PRECISION :: PgT
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average density
      DOUBLE PRECISION :: ROPGA, ROGA
! Average viscosity
      DOUBLE PRECISION :: MUGA
! Average coefficients
      DOUBLE PRECISION ::  MUoX
! Average U_g
      DOUBLE PRECISION Ugt
! Source terms (Surface)
      DOUBLE PRECISION Sdp
! Source terms (Volumetric)
      DOUBLE PRECISION V0, Vbf
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION HYS_drag, avgDrag
! jackson terms: local stress tensor quantity
      DOUBLE PRECISION :: ltau_w_g
!---------------------------------------------------------------------//

     integer err
     integer new_kt, new_ip, new_im, new_jp, new_jm, new_kp, new_km
     integer new_ijkt, new_imjk,new_ipjk,new_ijmk,new_ijpk,new_ijkm, new_ijkp, new_imjkp,new_ijmkp


! Set reference phase to gas
      M = 0

      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         ! Original
         ! I = I_OF(IJK)
         ! J = J_OF(IJK)
         ! K = K_OF(IJK)

          IJK = FUNIJK(i,j,k)

         IJKT = TOP_OF(IJK)

         IMJK = IM_OF(IJK)
         IPJK = IP_OF(IJK)

         IJMK = JM_OF(IJK)
         IJPK = JP_OF(IJK)

         IJKM = KM_OF(IJK)
         IJKP = KP_OF(IJK)

         IMJKP = KP_OF(IMJK)
         IJMKP = KP_OF(IJMK)

         ! New
          New_KT = NEW_TOP_OF(I,j,k)
          New_IJKT = FUNIJK(i,j,New_KT)

          New_IM = IMINUS(I,J,K)
          New_IMJK = FUNIJK(NEW_IM,j,k)

          New_IP = IPLUS(I,J,K)
          New_IPJK = FUNIJK(NEW_IP,j,k)

          New_JM = JMINUS(I,J,K)
          New_IJMK = FUNIJK(i,NEW_JM,k)

          New_JP = JPLUS(I,J,K)
          New_IJPK = FUNIJK(i,NEW_JP,k)

          New_KM = KMINUS(I,J,K)
          New_IJKM = FUNIJK(i,j,New_KM)

          New_KP = KPLUS(I,J,K)
          New_IJKP = FUNIJK(i,j,New_KP)

          New_IMJKP = FUNIJK(New_IM,j,KPLUS(New_IM,j,k))
          New_IJMKP = FUNIJK(i,New_JM,KPLUS(i,New_JM,k))

         if (.not. WALL_AT(ijk)) then
          err = 0
          err = max(abs(ijkt-new_ijkt),err)
          err = max(abs(ijkm-new_ijkm),err)
          err = max(abs(ipjk-new_ipjk),err)
          err = max(abs(imjk-new_imjk),err)
          err = max(abs(ijpk-new_ijpk),err)
          err = max(abs(ijmk-new_ijmk),err)
          err = max(abs(ijkp-new_ijkp),err)
          err = max(abs(ijkm-new_ijkm),err)
          err = max(abs(imjkp-new_imjkp),err)
          err = max(abs(ijmkp-new_ijmkp),err)

          if(err /= 0) then
             write(*,*)'ERR      AT I,j,k        ' ,i,j,FUNIJK(i,j,k)
             write(*,*)'ERR ',i,j,k,err
             stop
          endif
          endif

         EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)

! Impermeable internal surface
         IF (IP_AT_T(IJK)) THEN
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
! set velocity equal to that of bottom or top cell if solids are
! present in those cells else set velocity equal to known value
            IF (EP_G(BOTTOM_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,B) = ONE
            ELSE IF (EP_G(TOP_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,T) = ONE
            ELSE
               B_M(IJK) = -W_G(IJK)
            ENDIF

! Cartesian grid implementation
         ELSEIF (BLOCKED_W_CELL_AT(IJK)) THEN
            A_M(IJK,E) = ZERO
            A_M(IJK,W) = ZERO
            A_M(IJK,N) = ZERO
            A_M(IJK,S) = ZERO
            A_M(IJK,T) = ZERO
            A_M(IJK,B) = ZERO
            A_M(IJK,0) = -ONE
            B_M(IJK) = ZERO

! Normal case
         ELSE

! Surface forces

! Pressure term
            PGT = P_G(IJKT)
            IF (CYCLIC_Z_PD) THEN
               IF (KMAP(K_OF(IJK)).EQ.KMAX1) PGT = P_G(IJKT) - DELP_Z
            ENDIF
            IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                SDP = -P_SCALE*EPGA*(PGT - P_G(IJK))*AXY(IJK)
            ELSE
                SDP = -P_SCALE*EPGA*(PGT * A_WPG_T(IJK) - P_G(IJK) * A_WPG_B(IJK) )
            ENDIF

            IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
! Volumetric forces
               ROPGA = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K)
               ROGA = AVG_Z(RO_G(IJK),RO_G(IJKT),K)

! Previous time step
               V0 = AVG_Z(ROP_GO(IJK),ROP_GO(IJKT),K)*ODT

            ELSE
! Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + VOL(IJKT)*ROP_G(IJKT))/&
                  (VOL(IJK) + VOL(IJKT))
               ROGA  = (VOL(IJK)*RO_G(IJK) + VOL(IJKT)*RO_G(IJKT))/&
                  (VOL(IJK) + VOL(IJKT))
! Previous time step
               V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IJKT)*ROP_GO(IJKT))*&
                  ODT/(VOL(IJK) + VOL(IJKT))
            ENDIF


! Body force
            VBF = ROPGA*GRAVITY_Z


            ltau_w_g = tau_w_g(ijk)

! Collect the terms

            A_M(IJK,0) = -(A_M(IJK,E)+A_M(IJK,W)+&
               A_M(IJK,N)+A_M(IJK,S)+A_M(IJK,T)+A_M(IJK,B)+&
               V0*VOL_W(IJK))

            B_M(IJK) = B_M(IJK) - ( SDP + lTAU_W_G  + &
               ( (V0)*W_GO(IJK) + VBF)*VOL_W(IJK) )

         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)
      ENDDO   ! end do loop over ijk
      ENDDO   ! end do loop over ijk
      ENDDO   ! end do loop over ijk

! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_W_G(A_M, B_M)
! modifications for bc
      CALL SOURCE_W_G_BC (A_M, B_M)
! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_W_G_BC(A_M, B_M)

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
      USE indices
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
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
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = -ONE
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
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
            IF (NS_WALL_AT(IJK)) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and oppostive to the adjacent fluid cell value)
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = -ONE
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
! Setting the wall velocity equal to the adjacent fluid velocity (set
! the boundary cell value equal to adjacent fluid cell value)
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

! west yz plane
      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E) = -ONE
               A_M(IJK,W) = ZERO
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
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

! east yz plane
      I1 = IMAX2
      DO K1 = kmin3,kmax3
         DO J1 = jmin3,jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E) = ZERO
               A_M(IJK,W) = -ONE
               A_M(IJK,N) = ZERO
               A_M(IJK,S) = ZERO
               A_M(IJK,T) = ZERO
               A_M(IJK,B) = ZERO
               A_M(IJK,0) = -ONE
               B_M(IJK) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
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
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        A_M(IJK,E) = ZERO
                        A_M(IJK,W) = ZERO
                        A_M(IJK,N) = ZERO
                        A_M(IJK,S) = ZERO
                        A_M(IJK,T) = ZERO
                        A_M(IJK,B) = ZERO
                        A_M(IJK,0) = -ONE
                        B_M(IJK) = ZERO
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           A_M(IJK,E) = -ONE
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           A_M(IJK,W) = -ONE
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           A_M(IJK,N) = -ONE
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           A_M(IJK,S) = -ONE
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
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        A_M(IJK,E) = ZERO
                        A_M(IJK,W) = ZERO
                        A_M(IJK,N) = ZERO
                        A_M(IJK,S) = ZERO
                        A_M(IJK,T) = ZERO
                        A_M(IJK,B) = ZERO
                        A_M(IJK,0) = -ONE
                        B_M(IJK) = ZERO
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           A_M(IJK,E) = ONE
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           A_M(IJK,W) = ONE
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           A_M(IJK,N) = ONE
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           A_M(IJK,S) = ONE
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
                        JM = JM1(J)
                        A_M(IJK,E) = ZERO
                        A_M(IJK,W) = ZERO
                        A_M(IJK,N) = ZERO
                        A_M(IJK,S) = ZERO
                        A_M(IJK,T) = ZERO
                        A_M(IJK,B) = ZERO
                        A_M(IJK,0) = -ONE
                        B_M(IJK) = ZERO
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,E) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_WW_G(L)
                           ELSE
                                 A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODX_E(I))
                                 A_M(IJK,E) = -(HALF*BC_HW_G(L)-ODX_E(I))
                                 B_M(IJK) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,W) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_WW_G(L)
                           ELSE
                                 A_M(IJK,W) = -(HALF*BC_HW_G(L)-ODX_E(IM))
                                 A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODX_E(IM))
                                 B_M(IJK) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,N) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_WW_G(L)
                           ELSE
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODY_N(J))
                              A_M(IJK,N) = -(HALF*BC_HW_G(L)-ODY_N(J))
                              B_M(IJK) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,S) = -HALF
                              A_M(IJK,0) = -HALF
                              B_M(IJK) = -BC_WW_G(L)
                           ELSE
                              A_M(IJK,S) = -(HALF*BC_HW_G(L)-ODY_N(JM))
                              A_M(IJK,0) = -(HALF*BC_HW_G(L)+ODY_N(JM))
                              B_M(IJK) = -BC_HW_G(L)*BC_WW_G(L)
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
                           A_M(IJK,E) = ZERO
                           A_M(IJK,W) = ZERO
                           A_M(IJK,N) = ZERO
                           A_M(IJK,S) = ZERO
                           A_M(IJK,T) = ZERO
                           A_M(IJK,B) = ONE
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
                           A_M(IJK,E) = ZERO
                           A_M(IJK,W) = ZERO
                           A_M(IJK,N) = ZERO
                           A_M(IJK,S) = ZERO
                           A_M(IJK,T) = ZERO
                           A_M(IJK,B) = ONE
                           A_M(IJK,0) = -ONE
                           B_M(IJK) = ZERO
                           IJKM = KM_OF(IJK)
                           A_M(IJKM,E) = ZERO
                           A_M(IJKM,W) = ZERO
                           A_M(IJKM,N) = ZERO
                           A_M(IJKM,S) = ZERO
                           A_M(IJKM,T) = ZERO
                           A_M(IJKM,B) = ONE
                           A_M(IJKM,0) = -ONE
                           B_M(IJKM) = ZERO
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
                           IJK = FUNIJK(I,J,K)
                           IJKP = KP_OF(IJK)
                           A_M(IJKP,E) = ZERO
                           A_M(IJKP,W) = ZERO
                           A_M(IJKP,N) = ZERO
                           A_M(IJKP,S) = ZERO
                           A_M(IJKP,T) = ONE
                           A_M(IJKP,B) = ZERO
                           A_M(IJKP,0) = -ONE
                           B_M(IJKP) = ZERO
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
                        B_M(IJK) = -W_G(IJK)
                        IF (BC_PLANE(L) == 'B') THEN
! if the fluid cell is on the bottom side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           IJKB = BOTTOM_OF(IJK)
                           A_M(IJKB,E) = ZERO
                           A_M(IJKB,W) = ZERO
                           A_M(IJKB,N) = ZERO
                           A_M(IJKB,S) = ZERO
                           A_M(IJKB,T) = ZERO
                           A_M(IJKB,B) = ZERO
                           A_M(IJKB,0) = -ONE
                           B_M(IJKB) = -W_G(IJKB)
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
      use indices
      use param1, only: small_number, zero
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
            if(.NOT.fluid_at(ijk)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

            B_M(IJK) = B_M(IJK) - pSource * &
               PS_W_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_W_G
