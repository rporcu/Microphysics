!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_W_g                                            C
!  Purpose: Cross terms in the gradient of stress in W_g momentum      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
!                                                                      C
!                                                                      C
!  Comments: This routine calculates the components of the gas phase   C
!  viscous stress tensor of the w-momentum equation that cannot be     C
!  cast in the form: mu.grad(w). These components are stored in the    C
!  passed variable, which is then used as a source of the w-momentum   C
!  equation.                                                           C
!                                                                      C
!  The following w component viscous stress tensor terms are           C
!  calculated here:                                                    C
!  > part of 1/x d/dz (tau_zz) xdxdydz =>                              C
!            1/x d/dz (lambda.trcD) xdxdydz=>                          C
!    delta (lambda.trcD)Ap |T-B                                        C
!  > part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently        C
!    part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>              C
!            1/x d/dx(mu.du/dz) xdxdydz =>                             C
!    delta (mu/x du/dz)Ayz |E-W                                        C
!  > part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently        C
!    part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>              C
!            1/x d/dx(mu.du/dz) xdxdydz =>                             C
!    delta (mu/x du/dz)Ayz |E-W                                        C
!  > part of 1/x d/dz (tau_zz) xdxdydz =>                              C
!            1/x d/dz (mu/x dw/dz) xdxdydz =>                          C
!    delta (mu/x dw/dz)Axy |T-B                                        C
!                                                                      C
!  To reconstitute the complete w-momentum gas phase viscous stress    C
!  tensor requires including the terms calculated in source_w_g        C
!  and the 'diffusional' components (i.e., those of the form           C
!  mu.grad(w)                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_W_G(lTAU_W_G)

! Modules
!---------------------------------------------------------------------//
      USE param
      USE param1
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE sendrecv
      USE compar
      USE functions
      USE cutcell
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_W_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_w_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JM, KP
      INTEGER :: IJK, IJKE, IJKW, IJKN, IJKS, IJKT
      INTEGER :: IJKNT, IJKST, IJKTE, IJKTW
      INTEGER :: IJKP, IMJK, IJMK, IJKM
      INTEGER :: IMJKP, IJMKP
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average gradients
      DOUBLE PRECISION :: duodz
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! Source terms (Volumetric)
      DOUBLE PRECISION :: Vxz
!---------------------------------------------------------------------//

      IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(5)==1)) THEN

!$omp  parallel do default(none) &
!$omp  private(I, J, K, IM, JM, KP, IJK,                               &
!$omp          IJKE, IJKW, IJKN, IJKS, IJKT,                           &
!$omp          IJKTE, IJKNT, IJKTW, IJKST,                             &
!$omp          IMJK, IJMK, IJKM, IJKP, IJMKP, IMJKP,                   &
!$omp          EPGA, SBV, SSX, SSY, SSZ, vxz, duodz)                   &
!$omp  shared(ijkstart3, ijkend3, i_of, j_of, k_of, im1, jm1, kp1,     &
!$omp         do_k, ltau_w_g,              &
!$omp         ep_g, MU_G, LAMBDA_G, trd_g, v_g, w_g, u_g,          &
!$omp         axy,  axy_w, ayz_w, axz_w, vol_w,                        &
!$omp         dy, dz, ox, ox_e, odz, odz_t)
         DO IJK = IJKSTART3, IJKEND3
            K = K_OF(IJK)
            IJKT = TOP_OF(IJK)
            EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)
            IF ( .NOT.IP_AT_T(IJK) .AND. EPGA>DIL_EP_S) THEN
               J = J_OF(IJK)
               I = I_OF(IJK)
               IM = IM1(I)
               JM = JM1(J)
               KP = KP1(K)

               IJKP = KP_OF(IJK)
               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IMJKP = KP_OF(IMJK)
               IJMKP = JM_OF(IJKP)

               IJKN = NORTH_OF(IJK)
               IJKS = SOUTH_OF(IJK)
               IJKE = EAST_OF(IJK)
               IJKW = WEST_OF(IJK)
               IJKNT = TOP_OF(IJKN)
               IJKST = TOP_OF(IJKS)
               IJKTE = EAST_OF(IJKT)
               IJKTW = WEST_OF(IJKT)


! Surface forces
! Bulk viscosity term
! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (lambda.trcD) xdxdydz=>
! delta (lambda.trcD)Ap |T-B : at (i, j, k+1 - k-1)
               SBV = (LAMBDA_G(IJKT)*TRD_G(IJKT)-&
                      LAMBDA_G(IJK)*TRD_G(IJK))*AXY(IJK)

! shear stress terms
! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>
!         1/x d/dx(mu.du/dz) xdxdydz =>
! delta (mu/x du/dz)Ayz |E-W : at (i+1/2-i-1/2, j, k+1/2)
               SSX = AVG_Z_H(AVG_X_H(MU_G(IJK),MU_G(IJKE),I),&
                             AVG_X_H(MU_G(IJKT),MU_G(IJKTE),I),K)*&
                        (U_G(IJKP)-U_G(IJK))*OX_E(I)*ODZ_T(K)*&
                        AYZ_W(IJK) - &
                     AVG_Z_H(AVG_X_H(MU_G(IJKW),MU_G(IJK),IM),&
                             AVG_X_H(MU_G(IJKTW),MU_G(IJKT),IM),K)*&
                        (U_G(IMJKP)-U_G(IMJK))*ODZ_T(K)*DY(J)*&
                        (HALF*(DZ(K)+DZ(KP)))
! DY(J)*HALF(DZ(k)+DZ(kp)) = oX_E(IM)*AYZ_W(IMJK), but avoids singularity

! part of d/dy (tau_zy) xdxdydz =>
!         d/dy (mu/x dv/dz) xdxdydz =>
! delta (mu/x dv/dz)Axz |N-S : at (i, j+1/2 - j-1/2, k+1/2)
               SSY = AVG_Z_H(AVG_Y_H(MU_G(IJK),MU_G(IJKN),J),&
                             AVG_Y_H(MU_G(IJKT),MU_G(IJKNT),J),K)*&
                        (V_G(IJKP)-V_G(IJK))*OX(I)*ODZ_T(K)*AXZ_W(IJK) -&
                     AVG_Z_H(AVG_Y_H(MU_G(IJKS),MU_G(IJK),JM),&
                             AVG_Y_H(MU_G(IJKST),MU_G(IJKT),JM),K)*&
                        (V_G(IJMKP)-V_G(IJMK))*OX(I)*ODZ_T(K)*AXZ_W(IJMK)

! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (mu/x dw/dz) xdxdydz =>
! delta (mu/x dw/dz)Axy |T-B : at (i, j, k+1 - k-1)
               SSZ = MU_G(IJKT)*(W_G(IJKP)-W_G(IJK))*OX(I)*ODZ(KP)*&
                        AXY_W(IJK) - &
                     MU_G(IJK)*(W_G(IJK)-W_G(IJKM))*OX(I)*ODZ(K)*&
                        AXY_W(IJKM)

! Add the terms
               lTAU_W_G(IJK) =  SBV + SSX + SSY + SSZ


            ELSE
               lTAU_W_G(IJK) = ZERO
            ENDIF   ! end if (.NOT. IP_AT_T(IJK) .AND. EPGA>DIL_EP_S)
         ENDDO   ! end do ijk
!$omp end parallel do

      ELSE
! if cartesian grid
         CALL CALC_CG_TAU_W_G(lTAU_W_G)
      ENDIF

      call send_recv(ltau_w_g,2)

      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE CALC_TAU_W_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_CG_Tau_W_g                                         C
!  Purpose: Cross terms in the gradient of stress in w_g momentum      C
!  based on cartesian grid cut cell.                                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_TAU_w_G(lTAu_w_G)

! Modules
!---------------------------------------------------------------------//
      USE param
      USE param1
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE sendrecv
      USE compar
      USE fun_avg
      USE functions

      USE bc
      USE quadric
      USE cutcell
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_W_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_w_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JM, KP
      INTEGER :: IJK, IJKE, IJKW, IJKN, IJKS, IJKT
      INTEGER :: IJKNT, IJKST, IJKTE, IJKTW
      INTEGER :: IJKP, IMJK, IJMK, IJKM
      INTEGER :: IMJKP, IJMKP
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! cartesian grid
      DOUBLE PRECISION :: DEL_H, Nx, Ny, Nz
      LOGICAL :: U_NODE_AT_ET, U_NODE_AT_EB, U_NODE_AT_WT, U_NODE_AT_WB
      LOGICAL :: V_NODE_AT_NT, V_NODE_AT_NB, V_NODE_AT_ST, V_NODE_AT_SB
      DOUBLE PRECISION :: dudz_at_E, dudz_at_W
      DOUBLE PRECISION :: dvdz_at_N, dvdz_at_S
      DOUBLE PRECISION :: Xi, Yi, Zi, Ui, Vi, Sx, Sy, Sz
      DOUBLE PRECISION :: MU_G_CUT, SSX_CUT, SSY_CUT
      DOUBLE PRECISION :: UW_g, VW_g, WW_g
      INTEGER :: BCV
      INTEGER :: BCT
!---------------------------------------------------------------------//

!$omp  parallel do default(none) &
!$omp  private(I, J, K, IM, JM, KP, IJK,                               &
!$omp          IJKE, IJKW, IJKN, IJKS, IJKT,                           &
!$omp          IJKTE, IJKNT, IJKTW, IJKST,                             &
!$omp          IMJK, IJMK, IJKM, IJKP, IJMKP, IMJKP,                   &
!$omp          EPGA, SBV, SSX, SSY, SSZ,                               &
!$omp          BCV, BCT, noc_wg, uw_g, vw_g, ww_g,                     &
!$omp          del_h, nx, ny, nz, xi, yi, zi, sx, sy, sz, vi, ui,      &
!$omp          cut_tau_wg, MU_G_cut, ssy_cut, ssx_cut,                &
!$omp          dudz_at_e, dudz_at_w, dvdz_at_S, dvdz_at_N,             &
!$omp          u_node_at_wb, u_node_at_wt, u_node_at_eb, u_node_at_et, &
!$omp          v_node_at_sb, v_node_at_st, v_node_at_nb, v_node_at_nt) &
!$omp  shared(ijkstart3, ijkend3, i_of, j_of, k_of, im1, jm1, kp1,     &
!$omp         do_k, ltau_w_g,                &
!$omp         ep_g, mu_g, LAMBDA_G, trd_g, v_g, w_g, u_g,    &
!$omp         axy,  axy_w, ayz_w, axz_w, vol, ox,                      &
!$omp         bc_type, bc_w_id, bc_hw_g, bc_uw_g, bc_vw_g, bc_ww_g,    &
!$omp         oneodz_t_u, oneodz_t_v, oneodz_t_w,                      &
!$omp         x_u, y_u, z_u, x_v, y_v, z_v, x_w, y_w, z_w,             &
!$omp         wall_u_at, wall_v_at, area_w_cut, cut_w_cell_at,         &
!$omp         blocked_u_cell_at, blocked_v_cell_at)
      DO IJK = IJKSTART3, IJKEND3
         K = K_OF(IJK)
         IJKT = TOP_OF(IJK)
         EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)
         IF ( .NOT.IP_AT_T(IJK) .AND. EPGA>DIL_EP_S) THEN
            J = J_OF(IJK)
            I = I_OF(IJK)
            IM = IM1(I)
            JM = JM1(J)
            KP = KP1(K)

            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            IMJKP = KP_OF(IMJK)
            IJMKP = JM_OF(IJKP)

            IJKN = NORTH_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKE = EAST_OF(IJK)
            IJKW = WEST_OF(IJK)
            IJKNT = TOP_OF(IJKN)
            IJKST = TOP_OF(IJKS)
            IJKTE = EAST_OF(IJKT)
            IJKTW = WEST_OF(IJKT)

! bulk viscosity term
            SBV =  (LAMBDA_G(IJKT)*TRD_G(IJKT)) * AXY_W(IJK) &
                  -(LAMBDA_G(IJK) *TRD_G(IJK) ) * AXY_W(IJKM)

! shear stress terms
            IF(.NOT.CUT_W_CELL_AT(IJK))   THEN
               SSX = AVG_Z_H(AVG_X_H(MU_G(IJK),MU_G(IJKE),I),&
                             AVG_X_H(MU_G(IJKT),MU_G(IJKTE),I),K)*&
                        (U_G(IJKP)-U_G(IJK))*ONEoDZ_T_U(IJK)*AYZ_W(IJK) -&
                     AVG_Z_H(AVG_X_H(MU_G(IJKW),MU_G(IJK),IM),&
                             AVG_X_H(MU_G(IJKTW),MU_G(IJKT),IM),K)*&
                        (U_G(IMJKP)-U_G(IMJK))*ONEoDZ_T_U(IMJK)*AYZ_W(IMJK)

               SSY = AVG_Z_H(AVG_Y_H(MU_G(IJK),MU_G(IJKN),J),&
                             AVG_Y_H(MU_G(IJKT),MU_G(IJKNT),J),K)*&
                        (V_G(IJKP)-V_G(IJK))*ONEoDZ_T_V(IJK)*AXZ_W(IJK) - &
                     AVG_Z_H(AVG_Y_H(MU_G(IJKS),MU_G(IJK),JM),&
                             AVG_Y_H(MU_G(IJKST),MU_G(IJKT),JM),K)*&
                        (V_G(IJMKP)-V_G(IJMK))*ONEoDZ_T_V(IJMK)*AXZ_W(IJMK)

               SSZ = MU_G(IJKT)*(W_G(IJKP)-W_G(IJK))*&
                        ONEoDZ_T_W(IJK)*AXY_W(IJK) - &
                     MU_G(IJK)*(W_G(IJK)-W_G(IJKM))*&
                        ONEoDZ_T_W(IJKM)*AXY_W(IJKM)

! cut cell modifications
!---------------------------------------------------------------------//
            ELSE
               BCV = BC_W_ID(IJK)

               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               SELECT CASE (BCT)
                  CASE (CG_NSW,CG_MI)
                     CUT_TAU_WG = .TRUE.
                     NOC_WG     = .TRUE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE (CG_FSW)
                     CUT_TAU_WG = .FALSE.
                     NOC_WG     = .FALSE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE(CG_PSW)
                     IF(BC_HW_G(BC_W_ID(IJK))==UNDEFINED) THEN   ! same as NSW
                        CUT_TAU_WG = .TRUE.
                        NOC_WG     = .TRUE.
                        UW_g = BC_UW_G(BCV)
                        VW_g = BC_VW_G(BCV)
                        WW_g = BC_WW_G(BCV)
                     ELSEIF(BC_HW_G(BC_W_ID(IJK))==ZERO) THEN   ! same as FSW
                        CUT_TAU_WG = .FALSE.
                        NOC_WG     = .FALSE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     ELSE                              ! partial slip
                        CUT_TAU_WG = .FALSE.
                        NOC_WG     = .FALSE.
                     ENDIF
                  CASE (NONE)
                     lTAU_W_G(IJK) = ZERO
                     CYCLE
               END SELECT

               IF(CUT_TAU_WG) THEN
                  MU_G_CUT = (VOL(IJK)*MU_G(IJK) + &
                     VOL(IJKP)*MU_G(IJKT))/(VOL(IJK) + VOL(IJKP))
               ELSE
                  MU_G_CUT = ZERO
               ENDIF

! SSX:
               U_NODE_AT_ET = ((.NOT.BLOCKED_U_CELL_AT(IJKP)).AND.&
                               (.NOT.WALL_U_AT(IJKP)))
               U_NODE_AT_EB = ((.NOT.BLOCKED_U_CELL_AT(IJK)).AND.&
                               (.NOT.WALL_U_AT(IJK)))
               U_NODE_AT_WT = ((.NOT.BLOCKED_U_CELL_AT(IMJKP)).AND.&
                               (.NOT.WALL_U_AT(IMJKP)))
               U_NODE_AT_WB = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.&
                               (.NOT.WALL_U_AT(IMJK)))

               IF(U_NODE_AT_ET.AND.U_NODE_AT_EB) THEN
                  Ui = HALF * (U_G(IJKP) + U_G(IJK))
                  Xi = HALF * (X_U(IJKP) + X_U(IJK))
                  Yi = HALF * (Y_U(IJKP) + Y_U(IJK))
                  Zi = HALF * (Z_U(IJKP) + Z_U(IJK))
                  Sx = X_U(IJKP) - X_U(IJK)
                  Sy = Y_U(IJKP) - Y_U(IJK)
                  Sz = Z_U(IJKP) - Z_U(IJK)
                  CALL GET_DEL_H(IJK, 'W_MOMENTUM', Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dudz_at_E =  (U_G(IJKP) - U_G(IJK)) * &
                     ONEoDZ_T_U(IJK)
                  IF(NOC_WG) dudz_at_E = dudz_at_E - ((Ui-UW_g) * &
                     ONEoDZ_T_U(IJK)/DEL_H*(Sx*Nx+Sy*Ny))
               ELSE
                  dudz_at_E =  ZERO
               ENDIF

               IF(U_NODE_AT_WT.AND.U_NODE_AT_WB) THEN
                  Ui = HALF * (U_G(IMJKP) + U_G(IMJK))
                  Xi = HALF * (X_U(IMJKP) + X_U(IMJK))
                  Yi = HALF * (Y_U(IMJKP) + Y_U(IMJK))
                  Zi = HALF * (Z_U(IMJKP) + Z_U(IMJK))
                  Sx = X_U(IMJKP) - X_U(IMJK)
                  Sy = Y_U(IMJKP) - Y_U(IMJK)
                  Sz = Z_U(IMJKP) - Z_U(IMJK)
                  CALL GET_DEL_H(IJK, 'W_MOMENTUM', Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dudz_at_W =  (U_G(IMJKP) - U_G(IMJK)) * &
                     ONEoDZ_T_U(IMJK)
                  IF(NOC_WG) dudz_at_W = dudz_at_W - ((Ui-UW_g) * &
                     ONEoDZ_T_U(IMJK)/DEL_H*(Sx*Nx+Sy*Ny))
               ELSE
                  dudz_at_W =  ZERO
               ENDIF

               IF(U_NODE_AT_EB) THEN
                  CALL GET_DEL_H(IJK, 'W_MOMENTUM', X_U(IJK), &
                     Y_U(IJK), Z_U(IJK), Del_H, Nx, Ny, Nz)
                  SSX_CUT = - MU_G_CUT * (U_G(IJK) - UW_g) / &
                     DEL_H * (Nz*Nx) * Area_W_CUT(IJK)
               ELSE
                  SSX_CUT =  ZERO
               ENDIF
               SSX = AVG_Z_H(AVG_X_H(MU_G(IJK),MU_G(IJKE),I),&
                             AVG_X_H(MU_G(IJKT),MU_G(IJKTE),I),K)*&
                        dudz_at_E*AYZ_W(IJK) - &
                     AVG_Z_H(AVG_X_H(MU_G(IJKW),MU_G(IJK),IM),&
                             AVG_X_H(MU_G(IJKTW),MU_G(IJKT),IM),K)*&
                        dudz_at_W*AYZ_W(IMJK) + SSX_CUT

! SSY:
               V_NODE_AT_NT = ((.NOT.BLOCKED_V_CELL_AT(IJKP)).AND.&
                               (.NOT.WALL_V_AT(IJKP)))
               V_NODE_AT_NB = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.&
                               (.NOT.WALL_V_AT(IJK)))
               V_NODE_AT_ST = ((.NOT.BLOCKED_V_CELL_AT(IJMKP)).AND.&
                               (.NOT.WALL_V_AT(IJMKP)))
               V_NODE_AT_SB = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.&
                               (.NOT.WALL_V_AT(IJMK)))

               IF(V_NODE_AT_NT.AND.V_NODE_AT_NB) THEN
                  Vi = HALF * (V_G(IJKP) + V_G(IJK))
                  Xi = HALF * (X_V(IJKP) + X_V(IJK))
                  Yi = HALF * (Y_V(IJKP) + Y_V(IJK))
                  Zi = HALF * (Z_V(IJKP) + Z_V(IJK))
                  Sx = X_V(IJKP) - X_V(IJK)
                  Sy = Y_V(IJKP) - Y_V(IJK)
                  Sz = Z_V(IJKP) - Z_V(IJK)
                  CALL GET_DEL_H(IJK, 'W_MOMENTUM', Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dvdz_at_N =  (V_G(IJKP) - V_G(IJK)) * &
                     ONEoDZ_T_V(IJK)
                  IF(NOC_WG) dvdz_at_N = dvdz_at_N - ((Vi-VW_g) * &
                     ONEoDZ_T_V(IJK)/DEL_H*(Sx*Nx+Sy*Ny))
               ELSE
                     dvdz_at_N =  ZERO
               ENDIF

               IF(V_NODE_AT_ST.AND.V_NODE_AT_SB) THEN
                  Vi = HALF * (V_G(IJMKP) + V_G(IJMK))
                  Xi = HALF * (X_V(IJMKP) + X_V(IJMK))
                  Yi = HALF * (Y_V(IJMKP) + Y_V(IJMK))
                  Zi = HALF * (Z_V(IJMKP) + Z_V(IJMK))
                  Sx = X_V(IJMKP) - X_V(IJMK)
                  Sy = Y_V(IJMKP) - Y_V(IJMK)
                  Sz = Z_V(IJMKP) - Z_V(IJMK)
                  CALL GET_DEL_H(IJK, 'W_MOMENTUM', Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dvdz_at_S =  (V_G(IJMKP) - V_G(IJMK)) * &
                     ONEoDZ_T_V(IJMK)
                  IF(NOC_WG) dvdz_at_S = dvdz_at_S - ((Vi-VW_g) * &
                     ONEoDZ_T_V(IJMK)/DEL_H*(Sx*Nx+Sy*Ny))
               ELSE
                  dvdz_at_S =  ZERO
               ENDIF

               IF(V_NODE_AT_NB) THEN
                  CALL GET_DEL_H(IJK, 'W_MOMENTUM', X_V(IJK), &
                     Y_V(IJK), Z_V(IJK), Del_H, Nx, Ny, Nz)
                  SSY_CUT = - MU_G_CUT * (V_G(IJK) - VW_g) / &
                     DEL_H * (Nz*Ny) * Area_W_CUT(IJK)
               ELSE
                  SSY_CUT =  ZERO
               ENDIF

               SSY = AVG_Z_H(AVG_Y_H(MU_G(IJK),MU_G(IJKN),J),&
                             AVG_Y_H(MU_G(IJKT),MU_G(IJKNT),J),K)*&
                        dvdz_at_N*AXZ_W(IJK) - &
                     AVG_Z_H(AVG_Y_H(MU_G(IJKS),MU_G(IJK),JM),&
                             AVG_Y_H(MU_G(IJKST),MU_G(IJKT),JM),K)*&
                        dvdz_at_S*AXZ_W(IJMK) + SSY_CUT

! SSZ:
               CALL GET_DEL_H(IJK, 'W_MOMENTUM', X_W(IJK), &
                  Y_W(IJK), Z_W(IJK), Del_H, Nx, Ny, Nz)

               SSZ = MU_G(IJKT)*(W_G(IJKP)-W_G(IJK))*&
                        ONEoDZ_T_W(IJK)*AXY_W(IJK) - &
                     MU_G(IJK)*(W_G(IJK)-W_G(IJKM))*&
                        OX(I)*ONEoDZ_T_W(IJKM)*AXY_W(IJKM) - &
                     MU_G_CUT * (W_g(IJK) - WW_g) / DEL_H * &
                        (Nz**2) * Area_W_CUT(IJK)

            ENDIF  ! end if/else cut_w_cell_at

! Add the terms
            lTAU_W_G(IJK) =  SBV + SSX + SSY + SSZ

         ELSE
            lTAU_W_G(IJK) = ZERO
         ENDIF   ! end if (.NOT. IP_AT_T(IJK) .AND. EPGA>DIL_EP_S)
      ENDDO   ! end do ijk
!$omp end parallel do

      RETURN
      END SUBROUTINE CALC_CG_TAU_W_G
