!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_V_g                                            C
!  Purpose: Cross terms in the gradient of stress in V_g momentum      c
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
!                                                                      C
!                                                                      C
!  Comments: This routine calculates the components of the gas phase   C
!  viscous stress tensor of the v-momentum equation that cannot be     C
!  cast in the form: mu.grad(v). These components are stored in the    C
!  passed variable, which is then used as a source of the v-momentum   C
!  equation.                                                           C
!                                                                      C
!  The following v component viscous stress tensor terms are           C
!  calculated here:                                                    C
!  > part of d/dy (tau_yy) xdxdydz =>                                  C
!            d/dy (lambda.trcD) xdxdydz =>                             C
!    delta (lambda.trcD)Ap|N-S                                         C
!  > part of 1/x d/dx(x.tau_xy) xdxdydz =>                             C
!            1/x d/dx (x.mu.du/dy) xdxdydz =>                          C
!    delta (x.mu.du/dy)Ayz |E-W                                        C
!  > part of d/dy (tau_xy) xdxdydz =>                                  C
!           d/dy (mu.dv/dy) xdxdydz =>                                 C
!    delta (mu.dv/dx)Axz |N-S                                          C
!  > part of 1/x d/dz (tau_xz) xdxdydz =>                              C
!            1/x d/dz (mu.dw/dy) xdxdydz =>                            C
!    delta (mu.dw/dx)Axy |T-B                                          C
!                                                                      C
!  To reconstitute the full v-momentum gas phase viscous stress        C
!  tensor would require including the the 'diffusional' components     C
!  (i.e., those of the form mu.grad(v)                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_V_G(lTAU_V_G)

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
! TAU_V_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_V_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JP, KM
      INTEGER :: IJK, IJKN, IJKE, IJKW, IJKT, IJKB
      INTEGER :: IJKTN, IJKBN, IJKNE, IJKNW
      INTEGER :: IJPK, IJMK, IMJK, IJKM
      INTEGER :: IMJPK, IJPKM
      INTEGER :: itmp, jtmp, ktmp
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
!---------------------------------------------------------------------//
!     NOTE -- triply nested functions seem to break things -- hence the
!             use of the *tmp variables below
!---------------------------------------------------------------------//

      IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(4)==1)) THEN

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
            IJKN = NORTH_OF(IJK)
            EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)
            IF ( .NOT.IP_AT_N(IJK) .AND. EPGA>DIL_EP_S) THEN

               JP = JP1(J)
               IM = IM1(I)
               KM = KM1(K)

               IJMK = FUNIJK(i,jminus(i,j,k),k)
               IMJK = FUNIJK(iminus(i,j,k),j,k)

               ktmp = kminus(i,j,k)
               IJKM  = FUNIJK(i,j,ktmp)
               IJPKM = FUNIJK(i,jplus(i,j,ktmp),ktmp)

               jtmp = jplus(i,j,k)
               IJPK  = FUNIJK(i,jtmp,k)
               IMJPK = FUNIJK(iminus(i,jtmp,k),jtmp,k)

               IJKE = FUNIJK(ieast(i,j,k),j,k)
               IJKW = FUNIJK(iwest(i,j,k),j,k)
               IJKT = FUNIJK(i,j,ktop(i,j,k))
               IJKB = FUNIJK(i,j,kbot(i,j,k))

               jtmp = jnorth(i,j,k)
               IJKNE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

               itmp = iwest(i,j,k)
               IJKNW = FUNIJK(itmp,jnorth(itmp,j,k),k)

               ktmp = ktop(i,j,k)
               IJKTN = FUNIJK(i,jnorth(i,j,ktmp),ktmp)

               ktmp = kbot(i,j,k)
               IJKBN = FUNIJK(i,jnorth(i,j,ktmp),ktmp)

! Surface forces at i, j+1/2, k
! bulk viscosity term
! part of d/dy (tau_yy) xdxdydz =>
!         d/dy (lambda.trcD) xdxdydz =>
! delta (lambda.trcD)Ap|N-S  : at (i, j+1 - j-1, k)
               SBV = (LAMBDA_G(IJKN)*TRD_G(IJKN)-&
                      LAMBDA_G(IJK)*TRD_G(IJK))*AXZ(IJK)

! shear stress terms at i, j+1/2, k

! part of 1/x d/dx(x.tau_xy) xdxdydz =>
!         1/x d/dx (x.mu.du/dy) xdxdydz =>
! delta (x.mu.du/dy)Ayz |E-W : at (i+1/2 - i-1/2, j+1/2, k)
               SSX = AVG_Y_H(AVG_X_H(MU_G(IJK),MU_G(IJKE),I),&
                             AVG_X_H(MU_G(IJKN),MU_G(IJKNE),I),J)*&
                     (U_G(IJPK)-U_G(IJK))*ODY_N(J)*AYZ_V(IJK) - &
                     AVG_Y_H(AVG_X_H(MU_G(IJKW),MU_G(IJK),IM),&
                             AVG_X_H(MU_G(IJKNW),MU_G(IJKN),IM),J)*&
                     (U_G(IMJPK)-U_G(IMJK))*ODY_N(J)*AYZ_V(IMJK)

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dy) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at (i, j+1 - j-1, k)
               SSY = MU_G(IJKN)*(V_G(IJPK)-V_G(IJK))*ODY(JP)*&
                        AXZ_V(IJK) - &
                     MU_G(IJK)*(V_G(IJK)-V_G(IJMK))*ODY(J)*&
                        AXZ_V(IJMK)

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dy) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at (i, j+1/2, k+1/2 - k-1/2)
               SSZ = AVG_Y_H(AVG_Z_H(MU_G(IJK),MU_G(IJKT),K),&
                             AVG_Z_H(MU_G(IJKN),MU_G(IJKTN),K),J)*&
                     (W_G(IJPK)-W_G(IJK))*ODY_N(J)*AXY_V(IJK) - &
                     AVG_Y_H(AVG_Z_H(MU_G(IJKB),MU_G(IJK),KM),&
                             AVG_Z_H(MU_G(IJKBN),MU_G(IJKN),KM),J)*&
                     (W_G(IJPKM)-W_G(IJKM))*ODY_N(J)*AXY_V(IJKM)

! Add the terms
               lTAU_V_G(IJK) = SBV + SSX + SSY + SSZ

            ELSE
               lTAU_V_G(IJK) = ZERO
            ENDIF   ! end if (.NOT. IP_AT_N(IJK) .AND. EPGA>DIL_EP_S)
         ENDDO
         ENDDO
         ENDDO

      ELSE
! cartesian grid case
         CALL CALC_CG_TAU_V_G(lTAU_V_G)
      ENDIF

      call send_recv(ltau_v_g,2)
      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE CALC_TAU_V_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_CG_Tau_v_g                                         C
!  Purpose: Cross terms in the gradient of stress in V_g momentum      C
!  based on cartesian grid cut cell.                                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_TAU_V_G(lTAU_V_G)

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

      USE cutcell

      USE bc
      USE quadric
      USE cutcell

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_V_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_V_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JP, KM
      INTEGER :: IJK, IJKN, IJKE, IJKW, IJKT, IJKB
      INTEGER :: IJKTN, IJKBN, IJKNE, IJKNW
      INTEGER :: IJPK, IJMK, IMJK, IJKM
      INTEGER :: IMJPK, IJPKM
! Average volume fraction
      DOUBLE PRECISION EPGA
! Source terms (Surface)
      DOUBLE PRECISION Sbv, Ssx, Ssy, Ssz

! Cartesian grid variables
      DOUBLE PRECISION :: DEL_H, Nx, Ny, Nz
      LOGICAL :: U_NODE_AT_NE, U_NODE_AT_NW, U_NODE_AT_SE, U_NODE_AT_SW
      LOGICAL :: W_NODE_AT_TN, W_NODE_AT_TS, W_NODE_AT_BN, W_NODE_AT_BS
      DOUBLE PRECISION :: dudy_at_E, dudy_at_W
      DOUBLE PRECISION :: dwdy_at_T, dwdy_at_B
      DOUBLE PRECISION :: Xi, Yi, Zi, Ui, Wi, Sx, Sy, Sz
      DOUBLE PRECISION :: MU_G_CUT, SSX_CUT, SSZ_CUT
      DOUBLE PRECISION :: UW_g, VW_g, WW_g
      INTEGER :: BCV
      INTEGER :: BCT
!---------------------------------------------------------------------//

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
         IJKN = NORTH_OF(IJK)
         EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)
         IF ( .NOT.IP_AT_N(IJK) .AND. EPGA>DIL_EP_S) THEN
            JP = JP1(J)
            IM = IM1(I)
            KM = KM1(K)

            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IMJK = IM_OF(IJK)
            IJKM = KM_OF(IJK)
            IMJPK = IM_OF(IJPK)
            IJPKM = JP_OF(IJKM)

            IJKE = EAST_OF(IJK)
            IJKNE = EAST_OF(IJKN)
            IJKW = WEST_OF(IJK)
            IJKNW = NORTH_OF(IJKW)
            IJKT = TOP_OF(IJK)
            IJKTN = NORTH_OF(IJKT)
            IJKB = BOTTOM_OF(IJK)
            IJKBN = NORTH_OF(IJKB)


! bulk viscosity term
            SBV =  (LAMBDA_G(IJKN)*TRD_G(IJKN)) * AXZ_V(IJK) - &
                   (LAMBDA_G(IJK) *TRD_G(IJK) ) * AXZ_V(IJMK)

! shear stress terms
            IF(.NOT.CUT_V_CELL_AT(IJK)) THEN
               SSX = AVG_Y_H(AVG_X_H(MU_G(IJK),MU_G(IJKE),I),&
                             AVG_X_H(MU_G(IJKN),MU_G(IJKNE),I),J)*&
                     (U_G(IJPK)-U_G(IJK))*ONEoDY_N_U(IJK)*AYZ_V(IJK) - &
                     AVG_Y_H(AVG_X_H(MU_G(IJKW),MU_G(IJK),IM),&
                             AVG_X_H(MU_G(IJKNW),MU_G(IJKN),IM),J)*&
                     (U_G(IMJPK)-U_G(IMJK))*ONEoDY_N_U(IMJK)*AYZ_V(IMJK)

               SSY = MU_G(IJKN)*(V_G(IJPK)-V_G(IJK))*&
                        ONEoDY_N_V(IJK)*AXZ_V(IJK) - &
                     MU_G(IJK)*(V_G(IJK)-V_G(IJMK))*&
                        ONEoDY_N_V(IJMK)*AXZ_V(IJMK)

               IF(DO_K) THEN
                  SSZ = AVG_Y_H(AVG_Z_H(MU_G(IJK),MU_G(IJKT),K),&
                                AVG_Z_H(MU_G(IJKN),MU_G(IJKTN),K),J)*&
                        (W_G(IJPK)-W_G(IJK))*ONEoDY_N_W(IJK)*AXY_V(IJK) - &
                        AVG_Y_H(AVG_Z_H(MU_G(IJKB),MU_G(IJK),KM),&
                                AVG_Z_H(MU_G(IJKBN),MU_G(IJKN),KM),J)*&
                        (W_G(IJPKM)-W_G(IJKM))*ONEoDY_N_W(IJKM)*AXY_V(IJKM)
               ELSE
                  SSZ = ZERO
               ENDIF

! cut cell modifications
!---------------------------------------------------------------------//
            ELSE

               BCV = BC_V_ID(IJK)
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               SELECT CASE (BCT)
                  CASE (CG_NSW,CG_MI)
                     CUT_TAU_VG = .TRUE.
                     NOC_VG     = .TRUE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE (CG_FSW)
                     CUT_TAU_VG = .FALSE.
                     NOC_VG     = .FALSE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE(CG_PSW)
                     IF(BC_HW_G(BC_V_ID(IJK))==UNDEFINED) THEN   ! same as NSW
                        CUT_TAU_VG = .TRUE.
                        NOC_VG     = .TRUE.
                        UW_g = BC_UW_G(BCV)
                        VW_g = BC_VW_G(BCV)
                        WW_g = BC_WW_G(BCV)
                     ELSEIF(BC_HW_G(BC_V_ID(IJK))==ZERO) THEN   ! same as FSW
                        CUT_TAU_VG = .FALSE.
                        NOC_VG     = .FALSE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     ELSE                              ! partial slip
                        CUT_TAU_VG = .FALSE.
                        NOC_VG     = .FALSE.
                     ENDIF
                  CASE (NONE)
                     lTAU_V_G(IJK) = ZERO
                     CYCLE
               END SELECT


               IF(CUT_TAU_VG) THEN
                  MU_G_CUT = (VOL(IJK)*MU_G(IJK) + &
                               VOL(IJPK)*MU_G(IJKN))/&
                              (VOL(IJK) + VOL(IJPK))
               ELSE
                  MU_G_CUT = ZERO
               ENDIF

! SSX:
               U_NODE_AT_NE = ((.NOT.BLOCKED_U_CELL_AT(IJPK)).AND.&
                               (.NOT.WALL_U_AT(IJPK)))
               U_NODE_AT_SE = ((.NOT.BLOCKED_U_CELL_AT(IJK)).AND.&
                               (.NOT.WALL_U_AT(IJK)))
               U_NODE_AT_NW = ((.NOT.BLOCKED_U_CELL_AT(IMJPK)).AND.&
                               (.NOT.WALL_U_AT(IMJPK)))
               U_NODE_AT_SW = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.&
                               (.NOT.WALL_U_AT(IMJK)))

               IF(U_NODE_AT_NE.AND.U_NODE_AT_SE) THEN
                  Ui = HALF * (U_G(IJPK) + U_G(IJK))
                  Xi = HALF * (X_U(IJPK) + X_U(IJK))
                  Yi = HALF * (Y_U(IJPK) + Y_U(IJK))
                  Zi = HALF * (Z_U(IJPK) + Z_U(IJK))
                  Sx = X_U(IJPK) - X_U(IJK)
                  Sy = Y_U(IJPK) - Y_U(IJK)
                  Sz = Z_U(IJPK) - Z_U(IJK)
                  CALL GET_DEL_H(IJK, 'V_MOMENTUM', Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dudy_at_E =  (U_G(IJPK) - U_G(IJK)) * ONEoDY_N_U(IJK)
                  IF(NOC_VG) dudy_at_E = dudy_at_E - ((Ui-UW_g) * &
                     ONEoDY_N_U(IJK)/DEL_H * (Sx*Nx+Sz*Nz) )
               ELSE
                  dudy_at_E =  ZERO
               ENDIF

               IF(U_NODE_AT_NW.AND.U_NODE_AT_SW) THEN
                  Ui = HALF * (U_G(IMJPK) + U_G(IMJK))
                  Xi = HALF * (X_U(IMJPK) + X_U(IMJK))
                  Yi = HALF * (Y_U(IMJPK) + Y_U(IMJK))
                  Zi = HALF * (Z_U(IMJPK) + Z_U(IMJK))
                  Sx = X_U(IMJPK) - X_U(IMJK)
                  Sy = Y_U(IMJPK) - Y_U(IMJK)
                  Sz = Z_U(IMJPK) - Z_U(IMJK)
                  CALL GET_DEL_H(IJK, 'V_MOMENTUM', Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dudy_at_W =  (U_G(IMJPK) - U_G(IMJK)) * ONEoDY_N_U(IMJK)
                  IF(NOC_VG) dudy_at_W = dudy_at_W - ((Ui-UW_g) * &
                     ONEoDY_N_U(IMJK)/DEL_H * (Sx*Nx+Sz*Nz) )
               ELSE
                  dudy_at_W =  ZERO
               ENDIF

               IF(U_NODE_AT_SE) THEN
                  CALL GET_DEL_H(IJK,'V_MOMENTUM', X_U(IJK), Y_U(IJK), &
                     Z_U(IJK), Del_H, Nx, Ny, Nz)
                  SSX_CUT = - MU_G_CUT * (U_G(IJK) - UW_g) / DEL_H * &
                     (Ny*Nx) * Area_V_CUT(IJK)
               ELSE
                  SSX_CUT =  ZERO
               ENDIF

               SSX = AVG_Y_H(AVG_X_H(MU_G(IJK),MU_G(IJKE),I),&
                             AVG_X_H(MU_G(IJKN),MU_G(IJKNE),I),J)*&
                     dudy_at_E*AYZ_V(IJK) - &
                     AVG_Y_H(AVG_X_H(MU_G(IJKW),MU_G(IJK),IM),&
                             AVG_X_H(MU_G(IJKNW),MU_G(IJKN),IM),J)*&
                     dudy_at_W*AYZ_V(IMJK) + SSX_CUT

! SSY:
               CALL GET_DEL_H(IJK, 'V_MOMENTUM', X_V(IJK), Y_V(IJK), &
                  Z_V(IJK), Del_H, Nx, Ny, Nz)

               SSY = MU_G(IJKN)*(V_G(IJPK)-V_G(IJK))*&
                        ONEoDY_N_V(IJK)*AXZ_V(IJK) - &
                     MU_G(IJK)*(V_G(IJK)-V_G(IJMK))*&
                        ONEoDY_N_V(IJMK)*AXZ_V(IJMK) - &
                     MU_G_CUT * (V_g(IJK) - VW_g) / DEL_H * &
                        (Ny**2) * Area_V_CUT(IJK)

! SSZ:
               IF(DO_K) THEN
                  W_NODE_AT_TN = ((.NOT.BLOCKED_W_CELL_AT(IJPK)).AND.&
                                  (.NOT.WALL_W_AT(IJPK)))
                  W_NODE_AT_TS = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.&
                                  (.NOT.WALL_W_AT(IJK)))
                  W_NODE_AT_BN = ((.NOT.BLOCKED_W_CELL_AT(IJPKM)).AND.&
                                  (.NOT.WALL_W_AT(IJPKM)))
                  W_NODE_AT_BS = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.&
                                  (.NOT.WALL_W_AT(IJKM)))

                  IF(W_NODE_AT_TN.AND.W_NODE_AT_TS) THEN
                        Wi = HALF * (W_G(IJPK) + W_G(IJK))
                        Xi = HALF * (X_W(IJPK) + X_W(IJK))
                        Yi = HALF * (Y_W(IJPK) + Y_W(IJK))
                        Zi = HALF * (Z_W(IJPK) + Z_W(IJK))
                        Sx = X_W(IJPK) - X_W(IJK)
                        Sy = Y_W(IJPK) - Y_W(IJK)
                        Sz = Z_W(IJPK) - Z_W(IJK)
                        CALL GET_DEL_H(IJK, 'V_MOMENTUM', Xi, Yi, Zi, &
                           Del_H, Nx, Ny, Nz)

                        dwdy_at_T = (W_G(IJPK)-W_G(IJK)) * &
                           ONEoDY_N_W(IJK)
                        IF(NOC_VG) dwdy_at_T = dwdy_at_T - ((Wi-WW_g) * &
                           ONEoDY_N_W(IJK)/DEL_H*(Sx*Nx+Sz*Nz))
                  ELSE
                     dwdy_at_T =  ZERO
                  ENDIF

                  IF(W_NODE_AT_BN.AND.W_NODE_AT_BS) THEN
                     Wi = HALF * (W_G(IJPKM) + W_G(IJKM))
                     Xi = HALF * (X_W(IJPKM) + X_W(IJKM))
                     Yi = HALF * (Y_W(IJPKM) + Y_W(IJKM))
                     Zi = HALF * (Z_W(IJPKM) + Z_W(IJKM))
                     Sx = X_W(IJPKM) - X_W(IJKM)
                     Sy = Y_W(IJPKM) - Y_W(IJKM)
                     Sz = Z_W(IJPKM) - Z_W(IJKM)
                     CALL GET_DEL_H(IJK, 'V_MOMENTUM', Xi, Yi, Zi,&
                        Del_H, Nx, Ny, Nz)

                     dwdy_at_B =  (W_G(IJPKM) - W_G(IJKM)) * &
                        ONEoDY_N_W(IJKM)
                     IF(NOC_VG) dwdy_at_B = dwdy_at_B - ((Wi-WW_g) * &
                        ONEoDY_N_W(IJKM)/DEL_H*(Sx*Nx+Sz*Nz))
                  ELSE
                     dwdy_at_B =  ZERO
                  ENDIF

                  IF(W_NODE_AT_TS) THEN
                     CALL GET_DEL_H(IJK, 'V_MOMENTUM', X_W(IJK), &
                        Y_W(IJK), Z_W(IJK), Del_H, Nx, Ny, Nz)
                     SSZ_CUT = - MU_G_CUT * (W_G(IJK) - WW_g) / &
                        DEL_H * (Ny*Nz) * Area_V_CUT(IJK)
                  ELSE
                     SSZ_CUT =  ZERO
                  ENDIF

                  SSZ = AVG_Y_H(AVG_Z_H(MU_G(IJK),MU_G(IJKT),K),&
                                AVG_Z_H(MU_G(IJKN),MU_G(IJKTN),K),J)*&
                           dwdy_at_T*AXY_V(IJK) - &
                        AVG_Y_H(AVG_Z_H(MU_G(IJKB),MU_G(IJK),KM),&
                                AVG_Z_H(MU_G(IJKBN),MU_G(IJKN),KM),J)*&
                           dwdy_at_B*AXY_V(IJKM) + SSZ_CUT
               ELSE
                  SSZ = ZERO
               ENDIF  ! DO_K

            ENDIF   ! end if/else cut_v_cell_at

! Add the terms
            lTAU_V_G(IJK) = SBV + SSX + SSY + SSZ

         ELSE
            lTAU_V_G(IJK) = ZERO
         ENDIF   ! end if/else (.NOT. IP_AT_N(IJK) .AND. EPGA>DIL_EP_S)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_CG_TAU_V_G
