!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_U_g                                            C
!  Purpose: Cross terms in the gradient of stress in U_g momentum      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
!                                                                      C
!                                                                      C
!  Comments: This routine calculates the components of the gas phase   C
!  viscous stress tensor of the u-momentum equation that cannot be     C
!  cast in the form: mu.grad(u). These components are stored in the    C
!  passed variable, which is then used as a source of the u-momentum   C
!  equation.                                                           C
!                                                                      C
!  The following u component viscous stress tensor terms are           C
!  calculated here:                                                    C
!  > Combined part of 1/x d/dx (x.tau_xx) xdxdydz and                  C
!                              -tau_zz/x xdxdydz =>                    C
!    1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>     C
!        d/dx (lambda.trcD) xdxdydz =>                                 C
!    delta (lambda.trcD)Ap|E-W                                         C
!  > part of 1/x d/dx(x.tau_xx) xdxdydz =>                             C
!            1/x d/dx (x.mu.du/dx) xdxdydz =>                          C
!    delta (mu du/dx)Ayz |E-W                                          C
!  > part of d/dy (tau_xy) xdxdydz =>                                  C
!           d/dy (mu.dv/dx) xdxdydz =>                                 C
!    delta (mu.dv/dx)Axz |N-S                                          C
!  > part of 1/x d/dz (tau_xz) xdxdydz =>                              C
!            1/x d/dz (mu.dw/dx) xdxdydz =>                            C
!    delta (mu.dw/dx)Axy |T-B                                          C
!                                                                      C
!  To reconstitute the complete u-momentum gas phase viscous stress    C
!  tensor would require including the term calculated in source_u_g    C
!  and the 'diffusional' components (i.e., those of the form           C
!  mu.grad(u)                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_U_G(lTAU_U_G)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3
      USE param1, only: zero, half
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc, only: dil_ep_s
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE functions
      USE cutcell
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_U_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_U_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IP, JM, KM
      INTEGER :: IJK, IJKE, IJKN, IJKS, IJKT, IJKB
      INTEGER :: IJKNE, IJKSE, IJKTE, IJKBE
      INTEGER :: IPJK, IMJK, IJMK, IJKM
      INTEGER :: IPJMK, IPJKM
      INTEGER :: itmp, jtmp, ktmp
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average viscosity
      DOUBLE PRECISION :: MU_Ge, MU_gbe, MUGA
! Average dW/Xdz
      DOUBLE PRECISION :: dWoXdz
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! Source terms (Volumetric)
      DOUBLE PRECISION :: Vtzb
!---------------------------------------------------------------------//
!     NOTE -- triply nested functions seem to break things -- hence the
!             use of the *tmp variables below
!---------------------------------------------------------------------//


      IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(3)==1)) THEN

        DO K = kstart3, kend3
          DO J = jstart3, jend3
            DO I = istart3, iend3

            IJK  = FUNIJK(i,j,k)
            IJKE = FUNIJK(ieast(i,j,k),j,k)

            EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

            IF (.NOT.IP_AT_E(IJK) .AND. EPGA>DIL_EP_S) THEN

               IP = IP1(I)
               JM = JM1(J)
               KM = KM1(K)

               IPJK = FUNIJK(iplus(i,j,k),j,k)
               IMJK = FUNIJK(iminus(i,j,k),j,k)
               IJMK = FUNIJK(i,jminus(i,j,k),k)
               IJKM = FUNIJK(i,j,kminus(i,j,k))

               itmp = iplus(i,j,k)
               IPJMK = FUNIJK(itmp,jminus(itmp,j,k),k)

               ktmp = kminus(i,j,k)
               IPJKM = FUNIJK(iplus(i,j,ktmp),j,ktmp)

               IJKE = FUNIJK(ieast(i,j,k),j,k)

               jtmp = jnorth(i,j,k)
               IJKN  = FUNIJK(i,jtmp,k)
               IJKNE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

               jtmp = jsouth(i,j,k)
               IJKS  = FUNIJK(i,jtmp,k)
               IJKSE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

               ktmp = ktop(i,j,k)
               IJKT  = FUNIJK(i,j,ktmp)
               IJKTE = FUNIJK(ieast(i,j,ktmp),j,ktmp)

               ktmp = kbot(i,j,k)
               IJKB  = FUNIJK(i,j,ktmp)
               IJKBE = FUNIJK(ieast(i,j,ktmp),j,ktmp)


! Surface forces at i+1/2, j, k
! bulk viscosity term
! combines part of 1/x d/dx (x.tau_xx) xdxdydz and -tau_zz/x xdxdydz =>
! combines 1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>
!              d/dx (lambda.trcD) xdxdydz
! delta (lambda.trcD)Ap |E-W : at (i+1 - i-1), j, k
               SBV = (LAMBDA_G(IJKE)*TRD_G(IJKE)-&
                      LAMBDA_G(IJK)*TRD_G(IJK))*AYZ(IJK)

! shear stress terms at i+1/2, j, k
! part of 1/x d/dx(x.tau_xx) xdxdydz =>
!         1/x d/dx (x.mu.du/dx) xdxdydz =>
! delta (mu du/dx)Ayz |E-W : at (i+1 - i-1), j, k
               SSX = MU_G(IJKE)*(U_G(IPJK)-U_G(IJK))*ODX(IP)*&
                        AYZ_U(IJK) - &
                     MU_G(IJK)*(U_G(IJK)-U_G(IMJK))*ODX(I)*&
                        AYZ_U(IMJK)

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dx) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at i+1/2, (j+1/2 - j-1/2), k
               SSY = AVG_X_H(AVG_Y_H(MU_G(IJK),MU_G(IJKN),J),&
                             AVG_Y_H(MU_G(IJKE),MU_G(IJKNE),J),I)*&
                        (V_G(IPJK)-V_G(IJK))*ODX_E(I)*AXZ_U(IJK) - &
                     AVG_X_H(AVG_Y_H(MU_G(IJKS),MU_G(IJK),JM),&
                             AVG_Y_H(MU_G(IJKSE),MU_G(IJKE),JM),I)*&
                        (V_G(IPJMK)-V_G(IJMK))*ODX_E(I)*AXZ_U(IJMK)

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dx) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at i+1/2, j, (k+1/2 - k-1/2)
               MU_GE = AVG_X_H(AVG_Z_H(MU_G(IJK),MU_G(IJKT),K),&
                                AVG_Z_H(MU_G(IJKE),MU_G(IJKTE),K),I)
               MU_GBE = AVG_X_H(AVG_Z_H(MU_G(IJKB),MU_G(IJK),KM),&
                                AVG_Z_H(MU_G(IJKBE),MU_G(IJKE),KM),I)
               SSZ = MU_GE*(W_G(IPJK)-W_G(IJK))*ODX_E(I)*&
                        AXY_U(IJK) - &
                     MU_GBE*(W_G(IPJKM)-W_G(IJKM))*ODX_E(I)*&
                        AXY_U(IJKM)


! Add the terms
               lTAU_U_G(IJK) = SBV + SSX + SSY + SSZ

            ELSE
               lTAU_U_G(IJK) = ZERO
            ENDIF   ! end if (.NOT.IP_AT_E(IJK) .AND. EPGA>DIL_EP_S)
         ENDDO
         ENDDO
         ENDDO

      ELSE
! if cartesian grid
         CALL CALC_CG_TAU_U_G(lTAU_U_G)
      ENDIF

      call send_recv(ltau_u_g,2)

      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE CALC_TAU_U_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_CG_Tau_U_g                                         C
!  Purpose: Cross terms in the gradient of stress in U_g momentum      C
!  based on cartesian grid cut cell.                                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_TAU_U_G(lTAu_U_G)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3
      USE param1, only: zero, half
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc, only: dil_ep_s
      USE geometry
      USE indices
      USE compar
      USE fun_avg
      USE functions

      USE bc
      USE quadric
      USE cutcell
      IMPLICIT NONE
! Dummy arguments
!---------------------------------------------------------------------//
! TAU_U_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_U_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IP, JM, KM
      INTEGER :: IJK, IJKE, IJKN, IJKS, IJKT, IJKB
      INTEGER :: IJKNE, IJKSE, IJKTE, IJKBE
      INTEGER :: IPJK, IMJK, IJMK, IJKM
      INTEGER :: IPJMK, IPJKM
      integer :: itmp, jtmp, ktmp
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average viscosity
      DOUBLE PRECISION :: MU_Ge, MU_gbe
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! Cartesian grid variables
      DOUBLE PRECISION :: DEL_H, Nx, Ny, Nz
      LOGICAL :: V_NODE_AT_NE, V_NODE_AT_NW, V_NODE_AT_SE, V_NODE_AT_SW
      LOGICAL :: W_NODE_AT_TE, W_NODE_AT_TW, W_NODE_AT_BE, W_NODE_AT_BW
      DOUBLE PRECISION :: dvdx_at_N, dvdx_at_S
      DOUBLE PRECISION :: dwdx_at_T, dwdx_at_B
      DOUBLE PRECISION :: Xi, Yi, Zi, Vi, Wi, Sx, Sy, Sz
      DOUBLE PRECISION :: MU_G_CUT, SSY_CUT, SSZ_CUT
      DOUBLE PRECISION :: UW_g, VW_g, WW_g
      INTEGER :: BCV
      INTEGER :: BCT
!---------------------------------------------------------------------//

        DO K = kstart3, kend3
         DO J = jstart3, jend3
          DO I = istart3, iend3

            IJK  = FUNIJK(i,j,k)
            IJKE = FUNIJK(ieast(i,j,k),j,k)

            EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)
            IF ( .NOT.IP_AT_E(IJK) .AND. EPGA>DIL_EP_S) THEN

               IP = IP1(I)
               JM = JM1(J)
               KM = KM1(K)
   
               IPJK = FUNIJK(iplus(i,j,k),j,k)
               IMJK = FUNIJK(iminus(i,j,k),j,k)
               IJMK = FUNIJK(i,jminus(i,j,k),k)
               IJKM = FUNIJK(i,j,kminus(i,j,k))
   
               itmp = iplus(i,j,k)
               IPJMK = FUNIJK(itmp,jminus(itmp,j,k),k)
   
               ktmp = kminus(i,j,k)
               IPJKM = FUNIJK(iplus(i,j,ktmp),j,ktmp)

               jtmp = jnorth(i,j,k)
               IJKN  = FUNIJK(i,jtmp,k)
               IJKNE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

               jtmp = jsouth(i,j,k)
               IJKS  = FUNIJK(i,jtmp,k)
               IJKSE = FUNIJK(ieast(i,jtmp,k),jtmp,k)

               ktmp = ktop(i,j,k)
               IJKT  = FUNIJK(i,j,ktmp)
               IJKTE = FUNIJK(ieast(i,j,ktmp),j,ktmp)

               ktmp = kbot(i,j,k)
               IJKB  = FUNIJK(i,j,ktmp)
               IJKBE = FUNIJK(ieast(i,j,ktmp),j,ktmp)

! bulk viscosity term
            SBV = (LAMBDA_G(IJKE)*TRD_G(IJKE)) * AYZ_U(IJK) - &
                  (LAMBDA_G(IJK) *TRD_G(IJK) ) * AYZ_U(IMJK)

! shear stress terms
            IF(.NOT.CUT_U_CELL_AT(IJK)) THEN
               SSX = MU_G(IJKE)*(U_G(IPJK)-U_G(IJK))*&
                        ONEoDX_E_U(IJK)*AYZ_U(IJK) - &
                     MU_G(IJK)*(U_G(IJK)-U_G(IMJK))*&
                        ONEoDX_E_U(IMJK)*AYZ_U(IMJK)

               SSY = AVG_X_H(AVG_Y_H(MU_G(IJK),MU_G(IJKN),J),&
                             AVG_Y_H(MU_G(IJKE),MU_G(IJKNE),J),I)*&
                     (V_G(IPJK)-V_G(IJK))*ONEoDX_E_V(IJK)*AXZ_U(IJK) - &
                     AVG_X_H(AVG_Y_H(MU_G(IJKS),MU_G(IJK),JM),&
                             AVG_Y_H(MU_G(IJKSE),MU_G(IJKE),JM),I)*&
                     (V_G(IPJMK)-V_G(IJMK))*ONEoDX_E_V(IJMK)*AXZ_U(IJMK)

               IF(DO_K) THEN
                  MU_GE = AVG_X_H(AVG_Z_H(MU_G(IJK),MU_G(IJKT),K),&
                                     AVG_Z_H(MU_G(IJKE),MU_G(IJKTE),K),I)
                  MU_GBE = AVG_X_H(AVG_Z_H(MU_G(IJKB),MU_G(IJK),KM),&
                                   AVG_Z_H(MU_G(IJKBE),MU_G(IJKE),KM),I)
                  SSZ = MU_GE*(W_G(IPJK)-W_G(IJK))*&
                           ONEoDX_E_W(IJK)*AXY_U(IJK) - &
                        MU_GBE*(W_G(IPJKM)-W_G(IJKM))*&
                           ONEoDX_E_W(IJKM)*AXY_U(IJKM)
               ELSE
                  SSZ = ZERO
               ENDIF

! cut cell modifications
!---------------------------------------------------------------------//
            ELSE

               BCV = BC_U_ID(IJK)
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               SELECT CASE (BCT)
                  CASE (CG_NSW,CG_MI)
                     CUT_TAU_UG = .TRUE.
                     NOC_UG     = .TRUE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE (CG_FSW)
                     CUT_TAU_UG = .FALSE.
                     NOC_UG     = .FALSE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE(CG_PSW)
                     IF(BC_HW_G(BC_U_ID(IJK))==UNDEFINED) THEN   ! same as NSW
                        CUT_TAU_UG = .TRUE.
                        NOC_UG     = .TRUE.
                        UW_g = BC_UW_G(BCV)
                        VW_g = BC_VW_G(BCV)
                        WW_g = BC_WW_G(BCV)
                     ELSEIF(BC_HW_G(BC_U_ID(IJK))==ZERO) THEN   ! same as FSW
                        CUT_TAU_UG = .FALSE.
                        NOC_UG     = .FALSE.
                     ELSE                              ! partial slip
                        CUT_TAU_UG = .FALSE.
                        NOC_UG     = .FALSE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     ENDIF
                  CASE (NONE)
                     lTAU_U_G(IJK) = ZERO
                     CYCLE
               END SELECT

               IF(CUT_TAU_UG) THEN
                  MU_G_CUT = (VOL(IJK)*MU_G(IJK) + &
                               VOL(IPJK)*MU_G(IJKE))/&
                              (VOL(IJK) + VOL(IPJK))
               ELSE
                  MU_G_CUT = ZERO
               ENDIF

! SSX:
               CALL GET_DEL_H(IJK,'U_MOMENTUM', X_U(IJK), &
                  Y_U(IJK), Z_U(IJK), Del_H, Nx, Ny, Nz)

               SSX = MU_G(IJKE)*(U_G(IPJK)-U_G(IJK))*&
                        ONEoDX_E_U(IJK)*AYZ_U(IJK) - &
                     MU_G(IJK)*(U_G(IJK)-U_G(IMJK))*&
                        ONEoDX_E_U(IMJK)*AYZ_U(IMJK) - &
                     MU_G_CUT * (U_g(IJK) - UW_g) / DEL_H * &
                        (Nx**2) * Area_U_CUT(IJK)

! SSY:
               V_NODE_AT_NE = ((.NOT.BLOCKED_V_CELL_AT(IPJK)).AND.&
                               (.NOT.WALL_V_AT(IPJK)))
               V_NODE_AT_NW = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.&
                               (.NOT.WALL_V_AT(IJK)))
               V_NODE_AT_SE = ((.NOT.BLOCKED_V_CELL_AT(IPJMK)).AND.&
                               (.NOT.WALL_V_AT(IPJMK)))
               V_NODE_AT_SW = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.&
                               (.NOT.WALL_V_AT(IJMK)))

               IF(V_NODE_AT_NE.AND.V_NODE_AT_NW) THEN
                  Vi = HALF * (V_G(IPJK) + V_G(IJK))
                  Xi = HALF * (X_V(IPJK) + X_V(IJK))
                  Yi = HALF * (Y_V(IPJK) + Y_V(IJK))
                  Zi = HALF * (Z_V(IPJK) + Z_V(IJK))
                  Sx = X_V(IPJK) - X_V(IJK)
                  Sy = Y_V(IPJK) - Y_V(IJK)
                  Sz = Z_V(IPJK) - Z_V(IJK)
                  CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dvdx_at_N = (V_G(IPJK)-V_G(IJK)) * ONEoDX_E_V(IJK)
                  IF(NOC_UG) dvdx_at_N = dvdx_at_N - ( (Vi-VW_g)*&
                     ONEoDX_E_V(IJK)/DEL_H * (Sy*Ny+Sz*Nz) )
               ELSE
                  dvdx_at_N =  ZERO
               ENDIF

               IF(V_NODE_AT_SE.AND.V_NODE_AT_SW) THEN
                  Vi = HALF * (V_G(IPJMK) + V_G(IJMK))
                  Xi = HALF * (X_V(IPJMK) + X_V(IJMK))
                  Yi = HALF * (Y_V(IPJMK) + Y_V(IJMK))
                  Zi = HALF * (Z_V(IPJMK) + Z_V(IJMK))
                  Sx = X_V(IPJMK) - X_V(IJMK)
                  Sy = Y_V(IPJMK) - Y_V(IJMK)
                  Sz = Z_V(IPJMK) - Z_V(IJMK)
                  CALL GET_DEL_H(IJK,'U_MOMENTUM', Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dvdx_at_S = (V_G(IPJMK)-V_G(IJMK)) * ONEoDX_E_V(IJMK)
                  IF(NOC_UG) dvdx_at_S = dvdx_at_S - ( (Vi-VW_g)*&
                     ONEoDX_E_V(IJMK)/DEL_H * (Sy*Ny+Sz*Nz) )
               ELSE
                  dvdx_at_S =  ZERO
               ENDIF

               IF(V_NODE_AT_NW) THEN
                  CALL GET_DEL_H(IJK,'U_MOMENTUM', X_V(IJK), &
                     Y_V(IJK), Z_V(IJK), Del_H, Nx, Ny, Nz)
                  SSY_CUT = -MU_G_CUT * (V_G(IJK)-VW_g)/DEL_H * &
                     (Nx*Ny) * Area_U_CUT(IJK)
               ELSE
                  SSY_CUT =  ZERO
               ENDIF

               SSY = AVG_X_H(AVG_Y_H(MU_G(IJK),MU_G(IJKN),J),&
                             AVG_Y_H(MU_G(IJKE),MU_G(IJKNE),J),I)*&
                        dvdx_at_N*AXZ_U(IJK) - &
                     AVG_X_H(AVG_Y_H(MU_G(IJKS),MU_G(IJK),JM),&
                             AVG_Y_H(MU_G(IJKSE),MU_G(IJKE),JM),I)*&
                        dvdx_at_S*AXZ_U(IJMK) + SSY_CUT

! SSZ:
               IF(DO_K) THEN
                  W_NODE_AT_TE = ((.NOT.BLOCKED_W_CELL_AT(IPJK)).AND.&
                                  (.NOT.WALL_W_AT(IPJK)))
                  W_NODE_AT_TW = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.&
                                  (.NOT.WALL_W_AT(IJK)))
                  W_NODE_AT_BE = ((.NOT.BLOCKED_W_CELL_AT(IPJKM)).AND.&
                                  (.NOT.WALL_W_AT(IPJKM)))
                  W_NODE_AT_BW = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.&
                                  (.NOT.WALL_W_AT(IJKM)))

                  IF(W_NODE_AT_TE.AND.W_NODE_AT_TW) THEN
                     Wi = HALF * (W_G(IPJK) + W_G(IJK))
                     Xi = HALF * (X_W(IPJK) + X_W(IJK))
                     Yi = HALF * (Y_W(IPJK) + Y_W(IJK))
                     Zi = HALF * (Z_W(IPJK) + Z_W(IJK))
                     Sx = X_W(IPJK) - X_W(IJK)
                     Sy = Y_W(IPJK) - Y_W(IJK)
                     Sz = Z_W(IPJK) - Z_W(IJK)
                     CALL GET_DEL_H(IJK, 'U_MOMENTUM', Xi, Yi, Zi, &
                        Del_H, Nx, Ny, Nz)

                     dwdx_at_T = (W_G(IPJK)-W_G(IJK)) * ONEoDX_E_W(IJK)
                     IF(NOC_UG) dwdx_at_T = dwdx_at_T - ( (Wi-WW_g) * &
                        ONEoDX_E_W(IJK)/DEL_H * (Sy*Ny+Sz*Nz) )
                  ELSE
                     dwdx_at_T =  ZERO
                  ENDIF

                  IF(W_NODE_AT_BE.AND.W_NODE_AT_BW) THEN
                     Wi = HALF * (W_G(IPJKM) + W_G(IJKM))
                     Xi = HALF * (X_W(IPJKM) + X_W(IJKM))
                     Yi = HALF * (Y_W(IPJKM) + Y_W(IJKM))
                     Zi = HALF * (Z_W(IPJKM) + Z_W(IJKM))
                     Sx = X_W(IPJKM) - X_W(IJKM)
                     Sy = Y_W(IPJKM) - Y_W(IJKM)
                     Sz = Z_W(IPJKM) - Z_W(IJKM)

                     CALL GET_DEL_H(IJK,'U_MOMENTUM', Xi, Yi, Zi,&
                        Del_H, Nx, Ny, Nz)

                     dwdx_at_B = (W_G(IPJKM)-W_G(IJKM)) * ONEoDX_E_W(IJKM)
                     IF(NOC_UG) dwdx_at_B = dwdx_at_B - ( (Wi-WW_g) *&
                        ONEoDX_E_W(IJKM)/DEL_H * (Sy*Ny+Sz*Nz))
                  ELSE
                     dwdx_at_B =  ZERO
                  ENDIF

                  IF(W_NODE_AT_TW) THEN
                     CALL GET_DEL_H(IJK,'U_MOMENTUM', X_W(IJK), &
                        Y_W(IJK), Z_W(IJK), Del_H, Nx, Ny, Nz)
                     SSZ_CUT = -MU_G_CUT * (W_G(IJK)-WW_g)/DEL_H * &
                        (Nx*Nz) * Area_U_CUT(IJK)
                  ELSE
                     SSZ_CUT =  ZERO
                  ENDIF

                  MU_GE = AVG_X_H(AVG_Z_H(MU_G(IJK),MU_G(IJKT),K),&
                                   AVG_Z_H(MU_G(IJKE),MU_G(IJKTE),K),I)
                  MU_GBE = AVG_X_H(AVG_Z_H(MU_G(IJKB),MU_G(IJK),KM),&
                                   AVG_Z_H(MU_G(IJKBE),MU_G(IJKE),KM),I)
                  SSZ = MU_GE*dwdx_at_T*AXY_U(IJK) - &
                        MU_GBE*dwdx_at_B*AXY_U(IJKM)  + SSZ_CUT
               ELSE
                 SSZ = ZERO
               ENDIF ! DO_K

            ENDIF  ! end if/else cut_u_cell_at

! Add the terms
            lTAU_U_G(IJK) = SBV + SSX + SSY + SSZ

         ELSE
            lTAU_U_G(IJK) = ZERO
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_CG_TAU_U_G
