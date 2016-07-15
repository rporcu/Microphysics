!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_W_g(A_m, B_m, IER)                           C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for W_g momentum eq.                                                C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CG_SOURCE_W_G(A_M, B_M)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      use ambm, only: e, w, s, n, t, b
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE bc
      USE compar
      USE sendrecv
      USE drag
      USE fun_avg
      USE functions
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
      USE quadric
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IJKT, IMJK, IJKP, IMJKP,&
                       IJKE, IJKW, IJKTE, IM, IPJK
!
!                      Phase index
      INTEGER          M
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER :: JM,IP,JP,IJMK,IJPK,IJKC,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJKM,KM,KP,IJMKP
      INTEGER :: IJKTN,IJKWT,IJKST
      DOUBLE PRECISION :: We,Ww,Wn,Ws,Wt,Wb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_G_E,MU_G_W,MU_G_N,MU_G_S,MU_G_T,MU_G_B,MU_G_CUT
      DOUBLE PRECISION :: WW_g
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!                       virtual (added) mass
      DOUBLE PRECISION ROP_MA, U_se, Usw, Wse, Wsw, Wsn, Wss, Wst, Wsb, Usc,Vsc,Vsn,Vss
! Wall function
      DOUBLE PRECISION :: W_F_Slip
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!-----------------------------------------------

      IF(CG_SAFE_MODE(5)==1) RETURN

!
      M = 0
      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN
!
!
!!!$omp  parallel do private( I, J, K, IJK, IJKT, ISV, Sdp, V0, Vpm, Vmt, Vbf, &
!!!$omp&  PGT, ROGA, IMJK, IJKP, IMJKP, IJKW, IJKTE, IJKTW, IM, IPJK,  &
!!!$omp&  CTE, CTW, SXZB, EPMUOX, VXZA, VXZB, UGT, VCOA, VCOB, IJKE,&
!!!$omp&  MUGA, ROPGA, EPGA, LINE)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IJKT = TOP_OF(IJK)
         EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)
         IF (IP_AT_T(IJK)) THEN
!
!        do nothing
!
!       dilute flow
         ELSE IF (EPGA <= DIL_EP_S) THEN
!
!        do nothing
!
!         ELSE
         ELSEIF(INTERIOR_CELL_AT(IJK)) THEN
!
            BCV = BC_W_ID(IJK)

            IF(BCV > 0 ) THEN
               BCT = BC_TYPE(BCV)
            ELSE
               BCT = 'NONE'
            ENDIF

            SELECT CASE (BCT)
               CASE ('CG_NSW')
                  NOC_WG = .TRUE.
                  WW_g = ZERO
                  MU_G_CUT =  (VOL(IJK)*MU_G(IJK) + VOL(IJKT)*MU_G(IJKT))/(VOL(IJK) + VOL(IJKT))

                  A_M(IJK,0,M) = A_M(IJK,0,M) - MU_G_CUT * Area_W_CUT(IJK)/DELH_W(IJK)

               CASE ('CG_FSW')
                  NOC_WG = .FALSE.
                  WW_g = ZERO
               CASE('CG_PSW')
                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     NOC_WG = .TRUE.
                     WW_g = BC_WW_G(BCV)
                     MU_G_CUT = (VOL(IJK)*MU_G(IJK) + VOL(IJKT)*MU_G(IJKT))/(VOL(IJK) + VOL(IJKT))
                     A_M(IJK,0,M) = A_M(IJK,0,M) - MU_G_CUT * Area_W_CUT(IJK)/DELH_W(IJK)
                     B_M(IJK,M) = B_M(IJK,M) - MU_G_CUT * WW_g * Area_W_CUT(IJK)/DELH_W(IJK)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     NOC_WG = .FALSE.
                     WW_g = ZERO
                  ELSE                              ! partial slip
                     NOC_WG = .FALSE.
                     WW_g = BC_WW_G(BCV)
                     MU_G_CUT = (VOL(IJK)*MU_G(IJK) + VOL(IJKT)*MU_G(IJKT))/(VOL(IJK) + VOL(IJKT))
                     A_M(IJK,0,M) = A_M(IJK,0,M) - MU_G_CUT * Area_W_CUT(IJK)*(BC_HW_G(BCV))
                     B_M(IJK,M) = B_M(IJK,M) - MU_G_CUT * WW_g * Area_W_CUT(IJK)*(BC_HW_G(BCV))

                  ENDIF
               CASE ('NONE')
                  NOC_WG = .FALSE.
            END SELECT

            IF(NOC_WG) THEN

               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IPJK = IP_OF(IJK)
               IJPK = JP_OF(IJK)
               IJKP = KP_OF(IJK)

               We = Theta_We_bar(IJK)  * W_g(IJK)  + Theta_We(IJK)  * W_g(IPJK)
               Ww = Theta_We_bar(IMJK) * W_g(IMJK) + Theta_We(IMJK) * W_g(IJK)

               Wn = Theta_Wn_bar(IJK)  * W_g(IJK)  + Theta_Wn(IJK)  * W_g(IJPK)
               Ws = Theta_Wn_bar(IJMK) * W_g(IJMK) + Theta_Wn(IJMK) * W_g(IJK)

               Wt = Theta_Wt_bar(IJK)  * W_g(IJK)  + Theta_Wt(IJK)  * W_g(IJKP)
               Wb = Theta_Wt_bar(IJKM) * W_g(IJKM) + Theta_Wt(IJKM) * W_g(IJK)

               IJKE = EAST_OF(IJK)

               ijkt = top_of(ijk)

               IF (WALL_AT(IJK)) THEN
                  IJKC = IJKT
               ELSE
                  IJKC = IJK
               ENDIF

               IP = IP1(I)
               IM = IM1(I)
               IJKN = NORTH_OF(IJK)
               IJKNE = EAST_OF(IJKN)

               JM = JM1(J)
               IPJMK = IP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKSE = EAST_OF(IJKS)

               KP = KP1(K)
               IJKT = TOP_OF(IJK)
               IJKE = EAST_OF(IJK)
               IJKP = KP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
               IJKTE = EAST_OF(IJKT)
               IJKW = WEST_OF(IJK)
               IJKWT = TOP_OF(IJKW)
               IJKS = SOUTH_OF(IJK)
               IJKST = TOP_OF(IJKS)

               MU_G_E = AVG_Z_H(AVG_X_H(MU_G(IJKC),MU_G(IJKE),I),&
                         AVG_X_H(MU_G(IJKT),MU_G(IJKTE),I),K)

               MU_G_W = AVG_Z_H(AVG_X_H(MU_G(IJKW),MU_G(IJKC),IM),&
                         AVG_X_H(MU_G(IJKWT),MU_G(IJKT),IM),K)

               MU_G_N = AVG_Z_H(AVG_Y_H(MU_G(IJKC),MU_G(IJKN),J),&
                         AVG_Y_H(MU_G(IJKT),MU_G(IJKTN),J),K)

               MU_G_S = AVG_Z_H(AVG_Y_H(MU_G(IJKS),MU_G(IJKC),JM),&
                         AVG_Y_H(MU_G(IJKST),MU_G(IJKT),JM),K)

               MU_G_T = MU_G(IJKT)
               MU_G_B = MU_G(IJKC)

               B_NOC =     MU_G_E * Ayz_W(IJK)  * (We-WW_g) * NOC_W_E(IJK)  &
                       -   MU_G_W * Ayz_W(IMJK) * (Ww-WW_g) * NOC_W_E(IMJK) &
                       +   MU_G_N * Axz_W(IJK)  * (Wn-WW_g) * NOC_W_N(IJK)  &
                       -   MU_G_S * Axz_W(IJMK) * (Ws-WW_g) * NOC_W_N(IJMK) &
                       +   MU_G_T * Axy_W(IJK)  * (Wt-WW_g) * NOC_W_T(IJK)  *2.0d0&
                       -   MU_G_B * Axy_W(IJKM) * (Wb-WW_g) * NOC_W_T(IJKM) *2.0D0

               B_M(IJK,M) = B_M(IJK,M)   +  B_NOC
            ENDIF


         ENDIF
      END DO

      RETURN
      END SUBROUTINE CG_SOURCE_W_G

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_W_g_BC(A_m, B_m, IER)                        C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for W_g momentum eq.                                                C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CG_SOURCE_W_G_BC(A_M, B_M)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      use ambm, only: e, w, s, n, t, b
      USE scales
      USE constant
      USE physprop
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

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
      USE quadric
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          IJK, IJKB
!
!                      Solids phase
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT

!-----------------------------------------------

      IF(CG_SAFE_MODE(5)==1) RETURN

!
      M = 0

      DO IJK = ijkstart3, ijkend3

         BCV = BC_W_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE(BCV)
         ELSE
            BCT = 'NONE'
         ENDIF

         SELECT CASE (BCT)

            CASE ('CG_NSW')

               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE

                  B_M(IJK,M) = ZERO

               ENDIF

            CASE ('CG_FSW')

               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE

!                  B_M(IJK,M) = - W_g(W_MASTER_OF(IJK))  ! Velocity of master node

                  B_M(IJK,M) = ZERO

                  IF(DABS(NORMAL_W(IJK,3))/=ONE) THEN

                     IF (W_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                        A_M(IJK,E,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                        A_M(IJK,W,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                        A_M(IJK,N,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                        A_M(IJK,S,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                        A_M(IJK,T,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                        A_M(IJK,B,M) = ONE
                     ENDIF

                  ENDIF

               ENDIF

            CASE ('CG_PSW')

               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE


                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = -BC_WW_G(BCV)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     B_M(IJK,M) = ZERO

                     IF(DABS(NORMAL_W(IJK,3))/=ONE) THEN

                        IF (W_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                           A_M(IJK,E,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                           A_M(IJK,W,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                           A_M(IJK,N,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                           A_M(IJK,S,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                           A_M(IJK,T,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                           A_M(IJK,B,M) = ONE
                        ENDIF

                     ENDIF

                  ELSE                              ! partial slip





                  ENDIF

               ENDIF


            CASE ('CG_MI')

               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE

               IF(BC_W_g(BCV)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_W_g(BCV)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_g(BCV)*NORMAL_W(IJK,3)
               ENDIF


               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJKB,E,M) = ZERO
                  A_M(IJKB,W,M) = ZERO
                  A_M(IJKB,N,M) = ZERO
                  A_M(IJKB,S,M) = ZERO
                  A_M(IJKB,T,M) = ZERO
                  A_M(IJKB,B,M) = ZERO
                  A_M(IJKB,0,M) = -ONE

                  IF(BC_W_g(BCV)/=UNDEFINED) THEN
                     B_M(IJKB,M) = - BC_W_g(BCV)
                  ELSE
                     B_M(IJKB,M) = - BC_VELMAG_g(BCV)*NORMAL_W(IJK,3)
                  ENDIF


               ENDIF

            CASE ('CG_PO')

               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJK,B,M) = ONE
                  A_M(IJK,0,M) = -ONE

               ENDIF

         END SELECT

         BCV = BC_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE(BCV)
         ELSE
            BCT = 'NONE'
         ENDIF

         SELECT CASE (BCT)

            CASE ('CG_MI')

               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE

               IF(BC_W_g(BCV)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_W_g(BCV)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_g(BCV)*NORMAL_S(IJK,3)
               ENDIF


               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJKB,E,M) = ZERO
                  A_M(IJKB,W,M) = ZERO
                  A_M(IJKB,N,M) = ZERO
                  A_M(IJKB,S,M) = ZERO
                  A_M(IJKB,T,M) = ZERO
                  A_M(IJKB,B,M) = ZERO
                  A_M(IJKB,0,M) = -ONE

                  IF(BC_W_g(BCV)/=UNDEFINED) THEN
                     B_M(IJKB,M) = - BC_W_g(BCV)
                  ELSE
                     B_M(IJKB,M) = - BC_VELMAG_g(BCV)*NORMAL_S(IJK,3)
                  ENDIF


               ENDIF

            CASE ('CG_PO')

               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJK,B,M) = ONE
                  A_M(IJK,0,M) = -ONE

               ENDIF

         END SELECT


      ENDDO

      RETURN

!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      END SUBROUTINE CG_SOURCE_W_G_BC

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3
!// 360 Check if i,j,k resides on current processor
