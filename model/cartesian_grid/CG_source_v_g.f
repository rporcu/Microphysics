!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_V_g(A_m, B_m, IER)                              C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for V_g momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_V_G(A_M, B_M)
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
      use matrix, only: e, w, s, n, t, b
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE indices
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
      INTEGER          I, J, K, IJK, IJKN
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
      DOUBLE PRECISION B_m(DIMENSION_3)

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER ::          IM,JP,KM
      INTEGER ::          IMJK,IPJK,IJMK,IJPK,IJKP,IJKM,IJKC,IJKE,IJKNE,IJKW,IJKWN,IMJPK
      INTEGER ::          IJKT,IJKTN,IJKB,IJKBN
      DOUBLE PRECISION :: Vn,Vs,Ve,Vw, Vt,Vb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_G_E,MU_G_W,MU_G_N,MU_G_S,MU_G_T,MU_G_B,MU_G_CUT
      DOUBLE PRECISION :: VW_g
      INTEGER :: BCV
      INTEGER :: BCT

!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      IF(CG_SAFE_MODE(4)==1) RETURN

!
      M = 0
      IF (.NOT.MOMENTUM_Y_EQ(0)) RETURN
!
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IJKN = NORTH_OF(IJK)
         EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)
         IF (IP_AT_N(IJK)) THEN
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
            BCV = BC_V_ID(IJK)

            IF(BCV > 0 ) THEN
               BCT = BC_TYPE_ENUM(BCV)
            ELSE
               BCT = NONE
            ENDIF

            SELECT CASE (BCT)
               CASE (CG_NSW)
                  NOC_VG = .TRUE.
                  VW_g = ZERO
                  MU_G_CUT = (VOL(IJK)*MU_G(IJK) + VOL(IJKN)*MU_G(IJKN))/(VOL(IJK) + VOL(IJKN))

                  A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_G_CUT * Area_V_CUT(IJK)/DELH_V(IJK)
               CASE (CG_FSW)
                  NOC_VG = .FALSE.
                  VW_g = ZERO
               CASE(CG_PSW)
                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     NOC_VG = .TRUE.
                     VW_g = BC_VW_G(BCV)
                     MU_G_CUT = (VOL(IJK)*MU_G(IJK) + VOL(IJKN)*MU_G(IJKN))/(VOL(IJK) + VOL(IJKN))
                     A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_G_CUT * Area_V_CUT(IJK)/DELH_V(IJK)
                     B_M(IJK) = B_M(IJK) - MU_G_CUT * VW_g * Area_V_CUT(IJK)/DELH_V(IJK)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     NOC_VG = .FALSE.
                     VW_g = ZERO
                  ELSE                              ! partial slip
                     NOC_VG = .FALSE.
                     VW_g = BC_VW_G(BCV)
                     MU_G_CUT = (VOL(IJK)*MU_G(IJK) + VOL(IJKN)*MU_G(IJKN))/(VOL(IJK) + VOL(IJKN))
                     A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_G_CUT * Area_V_CUT(IJK)*(BC_HW_G(BCV))
                     B_M(IJK) = B_M(IJK) - MU_G_CUT * VW_g * Area_V_CUT(IJK)*(BC_HW_G(BCV))
                  ENDIF
               CASE (NONE, CG_MI)
                  NOC_VG = .FALSE.
               CASE DEFAULT
                  STOP __LINE__
            END SELECT


            IF(NOC_VG) THEN

               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IPJK = IP_OF(IJK)
               IJPK = JP_OF(IJK)
               IJKP = KP_OF(IJK)

               Vn = Theta_Vn_bar(IJK)  * V_g(IJK)  + Theta_Vn(IJK)  * V_g(IJPK)
               Vs = Theta_Vn_bar(IJMK) * V_g(IJMK) + Theta_Vn(IJMK) * V_g(IJK)

               Ve = Theta_Ve_bar(IJK)  * V_g(IJK)  + Theta_Ve(IJK)  * V_g(IPJK)
               Vw = Theta_Ve_bar(IMJK) * V_g(IMJK) + Theta_Ve(IMJK) * V_g(IJK)


               IJKN = NORTH_OF(IJK)
               IF (WALL_AT(IJK)) THEN
                  IJKC = IJKN
               ELSE
                  IJKC = IJK
               ENDIF
               JP = JP1(J)
               IJKE = EAST_OF(IJK)
               IJKNE = EAST_OF(IJKN)

               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWN = NORTH_OF(IJKW)
               IMJPK = JP_OF(IMJK)


               KM = KM1(K)
               IJKT = TOP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)

               IJKB = BOTTOM_OF(IJK)
               IJKBN = NORTH_OF(IJKB)


               MU_G_E = AVG_Y_H(AVG_X_H(MU_G(IJKC),MU_G(IJKE),I),&
                         AVG_X_H(MU_G(IJKN),MU_G(IJKNE),I),J)

               MU_G_W = AVG_Y_H(AVG_X_H(MU_G(IJKW),MU_G(IJKC),IM),&
                         AVG_X_H(MU_G(IJKWN),MU_G(IJKN),IM),J)

               MU_G_N = MU_G(IJKN)

               MU_G_S = MU_G(IJKC)


               B_NOC =     MU_G_N * Axz_V(IJK)  * (Vn-VW_g) * NOC_V_N(IJK)   *2.0d0&
                       -   MU_G_S * Axz_V(IJMK) * (Vs-VW_g) * NOC_V_N(IJMK)  *2.0d0&
                       +   MU_G_E * Ayz_V(IJK)  * (Ve-VW_g) * NOC_V_E(IJK)   &
                       -   MU_G_W * Ayz_V(IMJK) * (Vw-VW_g) * NOC_V_E(IMJK)

               IF(DO_K) THEN

                  Vt = Theta_Vt_bar(IJK)  * V_g(IJK)  + Theta_Vt(IJK)  * V_g(IJKP)
                  Vb = Theta_Vt_bar(IJKM) * V_g(IJKM) + Theta_Vt(IJKM) * V_g(IJK)

                  MU_G_T = AVG_Y_H(AVG_Z_H(MU_G(IJKC),MU_G(IJKT),K),&
                            AVG_Z_H(MU_G(IJKN),MU_G(IJKTN),K),J)

                  MU_G_B = AVG_Y_H(AVG_Z_H(MU_G(IJKB),MU_G(IJKC),KM),&
                            AVG_Z_H(MU_G(IJKBN),MU_G(IJKN),KM),J)

                  B_NOC = B_NOC + MU_G_T * Axy_V(IJK)  * (Vt-VW_g) * NOC_V_T(IJK)   &
                                - MU_G_B * Axy_V(IJKM) * (Vb-VW_g) * NOC_V_T(IJKM)

               ENDIF

                  B_M(IJK) = B_M(IJK) + B_NOC

            ENDIF


         ENDIF
      END DO
!
      RETURN
      END SUBROUTINE CG_SOURCE_V_G


!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_V_g_BC(A_m, B_m, IER)                        C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for V_g momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_V_G_BC(A_M, B_M)
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
      use matrix, only: e, w, s, n, t, b
      USE scales
      USE constant
      USE physprop
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
      INTEGER          I,  J, K, JM, IJK,&
                       IM, KM, IJKS, IMJK
!
!                      Solids phase
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER :: BCV
      INTEGER :: BCT

!-----------------------------------------------

      IF(CG_SAFE_MODE(4)==1) RETURN

!
      M = 0

      DO IJK = ijkstart3, ijkend3

         BCV = BC_V_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE_ENUM(BCV)
         ELSE
            BCT = NONE
         ENDIF

         SELECT CASE (BCT)


            CASE (CG_NSW)
               IF(WALL_V_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE

                  B_M(IJK) = ZERO



               ENDIF


            CASE (CG_FSW)
               IF(WALL_V_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE

!                  B_M(IJK) = - V_g(V_MASTER_OF(IJK))    ! Velocity of master node

                  B_M(IJK) = ZERO

                  IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN

                     IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                        A_M(IJK,E,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                        A_M(IJK,W,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                        A_M(IJK,N,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                        A_M(IJK,S,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                        A_M(IJK,T,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                        A_M(IJK,B,M) = ONE
                     ENDIF

                  ENDIF

               ENDIF



            CASE (CG_PSW)
               IF(WALL_V_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE


                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK) = -BC_VW_G(BCV)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     B_M(IJK) = ZERO

                     IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN

                        IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                           A_M(IJK,E,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                           A_M(IJK,W,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                           A_M(IJK,N,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                           A_M(IJK,S,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                           A_M(IJK,T,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                           A_M(IJK,B,M) = ONE
                        ENDIF

                     ENDIF

                  ELSE                              ! partial slip

                     B_M(IJK) = ZERO

                     I = I_OF(IJK)
                     J = J_OF(IJK)
                     K = K_OF(IJK)

                     IM = I - 1
                     JM = J - 1
                     KM = K - 1

                     IMJK = FUNIJK(IM,J,K)

                     IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN

                        IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
!                          A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODX_E(I))
!                          A_M(IJK,E,M) = -(HALF*BC_HW_G(BCV)-ODX_E(I))
!                          B_M(IJK) = -BC_HW_G(BCV)*BC_VW_G(BCV)

                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDX_E_V(IJK))
                           A_M(IJK,E,M) = -(HALF*BC_HW_G(BCV)-ONEoDX_E_V(IJK))
                           B_M(IJK) = -BC_HW_G(BCV)*BC_VW_G(BCV)

                        ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
!                          A_M(IJK,W,M) = -(HALF*BC_HW_G(BCV)-ODX_E(IM))
!                          A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODX_E(IM))
!                          B_M(IJK) = -BC_HW_G(BCV)*BC_VW_G(BCV)

                           A_M(IJK,W,M) = -(HALF*BC_HW_G(BCV)-ONEoDX_E_V(IMJK))
                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDX_E_V(IMJK))
                           B_M(IJK) = -BC_HW_G(BCV)*BC_VW_G(BCV)


                        ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN

                              A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODY_N(J))
                              A_M(IJK,N,M) = -(HALF*BC_HW_G(BCV)-ODY_N(J))
                              B_M(IJK) = -BC_HW_G(BCV)*BC_UW_G(BCV)

!                           print*,'vg master at north'
                        ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN

                              A_M(IJK,S,M) = -(HALF*BC_HW_G(BCV)-ODY_N(JM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODY_N(JM))
                              B_M(IJK) = -BC_HW_G(BCV)*BC_UW_G(BCV)

!                           print*,'vg master at south'
                        ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                           A_M(IJK,T,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                           A_M(IJK,B,M) = ONE
                        ENDIF

                     ENDIF




                  ENDIF

               ENDIF




            CASE (CG_MI)
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE

               IF(BC_V_g(BCV)/=UNDEFINED) THEN
                  B_M(IJK) = - BC_V_g(BCV)
               ELSE
                  B_M(IJK) = - BC_VELMAG_g(BCV)*NORMAL_V(IJK,2)
               ENDIF


               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJKS,E,M) = ZERO
                  A_M(IJKS,W,M) = ZERO
                  A_M(IJKS,N,M) = ZERO
                  A_M(IJKS,S,M) = ZERO
                  A_M(IJKS,T,M) = ZERO
                  A_M(IJKS,B,M) = ZERO
                  A_M(IJKS,0,M) = -ONE

                  IF(BC_V_g(BCV)/=UNDEFINED) THEN
                     B_M(IJKS) = - BC_V_g(BCV)
                  ELSE
                     B_M(IJKS) = - BC_VELMAG_g(BCV)*NORMAL_V(IJK,2)
                  ENDIF


               ENDIF


            CASE (CG_PO)
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK) = ZERO

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJK,S,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK) = ZERO

               ENDIF


         END SELECT


         BCV = BC_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE_ENUM(BCV)
         ELSE
            BCT = NONE
         ENDIF

         SELECT CASE (BCT)


            CASE (CG_MI)
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE

               IF(BC_V_g(BCV)/=UNDEFINED) THEN
                  B_M(IJK) = - BC_V_g(BCV)
               ELSE
                  B_M(IJK) = - BC_VELMAG_g(BCV)*NORMAL_S(IJK,2)
               ENDIF


               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJKS,E,M) = ZERO
                  A_M(IJKS,W,M) = ZERO
                  A_M(IJKS,N,M) = ZERO
                  A_M(IJKS,S,M) = ZERO
                  A_M(IJKS,T,M) = ZERO
                  A_M(IJKS,B,M) = ZERO
                  A_M(IJKS,0,M) = -ONE

                  IF(BC_V_g(BCV)/=UNDEFINED) THEN
                     B_M(IJKS) = - BC_V_g(BCV)
                  ELSE
                     B_M(IJKS) = - BC_VELMAG_g(BCV)*NORMAL_S(IJK,2)
                  ENDIF


               ENDIF

            CASE (CG_PO)

               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK) = ZERO

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJK,S,M) = ONE
                  A_M(IJK,0,M) = -ONE

               ENDIF

         END SELECT

      ENDDO


      RETURN
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      END SUBROUTINE CG_SOURCE_V_G_BC

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3
!// 360 Check if i,j,k resides on current processor
