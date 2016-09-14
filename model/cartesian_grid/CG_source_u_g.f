!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_U_g(A_m, B_m, IER)                           C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for U_g momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_U_G(A_M, B_M)
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
      INTEGER          I, IJK, IJKE, IPJK, IJKM, IPJKM
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
      INTEGER :: J,K,JM,IP,KM
      INTEGER :: IMJK,IJMK,IJPK,IJKC,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJKP
      INTEGER :: IJKT,IJKTE,IJKB,IJKBE
      DOUBLE PRECISION :: Ue,Uw,Un,Us,Ut,Ub
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_G_E,MU_G_W,MU_G_N,MU_G_S,MU_G_T,MU_G_B,MU_G_CUT
      DOUBLE PRECISION :: UW_g
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT

!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      IF(CG_SAFE_MODE(3)==1) RETURN


!
      M = 0
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN
!
!!!$omp    parallel do private(I, IJK, IJKE, IJKM, IPJK, IPJKM,     &
!!!$omp&                  ISV, Sdp, V0, Vpm, Vmt, Vbf,              &
!!!$omp&                  Vcf, EPMUGA, VTZA, WGE, PGE, ROGA,        &
!!!$omp&                  MUGA, ROPGA, EPGA )
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)
         IJKM = KM_OF(IJK)
         IPJK = IP_OF(IJK)
         IPJKM = IP_OF(IJKM)
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)
         IF (IP_AT_E(IJK)) THEN
!
!        do nothing
!
!       dilute flow
         ELSE IF (EPGA <= DIL_EP_S) THEN
!
!        do nothing
!
!        ELSE
         ELSEIF(INTERIOR_CELL_AT(IJK)) THEN
!
            BCV = BC_U_ID(IJK)

            IF(BCV > 0 ) THEN
               BCT = BC_TYPE(BCV)
            ELSE
               BCT = 'NONE'
            ENDIF

            SELECT CASE (BCT)
               CASE ('CG_NSW')
                  NOC_UG = .TRUE.
                  UW_g = ZERO
                  MU_G_CUT =  (VOL(IJK)*MU_G(IJK) + VOL(IPJK)*MU_G(IJKE))/(VOL(IJK) + VOL(IPJK))

                  A_M(IJK,0,M) = A_M(IJK,0,M) - MU_G_CUT * Area_U_CUT(IJK)/DELH_U(IJK)
               CASE ('CG_FSW')
                  NOC_UG = .FALSE.
                  UW_g = ZERO
               CASE('CG_PSW')
                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     NOC_UG = .TRUE.
                     UW_g = BC_UW_G(BCV)
                     MU_G_CUT =  (VOL(IJK)*MU_G(IJK) + VOL(IPJK)*MU_G(IJKE))/(VOL(IJK) + VOL(IPJK))
                     A_M(IJK,0,M) = A_M(IJK,0,M) - MU_G_CUT * Area_U_CUT(IJK)/DELH_U(IJK)
                     B_M(IJK,M) = B_M(IJK,M) - MU_G_CUT * UW_g * Area_U_CUT(IJK)/DELH_U(IJK)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     NOC_UG = .FALSE.
                     UW_g = ZERO
                  ELSE                              ! partial slip
                     NOC_UG = .FALSE.
                     UW_g = BC_UW_G(BCV)
                     MU_G_CUT =  (VOL(IJK)*MU_G(IJK) + VOL(IPJK)*MU_G(IJKE))/(VOL(IJK) + VOL(IPJK))
                     A_M(IJK,0,M) = A_M(IJK,0,M) - MU_G_CUT * Area_U_CUT(IJK)*(BC_HW_G(BCV))
                     B_M(IJK,M) = B_M(IJK,M) - MU_G_CUT * UW_g * Area_U_CUT(IJK)*(BC_HW_G(BCV))
                  ENDIF
               CASE ('NONE', 'CG_MI')
                  NOC_UG = .FALSE.
            END SELECT

            IF(NOC_UG) THEN

               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)

               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IPJK = IP_OF(IJK)
               IJPK = JP_OF(IJK)
               IJKP = KP_OF(IJK)

               Ue = Theta_Ue_bar(IJK)  * U_g(IJK)  + Theta_Ue(IJK)  * U_g(IPJK)
               Uw = Theta_Ue_bar(IMJK) * U_g(IMJK) + Theta_Ue(IMJK) * U_g(IJK)

               Un = Theta_Un_bar(IJK)  * U_g(IJK)  + Theta_Un(IJK)  * U_g(IJPK)
               Us = Theta_Un_bar(IJMK) * U_g(IJMK) + Theta_Un(IJMK) * U_g(IJK)


               IJKE = EAST_OF(IJK)

               IF (WALL_AT(IJK)) THEN
                  IJKC = IJKE
               ELSE
                  IJKC = IJK
               ENDIF
               IP = IP1(I)
               IJKN = NORTH_OF(IJK)
               IJKNE = EAST_OF(IJKN)

               JM = JM1(J)
               IPJMK = IP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKSE = EAST_OF(IJKS)

               KM = KM1(K)

               IJKT = TOP_OF(IJK)
               IJKTE = EAST_OF(IJKT)

               IJKB = BOTTOM_OF(IJK)
               IJKBE = EAST_OF(IJKB)

               MU_G_E = MU_G(IJKE)
               MU_G_W = MU_G(IJKC)

               MU_G_N = AVG_X_H(AVG_Y_H(MU_G(IJKC),MU_G(IJKN),J),&
                                 AVG_Y_H(MU_G(IJKE),MU_G(IJKNE),J),I)

               MU_G_S = AVG_X_H(AVG_Y_H(MU_G(IJKS),MU_G(IJKC),JM),&
                                 AVG_Y_H(MU_G(IJKSE),MU_G(IJKE),JM),I)

               B_NOC =     MU_G_E * Ayz_U(IJK)  * (Ue-UW_g) * NOC_U_E(IJK)  *2.0d0&
                       -   MU_G_W * Ayz_U(IMJK) * (Uw-UW_g) * NOC_U_E(IMJK) *2.0d0&
                       +   MU_G_N * Axz_U(IJK)  * (Un-UW_g) * NOC_U_N(IJK)  &
                       -   MU_G_S * Axz_U(IJMK) * (Us-UW_g) * NOC_U_N(IJMK)

               IF(DO_K) THEN

                  Ut = Theta_Ut_bar(IJK)  * U_g(IJK)  + Theta_Ut(IJK)  * U_g(IJKP)
                  Ub = Theta_Ut_bar(IJKM) * U_g(IJKM) + Theta_Ut(IJKM) * U_g(IJK)

                  MU_G_T = AVG_X_H(AVG_Z_H(MU_G(IJKC),MU_G(IJKT),K),&
                            AVG_Z_H(MU_G(IJKE),MU_G(IJKTE),K),I)

                  MU_G_B = AVG_X_H(AVG_Z_H(MU_G(IJKB),MU_G(IJKC),KM),&
                            AVG_Z_H(MU_G(IJKBE),MU_G(IJKE),KM),I)

                  B_NOC = B_NOC  +   MU_G_T * Axy_U(IJK)  * (Ut-UW_g) * NOC_U_T(IJK)  &
                                 -   MU_G_B * Axy_U(IJKM) * (Ub-UW_g) * NOC_U_T(IJKM)
               ENDIF

               B_M(IJK,M) = B_M(IJK,M)   +  B_NOC


            ENDIF

         ENDIF

      END DO


      RETURN
      END SUBROUTINE CG_SOURCE_U_G
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_U_g_BC(A_m, B_m, IER)                        C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for U_g momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_U_G_BC(A_M, B_M)
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
      INTEGER          I,  J, K, IM, IJK,&
                       JM, KM, IJKW, IMJK, IJMK,IJKM
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

      IF(CG_SAFE_MODE(3)==1) RETURN

      M = 0

      DO IJK = ijkstart3, ijkend3

         BCV = BC_U_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE(BCV)
         ELSE
            BCT = 'NONE'
         ENDIF


         SELECT CASE (BCT)

            CASE ('CG_NSW')

               IF(WALL_U_AT(IJK)) THEN

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

              IF(WALL_U_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE

!                  B_M(IJK,M) = - U_g(U_MASTER_OF(IJK))  ! Velocity of master node

                  B_M(IJK,M) = ZERO

                  IF(DABS(NORMAL_U(IJK,1))/=ONE) THEN

                     IF (U_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                        A_M(IJK,E,M) = ONE
                     ELSEIF (U_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                        A_M(IJK,W,M) = ONE
                     ELSEIF (U_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                        A_M(IJK,N,M) = ONE
                     ELSEIF (U_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                        A_M(IJK,S,M) = ONE
                     ELSEIF (U_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                        A_M(IJK,T,M) = ONE
                     ELSEIF (U_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                        A_M(IJK,B,M) = ONE
                     ENDIF

                  ENDIF

               ENDIF


            CASE ('CG_PSW')

               IF(WALL_U_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE


                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = -BC_UW_G(BCV)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     B_M(IJK,M) = ZERO

                     IF(DABS(NORMAL_U(IJK,1))/=ONE) THEN

                        IF (U_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                           A_M(IJK,E,M) = ONE
                        ELSEIF (U_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                           A_M(IJK,W,M) = ONE
                        ELSEIF (U_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                           A_M(IJK,N,M) = ONE
                        ELSEIF (U_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                           A_M(IJK,S,M) = ONE
                        ELSEIF (U_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                           A_M(IJK,T,M) = ONE
                        ELSEIF (U_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                           A_M(IJK,B,M) = ONE
                        ENDIF

                     ENDIF

                  ELSE                              ! partial slip

                     B_M(IJK,M) = ZERO

                     I = I_OF(IJK)
                     J = J_OF(IJK)
                     K = K_OF(IJK)

                     IM = I - 1
                     JM = J - 1
                     KM = K - 1

                     IMJK = FUNIJK(IM,J,K)
                     IJMK = FUNIJK(I,JM,K)
                     IJKM = FUNIJK(I,J,KM)


                     IF(DABS(NORMAL_U(IJK,1))/=ONE) THEN

                        IF (U_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDX_E_U(IJK))
                           A_M(IJK,E,M) = -(HALF*BC_HW_G(BCV)-ONEoDX_E_U(IJK))
                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)

!                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODX_E(I))
!                           A_M(IJK,E,M) = -(HALF*BC_HW_G(BCV)-ODX_E(I))
!                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)
!                           print*,'ug master at east'
                        ELSEIF (U_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                           A_M(IJK,W,M) = -(HALF*BC_HW_G(BCV)-ONEoDX_E_U(IMJK))
                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDX_E_U(IMJK))
                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)

!                           A_M(IJK,W,M) = -(HALF*BC_HW_G(BCV)-ODX_E(IM))
!                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODX_E(IM))
!                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)
!                           print*,'ug master at west'
                        ELSEIF (U_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN

                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDY_N_U(IJK))
                           A_M(IJK,N,M) = -(HALF*BC_HW_G(BCV)-ONEoDY_N_U(IJK))
                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)

!                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODY_N(J))
!                           A_M(IJK,N,M) = -(HALF*BC_HW_G(BCV)-ODY_N(J))
!                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)
                        ELSEIF (U_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                           A_M(IJK,S,M) = -(HALF*BC_HW_G(BCV)-ONEoDY_N_U(IJMK))
                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDY_N_U(IJMK))
                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)

!                           A_M(IJK,S,M) = -(HALF*BC_HW_G(BCV)-ODY_N(JM))
!                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODY_N(JM))
!                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)
                        ELSEIF (U_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                              A_M(IJK,0,M)=-(HALF*BC_HW_G(BCV)+ONEoDZ_T_U(IJK))
                              A_M(IJK,T,M)=-(HALF*BC_HW_G(BCV)-ONEoDZ_T_U(IJK))
                              B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)

!                              A_M(IJK,0,M)=-(HALF*BC_HW_G(L)+ODZ_T(K)*OX_E(I))
!                              A_M(IJK,T,M)=-(HALF*BC_HW_G(L)-ODZ_T(K)*OX_E(I))
!                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                        ELSEIF (U_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                              A_M(IJK,B,M) = -(HALF*BC_HW_G(BCV)-ONEoDZ_T_U(IJKM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDZ_T_U(IJKM))
                              B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV)

!                              A_M(IJK,B,M) = -(HALF*BC_HW_G(L)-ODZ_T(KM)*OX_E(I&
!                                 ))
!                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(KM)*OX_E(I&
!                                 ))
!                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                        ENDIF

                     ENDIF

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

               IF(BC_U_g(BCV)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_U_g(BCV)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_g(BCV)*NORMAL_U(IJK,1)
               ENDIF


               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJKW,E,M) = ZERO
                  A_M(IJKW,W,M) = ZERO
                  A_M(IJKW,N,M) = ZERO
                  A_M(IJKW,S,M) = ZERO
                  A_M(IJKW,T,M) = ZERO
                  A_M(IJKW,B,M) = ZERO
                  A_M(IJKW,0,M) = -ONE

                  IF(BC_U_g(BCV)/=UNDEFINED) THEN
                     B_M(IJKW,M) = - BC_U_g(BCV)
                  ELSE
                     B_M(IJKW,M) = - BC_VELMAG_g(BCV)*NORMAL_U(IJK,1)
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

               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJK,W,M) = ONE
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

               IF(BC_U_g(BCV)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_U_g(BCV)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_g(BCV)*NORMAL_S(IJK,1)
               ENDIF


               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJKW,E,M) = ZERO
                  A_M(IJKW,W,M) = ZERO
                  A_M(IJKW,N,M) = ZERO
                  A_M(IJKW,S,M) = ZERO
                  A_M(IJKW,T,M) = ZERO
                  A_M(IJKW,B,M) = ZERO
                  A_M(IJKW,0,M) = -ONE

                  IF(BC_U_g(BCV)/=UNDEFINED) THEN
                     B_M(IJKW,M) = - BC_U_g(BCV)
                  ELSE
                     B_M(IJKW,M) = - BC_VELMAG_g(BCV)*NORMAL_S(IJK,1)
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


               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJK,W,M) = ONE
                  A_M(IJK,0,M) = -ONE

               ENDIF

         END SELECT

      ENDDO


      RETURN

!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      END SUBROUTINE CG_SOURCE_U_G_BC

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3
!// 360 Check if i,j,k resides on current processor
