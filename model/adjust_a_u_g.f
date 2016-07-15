!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_A_U_g(A_m, B_m, IER)                            C
!  Author: M. Syamlal                                 Date:  2-AUG-96  C
!                                                                      C
!  Purpose: Handle the special case of the center coefficient in       C
!  U_g momentum eq. becoming zero.                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ADJUST_A_U_G(A_M, B_M)

      USE param
      USE param1
      USE fldvar
      USE geometry
      USE run
      USE indices
      USE usr
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      use matrix, only: e, w, s, n, t, b

      IMPLICIT NONE

!                      Indices
      INTEGER          I, IP, IJK, IJKE, IMJK
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------

      M = 0
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN


      DO IJK = ijkstart3, ijkend3
         IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            IF (B_M(IJK,M) < ZERO) THEN
               IJKE = EAST_OF(IJK)
               IP = IP1(I_OF(IJK))
               IF (ROP_G(IJKE)*AYZ_U(IJK) > SMALL_NUMBER) THEN
                  B_M(IJK,M) = SQRT((-B_M(IJK,M)/(ROP_G(IJKE)*AVG_X_E(ONE,ZERO,&
                     IP)*AYZ_U(IJK))))
               ELSE
                  B_M(IJK,M) = ZERO
               ENDIF
            ELSE IF (B_M(IJK,M) > ZERO) THEN
               I = I_OF(IJK)
               IMJK = IM_OF(IJK)
               IF (ROP_G(IJK)*AYZ_U(IMJK) > SMALL_NUMBER) THEN
                  B_M(IJK,M) = SQRT(B_M(IJK,M)/(ROP_G(IJK)*AVG_X_E(ZERO,ONE,I)*&
                     AYZ_U(IMJK)))
               ELSE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDIF
         ENDIF
      END DO

      RETURN
      END SUBROUTINE ADJUST_A_U_G
