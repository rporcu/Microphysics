MODULE ur_facs

   use param, only: DIM_EQS

! Under relaxation factors for coefficient update:
!  [0]  every time step (explicit)
!  [1]  every iteration (implicit)
! (0,1) under-relaxed
   DOUBLE PRECISION :: UR_FAC(DIM_EQS)


   contains


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: UNDER_RELAX                                             !
!  Author: M. Syamlal                                 Date: 24-MAY-96  !
!                                                                      !
!  Purpose: Under-relax equation.                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE UNDER_RELAX(VAR, A_M, B_M, AXIS, EQ)

   use param, only: dimension_3
   use param1, only: one
   use compar, only: ijkstart3, ijkend3

   implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Variable
      DOUBLE PRECISION :: Var(DIMENSION_3)
! Septadiagonal matrix
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3)
!   Vector b_m
      DOUBLE PRECISION :: B_m(DIMENSION_3)
! Equation ID
      INTEGER :: EQ
! Axis ID: U, V, W, S (scalar)
      CHARACTER :: AXIS

! Local variables:
!---------------------------------------------------------------------//
! Loop index
      INTEGER :: IJK
! Functions of under-relaxation factor
      DOUBLE PRECISION :: f1, f2
! Center coefficient
      DOUBLE PRECISION :: Ap

!.......................................................................!
      F1 = ONE/UR_FAC(EQ)
      F2 = F1 - ONE

      if (axis.eq.'S') then

         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN
               AP = A_M(IJK,0)
               IF (AP /= (-ONE)) THEN
                  A_M(IJK,0) = AP*F1
                  B_M(IJK) = B_M(IJK) + AP*VAR(IJK)*F2
               ENDIF
            ENDIF
         END DO

      else if (axis.eq.'U') then
         DO IJK = ijkstart3, ijkend3
            IF (FLOW_AT_E(IJK)) THEN
               AP = A_M(IJK,0)
               IF (AP /= (-ONE)) THEN
                  A_M(IJK,0) = AP*F1
                  B_M(IJK) = B_M(IJK) + AP*VAR(IJK)*F2
               ENDIF
            ENDIF
         END DO

      else if (axis.eq.'V') then
         DO IJK = ijkstart3, ijkend3
            IF (FLOW_AT_N(IJK)) THEN
               AP = A_M(IJK,0)
               IF (AP /= (-ONE)) THEN
                  A_M(IJK,0) = AP*F1
                  B_M(IJK) = B_M(IJK) + AP*VAR(IJK)*F2
               ENDIF
            ENDIF
         END DO

      else if (axis.eq.'W') then
         DO IJK = ijkstart3, ijkend3
            IF (FLOW_AT_T(IJK)) THEN
               AP = A_M(IJK,0)
               IF (AP /= (-ONE)) THEN
                  A_M(IJK,0) = AP*F1
                  B_M(IJK) = B_M(IJK) + AP*VAR(IJK)*F2
               ENDIF
            ENDIF
         END DO

      endif

      RETURN

   contains

      include 'functions.inc'

   END SUBROUTINE UNDER_RELAX

END MODULE ur_facs
