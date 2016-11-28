   MODULE ADJUST_A
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_A_U_g(A_m, B_m, IER)                            C
!  Author: M. Syamlal                                 Date:  2-AUG-96  C
!                                                                      C
!  Purpose: Handle the special case of the center coefficient in       C
!  U_g momentum eq. becoming zero.                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   USE geometry, only: ayz, axz, axy
   USE fldvar, only: rop_g
   USE functions, only: iminus, jminus, kminus, ieast, jnorth, ktop, funijk

   contains

      double precision function denom_u_neg(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_u_neg = ROP_G(FUNIJK(ieast(i,j,k),j,k))*AYZ
      end function denom_u_neg

      double precision function denom_u_pos(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_u_pos = ROP_G(FUNIJK(i,j,k))*AYZ
      end function denom_u_pos

      double precision function denom_v_neg(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_v_neg = ROP_G(FUNIJK(i,jnorth(i,j,k),k))*AXZ
      end function denom_v_neg

      double precision function denom_v_pos(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_v_pos = ROP_G(FUNIJK(i,j,k))*AXZ
      end function denom_v_pos

      double precision function denom_w_neg(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_w_neg = ROP_G(FUNIJK(i,j,ktop(i,j,k)))*AXY
      end function denom_w_neg

      double precision function denom_w_pos(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_w_pos = ROP_G(FUNIJK(i,j,k))*AXY
      end function denom_w_pos

      SUBROUTINE ADJUST_A_G(axis, A_M, B_M)

         USE fun_avg, only: avg
         USE functions, only: ip1
         USE matrix, only: e, w, s, n, t, b
         USE param1, only: ONE, ZERO, small_number
      use functions, only: funijk
      use compar, only: istart2, iend2
      use compar, only: jstart2, jend2
      use compar, only: kstart2, kend2
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

!                      Indices
      INTEGER          IP, IJK
!
!                      Phase index
      INTEGER          M, I, J, K
!
      CHARACTER, INTENT(IN) :: axis

      DOUBLE PRECISION :: denominator, xxxm, xxxp

      abstract interface
         function denom (i,j,k)
            DOUBLE PRECISION :: denom
            integer, intent (in) :: i,j,k
         end function denom
      end interface

      procedure (denom), pointer :: denom_neg => null ()
      procedure (denom), pointer :: denom_pos => null ()
!-----------------------------------------------

      M = 0

      if (axis.eq.'U') then
         denom_neg => denom_u_neg
         denom_pos => denom_u_pos
      else if (axis.eq.'V') then
         denom_neg => denom_w_neg
         denom_pos => denom_w_pos
      else if (axis.eq.'W') then
         denom_neg => denom_w_neg
         denom_pos => denom_w_pos
      endif

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         IJK = FUNIJK(i,j,k)

         IF (ABS(A_M(I,J,K,0)) < SMALL_NUMBER) THEN
            A_M(I,J,K,E) = ZERO
            A_M(I,J,K,W) = ZERO
            A_M(I,J,K,N) = ZERO
            A_M(I,J,K,S) = ZERO
            A_M(I,J,K,T) = ZERO
            A_M(I,J,K,B) = ZERO
            A_M(I,J,K,0) = -ONE
            IF (B_M(I,J,K) < ZERO) THEN
               IP = IP1(I)

               denominator = denom_neg(i,j,k)
               xxxm = ONE
               xxxp = ZERO
            ELSE IF (B_M(I,J,K) > ZERO) THEN
               denominator = denom_pos(i,j,k)
               xxxm = ZERO
               xxxp = ONE
            ELSE
               denominator = ZERO
            ENDIF

            IF (denominator > SMALL_NUMBER) THEN
               B_M(I,J,K) = SQRT(ABS(B_M(I,J,K))/(denominator*AVG(xxxm,xxxp)))
            ELSE
               B_M(I,J,K) = ZERO
            ENDIF
         ENDIF
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE ADJUST_A_G
   end MODULE ADJUST_A
