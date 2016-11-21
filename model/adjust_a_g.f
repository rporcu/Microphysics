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
   USE geometry, only: ayz_u, axz_v, axy_w
   USE fldvar, only: rop_g
   USE functions, only: iminus, jminus, kminus, ieast, jnorth, ktop, funijk

   contains

      double precision function denom_u_neg(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_u_neg = ROP_G(FUNIJK(ieast(i,j,k),j,k))*AYZ_U(FUNIJK(i,j,k))
      end function denom_u_neg

      double precision function denom_u_pos(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_u_pos = ROP_G(FUNIJK(i,j,k))*AYZ_U(FUNIJK(iminus(i,j,k),j,k))
      end function denom_u_pos

      double precision function denom_v_neg(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_v_neg = ROP_G(FUNIJK(i,jnorth(i,j,k),k))*AXZ_V(FUNIJK(i,j,k))
      end function denom_v_neg

      double precision function denom_v_pos(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_v_pos = ROP_G(FUNIJK(i,j,k))*AXZ_V(FUNIJK(i,jminus(i,j,k),k))
      end function denom_v_pos

      double precision function denom_w_neg(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_w_neg = ROP_G(FUNIJK(i,j,ktop(i,j,k)))*AXY_W(FUNIJK(i,j,k))
      end function denom_w_neg

      double precision function denom_w_pos(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         denom_w_pos = ROP_G(FUNIJK(i,j,k))*AXY_W(FUNIJK(i,j,kminus(i,j,k)))
      end function denom_w_pos

      SUBROUTINE ADJUST_A_G(axis, A_M, B_M)

      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3, imap
      USE compar, only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar, only: istart1, iend1, jstart1, jend1, kstart1, kend1
         USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
         USE functions, only: ip1
         USE matrix, only: e, w, s, n, t, b
         USE param, only: dimension_3
         USE param1, only: ONE, ZERO, small_number
      use functions, only: funijk
      IMPLICIT NONE

!                      Indices
      INTEGER          IP, IJK
!
!                      Phase index
      INTEGER          M, I, J, K
!
      CHARACTER, INTENT(IN) :: axis
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)

      DOUBLE PRECISION :: denominator, xxxm, xxxp

      abstract interface
         function denom (i,j,k)
            DOUBLE PRECISION :: denom
            integer, intent (in) :: i,j,k
         end function denom
      end interface

      abstract interface
         function avg_iface(XXXm, XXXp, xL)
            DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
            INTEGER, INTENT(IN) :: xL
            DOUBLE PRECISION :: AVG_IFACE
         end function avg_iface
      end interface

      procedure (avg_iface), pointer :: avg => null ()
      procedure (denom), pointer :: denom_neg => null ()
      procedure (denom), pointer :: denom_pos => null ()
!-----------------------------------------------

      M = 0

      if (axis.eq.'U') then
         avg => avg_x_e
         denom_neg => denom_u_neg
         denom_pos => denom_u_pos
      else if (axis.eq.'V') then
         avg => avg_y_n
         denom_neg => denom_w_neg
         denom_pos => denom_w_pos
      else if (axis.eq.'W') then
         avg => avg_z_t
         denom_neg => denom_w_neg
         denom_pos => denom_w_pos
      endif

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         IJK = FUNIJK(i,j,k)

         IF (ABS(A_M(IJK,0)) < SMALL_NUMBER) THEN
            A_M(IJK,E) = ZERO
            A_M(IJK,W) = ZERO
            A_M(IJK,N) = ZERO
            A_M(IJK,S) = ZERO
            A_M(IJK,T) = ZERO
            A_M(IJK,B) = ZERO
            A_M(IJK,0) = -ONE
            IF (B_M(IJK) < ZERO) THEN
               IP = IP1(I)

               denominator = denom_neg(i,j,k)
               xxxm = ONE
               xxxp = ZERO
            ELSE IF (B_M(IJK) > ZERO) THEN
               IP = I
               denominator = denom_pos(i,j,k)
               xxxm = ZERO
               xxxp = ONE
            ELSE
               denominator = ZERO
            ENDIF

            IF (denominator > SMALL_NUMBER) THEN
               B_M(IJK) = SQRT(ABS(B_M(IJK))/(denominator*AVG(xxxm,xxxp,IP)))
            ELSE
               B_M(IJK) = ZERO
            ENDIF
         ENDIF
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE ADJUST_A_G
   end MODULE ADJUST_A
