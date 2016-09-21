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
   USE functions, only: i_of, im_of, jm_of, km_of, east_of, north_of, top_of

   contains

      double precision function denom_u_neg(ijk)
         implicit none
         integer, intent(in) :: ijk
         denom_u_neg = ROP_G(EAST_OF(IJK))*AYZ_U(IJK)
      end function denom_u_neg

      double precision function denom_u_pos(ijk)
         implicit none
         integer, intent(in) :: ijk
         denom_u_pos = ROP_G(ijk)*AYZ_U(IM_OF(IJK))
      end function denom_u_pos

      double precision function denom_v_neg(ijk)
         implicit none
         integer, intent(in) :: ijk
         denom_v_neg = ROP_G(NORTH_OF(IJK))*AXZ_V(IJK)
      end function denom_v_neg

      double precision function denom_v_pos(ijk)
         implicit none
         integer, intent(in) :: ijk
         denom_v_pos = ROP_G(ijk)*AXZ_V(JM_OF(IJK))
      end function denom_v_pos

      double precision function denom_w_neg(ijk)
         implicit none
         integer, intent(in) :: ijk
         denom_w_neg = ROP_G(TOP_OF(IJK))*AXY_W(IJK)
      end function denom_w_neg

      double precision function denom_w_pos(ijk)
         implicit none
         integer, intent(in) :: ijk
         denom_w_pos = ROP_G(ijk)*AXY_W(KM_OF(IJK))
      end function denom_w_pos

      SUBROUTINE ADJUST_A_G(axis, A_M, B_M)

         USE compar, only: ijkstart3, ijkend3
         USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
         USE indices, only: ip1
         USE matrix, only: e, w, s, n, t, b
         USE param, only: dimension_3, dimension_m
         USE param1, only: ONE, ZERO, small_number

      IMPLICIT NONE

!                      Indices
      INTEGER          IP, IJK
!
!                      Phase index
      INTEGER          M
!
      CHARACTER, INTENT(IN) :: axis
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)

      DOUBLE PRECISION :: denominator, xxxm, xxxp

      abstract interface
         function denom (ijk)
            DOUBLE PRECISION :: denom
            integer, intent (in) :: ijk
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

      DO IJK = ijkstart3, ijkend3
         IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN
            A_M(IJK,E,M) = ZERO
            A_M(IJK,W,M) = ZERO
            A_M(IJK,N,M) = ZERO
            A_M(IJK,S,M) = ZERO
            A_M(IJK,T,M) = ZERO
            A_M(IJK,B,M) = ZERO
            A_M(IJK,0,M) = -ONE
            IF (B_M(IJK) < ZERO) THEN
               IP = IP1(I_OF(IJK))

               denominator = denom_neg(ijk)
               xxxm = ONE
               xxxp = ZERO
            ELSE IF (B_M(IJK) > ZERO) THEN
               IP = I_OF(IJK)
               denominator = denom_pos(ijk)
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

      RETURN
      END SUBROUTINE ADJUST_A_G
   end MODULE ADJUST_A
