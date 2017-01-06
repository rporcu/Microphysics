   MODULE ADJUST_A

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_A_U_g(A_m, B_m, IER)                            C
!  Author: M. Syamlal                                 Date:  2-AUG-96  C
!                                                                      C
!  Purpose: Handle the special case of the center coefficient in       C
!  U_g momentum eq. becoming zero.                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

   contains

      subroutine adjust_a_g(axis, slo, shi, lo, hi, A_M, B_M, ROP_G, dx, dy, dz)

         USE functions, only: avg
         USE functions, only: ip1
         USE matrix, only: e, w, s, n, t, b
         USE param1, only: ONE, ZERO, small_number

         USE functions, only: iminus, jminus, kminus, ieast, jnorth, ktop

         use compar, only: istart2,iend2,jstart2,jend2,kstart2,kend2

         implicit none

         integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
!---------------------------------------------------------------------//
      CHARACTER, INTENT(IN) :: axis

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3), -3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), INTENT(IN   ) :: ROP_G&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), INTENT(IN   ) :: dx, dy, dz
!---------------------------------------------------------------------//

      INTEGER          IP
      INTEGER          I, J, K

      real(c_real) :: denominator, xxxm, xxxp
      real(c_real) :: axy, ayz, axz

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

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

               if (axis .eq. 'U') then
                  denominator = ROP_G(ieast(i,j,k),j,k)*AYZ
               else if (axis .eq. 'V') then
                  denominator = ROP_G(i,jnorth(i,j,k),k)*AXZ
               else if (axis .eq. 'W') then
                  denominator = ROP_G(i,j,ktop(i,j,k))*AXY
               end if

               xxxm = ONE
               xxxp = ZERO

            ELSE IF (B_M(I,J,K) > ZERO) THEN

               if (axis .eq. 'U') then
                  denominator = ROP_G(i,j,k)*AYZ
               else if (axis .eq. 'V') then
                  denominator = ROP_G(i,j,k)*AXZ
               else if (axis .eq. 'W') then
                  denominator = ROP_G(i,j,k)*AXY
               end if

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

      end subroutine adjust_a_g

   end MODULE ADJUST_A
