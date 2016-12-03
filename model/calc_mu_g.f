!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MU_g                                               C
!  Author: W. Sams/M. Syamlal                         Date: 18-JUL-94  C
!                                                                      C
!  Purpose: Calculate the molecular viscosity.                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MU_G(lambda_g,mu_g,mu_g0)

      USE param1   , only: zero, undefined
      USE functions, only: fluid_at
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3

      IMPLICIT NONE

      double precision, intent(inout) :: lambda_g (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(inout) ::     mu_g (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in   ) ::     mu_g0
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Cell indices
      INTEGER :: I,J,K

      DO K = kstart3, kend3
        DO J = jstart3, jend3
           DO I = istart3, iend3

! Gas viscosity   (in Poise or Pa.s)
! Calculating gas viscosity using Sutherland's formula with
! Sutherland's constant (C) given by Vogel's equation C = 1.47*Tb.
! For air  C = 110 (Tb=74.82)
!         mu = 1.71*10-4 poise at T = 273K

         IF (fluid_at(i,j,k)) THEN

            IF (mu_g0 == UNDEFINED) MU_G(I,J,K) = 1.7D-5 * &
               (293.15d0/273.0D0)**1.5D0 * (383.D0/(293.15d0+110.D0))

            LAMBDA_G(i,j,k) = -F2O3*MU_G(I,J,K)

         ELSE

            MU_G(I,J,K)  = ZERO
            LAMBDA_G(i,j,k) = ZERO

         ENDIF   ! end if (fluid_at(i,j,k))

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_MU_G
