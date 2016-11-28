!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MU_g                                               C
!  Author: W. Sams/M. Syamlal                         Date: 18-JUL-94  C
!                                                                      C
!  Purpose: Calculate the molecular viscosity.                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MU_G()

      USE param
      USE param1
      USE physprop
      USE geometry
      USE fldvar
      USE constant
      USE toleranc
      USE compar
      USE functions
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

      IMPLICIT NONE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Cell indices
      INTEGER :: I,J,K, IJK

      DO K = kstart3, kend3
        DO J = jstart3, jend3
           DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)

! Gas viscosity   (in Poise or Pa.s)
! Calculating gas viscosity using Sutherland's formula with
! Sutherland's constant (C) given by Vogel's equation C = 1.47*Tb.
! For air  C = 110 (Tb=74.82)
!         mu = 1.71*10-4 poise at T = 273K
         IF (fluid_at(i,j,k)) THEN

            IF (MU_G0 == UNDEFINED) MU_G(IJK) = to_SI*1.7D-4 * &
               (293.15d0/273.0D0)**1.5D0 * (383.D0/(293.15d0+110.D0))

            LAMBDA_G(i,j,k) = -F2O3*MU_G(IJK)

         ELSE
            MU_G(IJK)  = ZERO
            LAMBDA_G(i,j,k) = ZERO
         ENDIF   ! end if (fluid_at(i,j,k))

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_MU_G
