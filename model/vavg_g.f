MODULE VAVG_MOD
CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VAVG_g                                                 C
!  Purpose: Volume average                                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-APR-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      DOUBLE PRECISION FUNCTION VAVG_G (vel_g, vol)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
         USE fldvar, ONLY: EP_G
         USE functions, ONLY: IS_ON_myPE_wobnd, FUNIJK
         USE functions, ONLY: IS_ON_myPE_wobnd, fluid_cell
         USE mpi_utility, ONLY: GLOBAL_ALL_SUM
         USE param1, ONLY: ZERO
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

      IMPLICIT NONE

! Cell volume
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: VOL

! component of gas velocity
      DOUBLE PRECISION, INTENT(IN) ::  VEL_G(:)

!     Integral of U_g*EP_g for entire volume
      DOUBLE PRECISION :: sum_g

!                      Indices
      INTEGER :: I,J,K, IJK

!                      Total volume of computational cells
      DOUBLE PRECISION SUM_VOL

!  Integrate the velocity values for the whole domain,
!
      SUM_G = ZERO
      SUM_VOL = ZERO

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
      IF(.NOT.IS_ON_myPE_wobnd(I,J,K)) CYCLE
         IF (fluid_cell(i,j,k)) THEN
            SUM_VOL = SUM_VOL + VOL(IJK)
            SUM_G = SUM_G + vel_G(IJK)*EP_G(IJK)*VOL(IJK)
         ENDIF
      END DO
      END DO
      END DO

      CALL GLOBAL_ALL_SUM(SUM_VOL)
      CALL GLOBAL_ALL_SUM(SUM_G)

      VAVG_G = SUM_G/SUM_VOL

      RETURN
      END FUNCTION VAVG_G

      DOUBLE PRECISION FUNCTION VAVG_FLUX_G (FLUX_G, A_FACE)

      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE functions, ONLY: is_on_mype_wobnd, fluid_cell
      USE functions, ONLY: funijk
      USE mpi_utility, ONLY: global_all_sum
      USE param1

      IMPLICIT NONE

! face area - scalar cell
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: A_FACE

! gas mass flux
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) ::  Flux_g

      INTEGER :: I,J,K, IJK

! Integral of U_g*ROP_g*Area
      DOUBLE PRECISION :: SUM_G

!  Total volume of computational cells
      DOUBLE PRECISION :: SUM_AREA

!  Integrate the velocity values for the whole domain,

      SUM_G = ZERO
      SUM_AREA = ZERO

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
         IF(.NOT.IS_ON_myPE_wobnd(I,J,K)) CYCLE
         IF (fluid_cell(i,j,k)) THEN
            SUM_G = SUM_G + Flux_g(IJK)
            SUM_AREA = SUM_AREA + A_FACE(IJK)
         ENDIF
      END DO
      END DO
      END DO

      CALL GLOBAL_ALL_SUM(SUM_AREA)
      CALL GLOBAL_ALL_SUM(SUM_G)

      VAVG_Flux_G = SUM_G/SUM_AREA

      RETURN
      END FUNCTION VAVG_Flux_G
END MODULE VAVG_MOD
