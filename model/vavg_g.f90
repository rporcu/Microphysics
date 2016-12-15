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
      DOUBLE PRECISION FUNCTION VAVG_G (vel_g, ep_g, vol, flag)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE param1   , ONLY: ZERO

      IMPLICIT NONE

      ! Cell volume
      DOUBLE PRECISION, INTENT(IN   ) :: VOL

      ! component of gas velocity
      DOUBLE PRECISION, INTENT(IN   ) :: vel_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      integer, intent(in   ) :: flag &
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

!     Integral of U_g*EP_g for entire volume
      DOUBLE PRECISION :: sum_g

!                      Indices
      INTEGER :: I,J,K

!                      Total volume of computational cells
      DOUBLE PRECISION SUM_VOL

!  Integrate the velocity values for the whole domain,
!
      SUM_G = ZERO
      SUM_VOL = ZERO

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IF (1.eq.flag(i,j,k,1)) THEN
            SUM_VOL = SUM_VOL + VOL
            SUM_G = SUM_G + vel_G(I,J,K)*EP_G(I,J,K)*VOL
         ENDIF
      END DO
      END DO
      END DO

      ! CALL GLOBAL_ALL_SUM(SUM_VOL)
      ! CALL GLOBAL_ALL_SUM(SUM_G)

      VAVG_G = SUM_G/SUM_VOL

      RETURN
      END FUNCTION VAVG_G

      DOUBLE PRECISION FUNCTION VAVG_FLUX_G (FLUX_G, A_FACE, flag)

      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE param1, only: zero

      IMPLICIT NONE

! face area - scalar cell
      DOUBLE PRECISION, INTENT(IN) :: A_FACE

! gas mass flux
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) ::  Flux_g
      integer, intent(in   ) :: flag &
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

      INTEGER :: I,J,K

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

         IF (1.eq.flag(i,j,k,1)) THEN
            SUM_G = SUM_G + Flux_g(I,J,K)
            SUM_AREA = SUM_AREA + A_FACE
         ENDIF
      END DO
      END DO
      END DO

      ! CALL GLOBAL_ALL_SUM(SUM_AREA)
      ! CALL GLOBAL_ALL_SUM(SUM_G)

      VAVG_Flux_G = SUM_G/SUM_AREA

      RETURN
      END FUNCTION VAVG_Flux_G
END MODULE VAVG_MOD
