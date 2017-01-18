MODULE VAVG_MOD

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

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
      real(c_real) function vavg_g (slo, shi, vel_g, vlo, vhi, ep_g, vol, flag)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param1   , ONLY: ZERO

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),vlo(3),vhi(3)

      ! Cell volume
      real(c_real), INTENT(IN   ) :: VOL

      ! component of gas velocity
      real(c_real), INTENT(IN   ) :: vel_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), INTENT(IN   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer, intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

!     Integral of U_g*EP_g for entire volume
      real(c_real) :: sum_g

!                      Indices
      INTEGER :: I,J,K

!                      Total volume of computational cells
      real(c_real) SUM_VOL

!  Integrate the velocity values for the whole domain,
!
      SUM_G = ZERO
      SUM_VOL = ZERO

      DO K = slo(3),shi(3)
      DO J = slo(2),shi(2)
      DO I = slo(1),shi(1)
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

      end function vavg_g

      real(c_real) function vavg_flux_g (slo, shi, FLUX_G, A_FACE, flag)

      USE param1, only: zero

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      ! face area - scalar cell
      real(c_real), INTENT(IN) :: A_FACE

      ! gas mass flux
      real(c_real), intent(in) ::  Flux_g(:,:,:)

      integer    , intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      INTEGER :: I,J,K

      ! Integral of U_g*ROP_g*Area
      real(c_real) :: SUM_G

      !  Total volume of computational cells
      real(c_real) :: SUM_AREA

      !  Integrate the velocity values for the whole domain,

      SUM_G = ZERO
      SUM_AREA = ZERO

        DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
        DO I = slo(1),shi(1)
         IF (1.eq.flag(i,j,k,1)) THEN
            SUM_G = SUM_G + Flux_g(I,J,K)
            SUM_AREA = SUM_AREA + A_FACE
         ENDIF
      END DO
      END DO
      END DO

      ! CALL GLOBAL_ALL_SUM(SUM_AREA)
      ! CALL GLOBAL_ALL_SUM(SUM_G)

      vavg_flux_g = SUM_G/SUM_AREA

      end function vavg_flux_g
END MODULE VAVG_MOD
