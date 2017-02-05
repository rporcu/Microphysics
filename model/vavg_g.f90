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

      vavg_flux_g = SUM_G/SUM_AREA

      end function vavg_flux_g
END MODULE VAVG_MOD
