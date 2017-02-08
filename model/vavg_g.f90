module vavg_mod

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   contains
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
      real(c_real) function vavg_flux_g (dimlo, dimhi, flux_g, A_FACE)

      USE param1, only: zero

      IMPLICIT NONE

      integer     , intent(in   ) :: dimlo(3),dimhi(3)

      ! face area - scalar cell
      real(c_real), intent(IN) :: A_FACE

      ! gas mass flux
      real(c_real), intent(in) ::  flux_g(dimlo(1):dimhi(1),dimlo(2):dimhi(2),dimlo(3):dimhi(3))

      INTEGER :: I,J,K

      ! Integral of U_g*ROP_g*Area
      real(c_real) :: SUM_G

      !  Total volume of computational cells
      real(c_real) :: SUM_AREA

      !  Integrate the velocity values for the whole domain,

      sum_g = zero
      sum_area = zero

      DO K = dimlo(3),dimhi(3)
      DO J = dimlo(2),dimhi(2)
      DO I = dimlo(1),dimhi(1)
            sum_g = sum_g + flux_g(I,J,K)
            sum_area = sum_area + A_FACE
      END DO
      END DO
      END DO

      vavg_flux_g = sum_g/sum_area

      end function vavg_flux_g

end module vavg_mod
