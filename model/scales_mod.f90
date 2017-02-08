!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: scales_mod.f                                           C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      MODULE scales

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

! reference pressure
      real(c_real) :: P_ref

! pressure scale
      real(c_real) :: P_scale

      CONTAINS

      real(c_real) FUNCTION SCALE_PRESSURE(XXX)
      IMPLICIT NONE
      real(c_real), intent(IN) :: XXX
      SCALE_PRESSURE   = (XXX - P_ref) / P_scale
      END FUNCTION SCALE_PRESSURE

      real(c_real) FUNCTION UNSCALE_PRESSURE(XXX)
      IMPLICIT NONE
      real(c_real), intent(IN) :: XXX
      UNSCALE_PRESSURE = (XXX * P_scale + P_ref)
      END FUNCTION UNSCALE_PRESSURE

      END MODULE scales
