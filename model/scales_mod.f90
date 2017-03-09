!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: scales_mod.f                                           C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      module scales

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      ! reference pressure
      real(c_real) :: P_ref

      ! pressure scale
      real(c_real) :: P_scale

      contains

      real(c_real) function scale_pressure(XXX)
      IMPLICIT NONE
      real(c_real), intent(IN) :: XXX
      scale_pressure   = (XXX - P_ref) / P_scale
      END function scale_pressure

      real(c_real) function UNscale_pressure(XXX)
      IMPLICIT NONE
      real(c_real), intent(IN) :: XXX
      UNscale_pressure = (XXX * P_scale + P_ref)
      end function UNscale_pressure

      end module scales
