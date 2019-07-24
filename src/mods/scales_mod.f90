!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: scales_mod.f                                           C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      module scales

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      implicit none

      ! reference pressure
      real(rt) :: P_ref

      ! pressure scale
      real(rt) :: P_scale

      contains

      real(rt) function get_p_ref() bind(C)
      IMPLICIT NONE
      get_p_ref = p_ref
      END function get_p_ref

      real(rt) function get_p_scale() bind(C)
      IMPLICIT NONE
      get_p_scale = p_scale
      END function get_p_scale

      real(rt) function scale_pressure(XXX)
      IMPLICIT NONE
      real(rt), intent(IN) :: XXX
      scale_pressure   = (XXX - P_ref) / P_scale
      END function scale_pressure

      real(rt) function UNscale_pressure(XXX)
      IMPLICIT NONE
      real(rt), intent(IN) :: XXX
      UNscale_pressure = (XXX * P_scale + P_ref)
      end function UNscale_pressure

      end module scales
