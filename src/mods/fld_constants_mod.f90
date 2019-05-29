!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: fld_constants                                               C
!  Purpose: Common block containing field variable constants           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
module fld_const

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   ! Specified constant gas density
   real(rt) :: ro_g0

   ! Specified constant gas viscosity
   real(rt) :: mu_g0

   ! Average molecular weight of gas
   real(rt) :: mw_avg

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutines: get_ro_g0, get_mu_g0, get_mw_avg                       C
!  Purpose: Getters for the quantities defined in this module          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  real(rt) function get_ro_g0() bind(C)
    get_ro_g0 = ro_g0
    return
  end function get_ro_g0

  real(rt) function get_mu_g0() bind(C)
    get_mu_g0 = mu_g0
    return
  end function get_mu_g0

  real(rt) function get_mw_avg() bind(C)
    get_mw_avg = mw_avg
    return
  end function get_mw_avg

end module fld_const
