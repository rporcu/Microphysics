!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EOS                                                    C
!  Purpose: Equation of state for gas and initial solids density       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
module eos

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: EOSG                                                      C
!  Purpose: Equation of state for gas                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  real(rt) function eosg (mw, pg, tg)

! Global Variables:
!---------------------------------------------------------------------//
    use constant, only: gas_const
    use scales, only: unscale_pressure

    implicit none

! Dummy arguments
!---------------------------------------------------------------------//
    real(rt), intent(in) :: mw, pg, tg

    eosg = unscale_pressure(pg)*mw/(gas_const*tg)
    return
  end function eosg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: dROodP_g                                                  C
!  Purpose: derivative of gas density w.r.t pressure                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-AUG-96  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  real(rt) function droodp_g (rog, pg)

! Global Variables:
!---------------------------------------------------------------------//
    use scales, only: p_ref
    implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! gas density and pressure
    real(rt), intent(in) :: rog, pg

    droodp_g = rog/(pg + p_ref)
    return
  end function droodp_g

end module eos
