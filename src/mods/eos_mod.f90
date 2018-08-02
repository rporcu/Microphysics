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

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function: sutherland                                                !
!                                                                      !
!  Purpose: Compute a default value of gas viscosity where gas is      !
!  assumed to be air                                                   !
!                                                                      !
!  Literature/Document References:                                     !
!     Perry, R. H., and Chilton, C. H., Chemical Engineers' Handbook,  !
!        5th Edition, McGraw-Hill Inc., 1973, pp. 248, eqn. 3-133.     !
!     Arnold, J. H., Vapor viscosities and the Sutherland equation,    !
!        Journal of Chemical Physics, 1 (2), 1933, pp. 170-176.        !
!     Sutherland, W., The Viscosity of Gases and Molecular Force,      !
!        Phil. Mag. 5:507-531, 1893.                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  real(rt) function sutherland(tg)

   implicit none

   ! Dummy arguments
   !---------------------------------------------------------------------//
   ! gas temperature
   real(rt), intent(in) :: tg

   ! Gas viscosity   (in Poise or Pa.s)
   ! Calculating gas viscosity using Sutherland's formula with
   ! Sutherland's constant (C) given by Vogel's equation C = 1.47*Tb.
   ! For air  C = 110 (Tb=74.82)
   !         mu = 1.71*10-4 poise at T = 273K

   sutherland = 1.7d-5 * (tg/273.0d0)**1.5d0 * (383.d0/(tg+110.d0))
   return
 end function sutherland

end module eos
