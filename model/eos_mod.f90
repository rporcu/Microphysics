!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EOS                                                    C
!  Purpose: Equation of state for gas and initial solids density       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE eos

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: EOSG                                                      C
!  Purpose: Equation of state for gas                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      real(c_real) FUNCTION EOSG (MW, PG, TG)

! Global Variables:
!---------------------------------------------------------------------//
      USE constant, only: gas_const
      USE scales, only: unscale_pressure
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      real(c_real), INTENT(IN) :: MW, PG, TG

      EOSG = UNSCALE_PRESSURE(PG)*MW/(GAS_CONST*TG)
      RETURN
      END FUNCTION EOSG


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: dROodP_g                                                  C
!  Purpose: derivative of gas density w.r.t pressure                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-AUG-96  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      real(c_real) FUNCTION DROODP_G (ROG, PG)

! Global Variables:
!---------------------------------------------------------------------//
      USE scales, only: p_ref
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! gas density and pressure
      real(c_real), INTENT(IN) :: ROG, PG

      DROODP_G = ROG/(PG + P_REF)
      RETURN
      END FUNCTION DROODP_G



      END MODULE eos
