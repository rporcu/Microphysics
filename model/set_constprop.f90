MODULE SET_CONSTPROP_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_constprop                                           C
!  Purpose: This routine serves two purposes:                          C
!    1) initializes various variables everywhere in the domain with    C
!       a zero value. the real need for this is unclear. undefined     C
!       may be a better approach...                                    C
!    2) if defined, sets physical properties to their specified        C
!       constant value in the fluid domain. cannot set in flow         C
!       boundaries or later checks will also complain (may be          C
!       overly strict check)                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_CONSTPROP(ro_g,lambda_g,mu_g,flag)

! Modules
!-----------------------------------------------

      use compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use fld_const, only: ro_g0, mu_g0
      use param1   , only: zero, is_defined

      implicit none

      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: lambda_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer, intent(in   ) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

! First, initialize certain transport coefficients, physical
! properties, and other variable types everywhere in the
! domain to zero. Then, set these to their specified value
! (if defined) in fluid cells. Some are also set in flow cells.
! NOTE: DO NOT simply zero existing field variables.

      RO_g = ZERO
      MU_g = ZERO
      LAMBDA_G = ZERO

! Set specified constant physical properties values

      if (IS_DEFINED(ro_g0)) then
         ! Fluid and inflow/outflow cells: FLAG < 100
         where (flag < 100) ro_g = ro_g0
      end if

      if (IS_DEFINED(mu_g0)) then
         ! Strictly Fluid cells: FLAG = 1
         where (flag .eq. 1) mu_g = mu_g0
         where (flag .eq. 1) lambda_g = -(2.0d0/3.0d0)*mu_g0
      end if

      END SUBROUTINE SET_CONSTPROP
END MODULE SET_CONSTPROP_MODULE
