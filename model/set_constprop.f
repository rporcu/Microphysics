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
      SUBROUTINE SET_CONSTPROP

! Modules
!-----------------------------------------------
      USE param1, only: zero, half, one, undefined

      USE fldvar, only: ro_g

      USE fldvar, only: LAMBDA_G

      USE fldvar, only: ro_g0, mu_g0, mu_g

      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use functions, only: wall_at, fluid_cell, funijk

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: i,j,k,ijk
!-----------------------------------------------

! First, initialize certain transport coefficients, physical
! properties, and other variable types everywhere in the
! domain to zero. Then, set these to their specified value
! (if defined) in fluid cells. Some are also set in flow cells.
! NOTE: DO NOT simply zero existing field variables.
      RO_g = ZERO
      MU_g = ZERO
      LAMBDA_G = ZERO

! Set specified constant physical properties values
      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3

           ijk = funijk(i,j,k)

           IF (.NOT.WALL_AT(IJK)) THEN
! Fluid and inflow/outflow cells: FLAG < 100
            IF (RO_G0 /= UNDEFINED) RO_G(IJK) = RO_G0
           ENDIF

           IF (fluid_cell(i,j,k)) THEN
! Strictly Fluid cells: FLAG = 1
            IF (MU_G0 /= UNDEFINED) THEN
               MU_G(IJK) = MU_G0
               LAMBDA_G(IJK) = -(2.0d0/3.0d0)*MU_G0
              ENDIF
           ENDIF

          end do
        end do
      end do

      END SUBROUTINE SET_CONSTPROP
