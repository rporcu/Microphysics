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
      use functions, only: fluid_at
      use param1   , only: zero, half, one, undefined

      implicit none

      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: lambda_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer, intent(in   ) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,0:4)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: i,j,k
!--------------------

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

           IF (flag(i,j,k,1)<100) THEN
! Fluid and inflow/outflow cells: FLAG < 100
            IF (RO_G0 /= UNDEFINED) RO_G(I,J,K) = RO_G0
           ENDIF

           IF (fluid_at(i,j,k)) THEN
! Strictly Fluid cells: FLAG = 1
            IF (MU_G0 /= UNDEFINED) THEN
               MU_G(i,j,k) = MU_G0
               LAMBDA_G(i,j,k) = -(2.0d0/3.0d0)*MU_G0
              ENDIF
           ENDIF

          end do
        end do
      end do

      END SUBROUTINE SET_CONSTPROP
