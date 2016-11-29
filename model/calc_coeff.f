!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      !
!  Subroutine: CALC_COEFF_ALL                                          !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, exchange rates, and reaction rates.  !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_COEFF_ALL(FLAG, IER)

! Global variables:
!-----------------------------------------------------------------------
      ! Flag for explcit coupling between the fluid and particles.
      use discretelement, only: DES_EXPLICITLY_COUPLED

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! FLAG = 0, overwrite the coeff arrays, (e.g. start of a time step)
! FLAG = 1, do not overwrite
      INTEGER, intent(in) :: FLAG
! Error index
      INTEGER, intent(inout) :: IER
!-----------------------------------------------

      ! Calculate all physical properties, transport properties, 
      ! and exchange rates.
      CALL CALC_COEFF(IER, 2)

      IF (DES_EXPLICITLY_COUPLED) CALL CALC_DRAG_DES_EXPLICIT

      END SUBROUTINE CALC_COEFF_ALL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_COEFF                                              !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!  Reviewer:                                          Date:            !
!                                                                      !
!                                                                      !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_COEFF(IER, pLevel)

      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DES_CONTINUUM_COUPLED

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Error index
      INTEGER, intent(inout) :: IER
! Level to calculate physical properties.
! 0) Only density
! 1) Everything but density
! 2) All physical properties
      INTEGER, intent(in) :: pLevel
!-----------------------------------------------------------------------

! Calculate physical properties: (density, specific heat, diameter)
      CALL PHYSICAL_PROP(IER, pLevel)

! Calculate interphase coeffs: (momentum and energy)
      IF (DES_CONTINUUM_COUPLED .AND. .NOT.DES_EXPLICITLY_COUPLED) &
         CALL CALC_DRAG_DES_2FLUID

      END SUBROUTINE CALC_COEFF

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_TRD_AND_TAU                                        !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!  Reviewer:                                          Date:            !
!                                                                      !
!                                                                      !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_TRD_AND_TAU()

! Stress tensor cross terms.
      USE fldvar, only : TAU_U_G, TAU_V_G, TAU_W_G

      implicit none

!-----------------------------------------------------------------------

! Calculate the trace of the stress tensor (gas phase; m=0)
      CALL CALC_TRD_G

! Calculate the cross terms of the stress tensor (gas phase; m=0)
      CALL CALC_TAU_U_G (TAU_U_G)
      CALL CALC_TAU_V_G (TAU_V_G)
      CALL CALC_TAU_W_G (TAU_W_G)


      RETURN
      END SUBROUTINE CALC_TRD_AND_TAU
