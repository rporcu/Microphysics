!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DRAG_DES                                           !
!                                                                      !
!  Purpose: This subroutine is called from DES routines. It calls      !
!  functions that calcultate the drag force acting on particles. No    !
!  field variables are updated.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES

      use discretelement, only: DES_CONTINUUM_COUPLED

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_GARG

      use discretelement, only: DES_EXPLICITLY_COUPLED

      use discretelement, only: DRAG_FC, FC, MAX_PIP
      use functions, only: IS_NORMAL

      IMPLICIT NONE

      INTEGER :: II

! Apply the drag force calculated by the gas phase.
      IF(DES_EXPLICITLY_COUPLED) THEN

         IF(DES_CONTINUUM_COUPLED) THEN
            DO II = 1, MAX_PIP
               IF(IS_NORMAL(II)) &
                  FC(II,:) = FC(II,:) + DRAG_FC(II,:)
            ENDDO
         ENDIF

      ELSE

! Calculate gas-solids drag force on particle
         IF(DES_CONTINUUM_COUPLED) THEN
            SELECT CASE(DES_INTERP_SCHEME_ENUM)
            CASE(DES_INTERP_GARG) ; CALL DRAG_GS_DES0
            CASE DEFAULT; CALL DRAG_GS_DES1
            END SELECT
         ENDIF

      ENDIF

      RETURN
      END SUBROUTINE CALC_DRAG_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_DRAG_GS                                             !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side. It !
!  calls the correct routine for calculating the gas drag force.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES_2FLUID

      use discretelement, only: DES_CONTINUUM_COUPLED

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_GARG

      IMPLICIT NONE


! Calculate gas-solids drag force.
      IF(DES_CONTINUUM_COUPLED) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_GARG) ; CALL DRAG_GS_GAS0
         CASE DEFAULT; CALL DRAG_GS_GAS1
         END SELECT
      ENDIF

      RETURN
      END SUBROUTINE CALC_DRAG_DES_2FLUID


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DRAG_DES_EXPLICIT                                  !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side.    !
!  Moreover, it is only called once per time step so that the total    !
!  drag force is calculated explicitly for the fluid phase and DEM     !
!  particles. This is done to ensure momentum conservation and reduce  !
!  computational overhead.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES_EXPLICIT

      use discretelement, only: DES_CONTINUUM_COUPLED

      IMPLICIT NONE

! Bin particles to the fluid grid.
      CALL PARTICLES_IN_CELL
! Calculate interpolation weights
      CALL CALC_INTERP_WEIGHTS
! Calculate mean fields (EPg).
      CALL COMP_MEAN_FIELDS

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) CALL DRAG_GS_GAS1

      RETURN
      END SUBROUTINE CALC_DRAG_DES_EXPLICIT
