!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_PG_GRAD                                            !
!  Purpose: Calculate cell centered pressure force exerted on the      !
!           particles in the cell by the gas/fluid phase               !
!           Note that P_force is evaluated as -dp/dx                   !
!                                                                      !
!  Notes: This pressure force only needs to be calculated once during  !
!         the DEM loop (at the beginning) since the gas/fluid phase    !
!         is essentially static at that point (i.e., gas field is not  !
!         updated during DEM loop                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PG_GRAD(p_g)

      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3

! Particle volume.
      use discretelement, only: PVOL
! Gas pressure force by fluid cell
      use discretelement, only: P_FORCE
! Particle drag force
      use discretelement, only: DRAG_FC
! Flag for 3D simulatoins.
      use geometry, only: DO_K
! Loop bounds for fluid grid
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flags for cyclic BC with pressure drop
      use geometry, only: CYCLIC_X_PD, CYCLIC_Y_PD, CYCLIC_Z_PD
! Specified pressure drop
      use bc, only: DELP_X, DELP_Y, DELP_Z
! Domain length
      use geometry, only: XLENGTH, YLENGTH, ZLENGTH
! Gas phase pressure

      use discretelement, only: MAX_PIP, PIJK, DES_EXPLICITLY_COUPLED

      use functions, only: funijk
      use functions, only: fluid_at

      use functions, only: IS_NONEXISTENT
      use functions, only: IS_ENTERING, IS_ENTERING_GHOST
      use functions, only: IS_EXITING, IS_EXITING_GHOST

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO

      implicit none

      DOUBLE PRECISION, INTENT(INOUT) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, I, J, K
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lPF(3)
! Loop bound for
      INTEGER :: LP_BND
! mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION :: cPG(3)
!......................................................................!

! Calculate the gas phase pressure gradient. (dP/dx)
      CALL CALC_GRAD_DES(P_G, P_FORCE)

! Add in cyclic BC pressure drop.
      cPG(1) = merge(DELP_X/XLENGTH, ZERO, CYCLIC_X_PD)
      cPG(2) = merge(DELP_Y/YLENGTH, ZERO, CYCLIC_Y_PD)
      cPG(3) = merge(DELP_Z/ZLENGTH, ZERO, CYCLIC_Z_PD)

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
         P_FORCE(:,IJK) = cPG - P_FORCE(:,IJK)

      ENDDO
      ENDDO
      ENDDO

      IF(DES_EXPLICITLY_COUPLED) THEN

! Loop bounds for interpolation.
         LP_BND = merge(27,9,DO_K)

! Calculate the gas phase forces acting on each particle.

         DO NP=1,MAX_PIP

            IF(IS_NONEXISTENT(NP) .or.                              &
               IS_ENTERING(NP) .or. IS_ENTERING_GHOST(NP) .or.      &
               IS_EXITING(NP)  .or. IS_EXITING_GHOST(NP)) CYCLE

            i = PIJK(NP,1)
            j = PIJK(NP,2)
            k = PIJK(NP,3)
            if (.NOT.fluid_at(i,j,k)) CYCLE

             lPF = P_FORCE(:,PIJK(NP,4))

! Include gas pressure and gas-solids drag
            DRAG_FC(NP,:) = DRAG_FC(NP,:) + lPF*PVOL(NP)

         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CALC_PG_GRAD
