!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_PG_GRAD                                            !
!  Purpose: Calculate cell centered pressure force exerted on the      !
!           particles in the cell by the gas/fluid phase               !
!           Note that gradPg is evaluated as -dp/dx                    !
!                                                                      !
!  Notes: This pressure force only needs to be calculated once during  !
!         the DEM loop (at the beginning) since the gas/fluid phase    !
!         is essentially static at that point (i.e., gas field is not  !
!         updated during DEM loop                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PG_GRAD(p_g, gradPg)

      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      use discretelement, only: entering_particle, exiting_particle, entering_ghost, exiting_ghost, particle_state
      use discretelement, only: nonexistent

! Particle volume.
      use discretelement, only: PVOL
! Particle drag force
      use discretelement, only: DRAG_FC
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

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO

      implicit none

      DOUBLE PRECISION, INTENT(in   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(out  ) :: gradPg&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, I, J, K
! mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION :: cPG(3)
!......................................................................!

! Calculate the gas phase pressure gradient. (dP/dx)
      CALL CALC_GRAD_DES(P_G, gradPg)

! Add in cyclic BC pressure drop.
      cPG(1) = merge(DELP_X/XLENGTH, ZERO, CYCLIC_X_PD)
      cPG(2) = merge(DELP_Y/YLENGTH, ZERO, CYCLIC_Y_PD)
      cPG(3) = merge(DELP_Z/ZLENGTH, ZERO, CYCLIC_Z_PD)

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         gradPg(I,J,K,:) = cPG - gradPg(I,J,K,:)

      ENDDO
      ENDDO
      ENDDO

      IF(DES_EXPLICITLY_COUPLED) THEN

! Calculate the gas phase forces acting on each particle.
         DO NP=1,MAX_PIP

            IF(NONEXISTENT==PARTICLE_STATE(NP) .or.                              &
               ENTERING_PARTICLE==PARTICLE_STATE(NP) .or. ENTERING_GHOST==PARTICLE_STATE(NP) .or.      &
               EXITING_PARTICLE==PARTICLE_STATE(NP)  .or. EXITING_GHOST==PARTICLE_STATE(NP)) CYCLE

            i = PIJK(NP,1)
            j = PIJK(NP,2)
            k = PIJK(NP,3)
            if (.NOT.fluid_at(i,j,k)) CYCLE

! Include gas pressure and gas-solids drag
            DRAG_FC(NP,:) = DRAG_FC(NP,:) + gradPg(i,j,k,:)*PVOL(NP)

         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CALC_PG_GRAD
