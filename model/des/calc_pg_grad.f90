MODULE CALC_PG_GRAD_MODULE
   CONTAINS

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
      SUBROUTINE CALC_PG_GRAD(p_g, gradPg,  particle_state, des_pos_new,&
         pvol, drag_fc, flag)

      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      use calc_grad_des_module, only: calc_grad_des
      use discretelement, only: entering_particle, entering_ghost
      use discretelement, only: exiting_particle, exiting_ghost
      use discretelement, only: nonexistent

! Loop bounds for fluid grid
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flags for cyclic BC with pressure drop
      use geometry, only: CYCLIC_X_PD, CYCLIC_Y_PD, CYCLIC_Z_PD
      use geometry, only: dx, dy, dz
! Specified pressure drop
      use bc, only: DELP_X, DELP_Y, DELP_Z
! Domain length
      use geometry, only: XLENGTH, YLENGTH, ZLENGTH
! Gas phase pressure

      use discretelement, only: MAX_PIP, DES_EXPLICITLY_COUPLED

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO

      implicit none

      double precision, intent(in   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(out  ) :: gradpg&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      integer         , intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

      double precision, intent(in   ) :: pvol(:)
      double precision, intent(in   ) :: des_pos_new(:,:)
      double precision, intent(inout) :: drag_fc(:,:)
      integer         , intent(in   ) :: particle_state(:)

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, I, J, K
! mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION :: cPG(3)
! One over cell volume
      double precision :: Oodx, Oody, Oodz
!......................................................................!

! Calculate the gas phase pressure gradient. (dP/dx)
      CALL CALC_GRAD_DES(P_G, gradPg, flag)

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

         Oodx = 1.0d0/dx
         Oody = 1.0d0/dy
         Oodz = 1.0d0/dz

! Calculate the gas phase forces acting on each particle.
         DO NP=1,MAX_PIP

            if(nonexistent==particle_state(np) .or.        &
               entering_particle==particle_state(np) .or.  &
               entering_ghost==particle_state(np) .or.     &
               exiting_particle==particle_state(np)  .or.  &
               exiting_ghost==particle_state(np)) cycle

! Fluid cell containing the particle
            i = floor(des_pos_new(np,1)*Oodx) + 1
            j = floor(des_pos_new(np,2)*Oody) + 1
            k = floor(des_pos_new(np,3)*Oodz) + 1

            if (.NOT.1.eq.flag(i,j,k,1)) CYCLE

! Include gas pressure and gas-solids drag
            DRAG_FC(NP,:) = DRAG_FC(NP,:) + gradPg(i,j,k,:)*PVOL(NP)

         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CALC_PG_GRAD

END MODULE CALC_PG_GRAD_MODULE
