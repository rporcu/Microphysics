MODULE CALC_PG_GRAD_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

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
      SUBROUTINE CALC_PG_GRAD(slo, shi, lo, hi, &
                              p_g, gradPg,  particle_state, des_pos_new,&
                              pvol, drag_fc, flag, dx, dy, dz)

      use calc_grad_des_module, only: calc_grad_des
      use discretelement, only: entering_particle, entering_ghost
      use discretelement, only: exiting_particle, exiting_ghost
      use discretelement, only: nonexistent

      ! Flags for cyclic BC with pressure drop
      use geometry, only: CYCLIC_X_PD, CYCLIC_Y_PD, CYCLIC_Z_PD

      ! Specified pressure drop
      use bc, only: DELP_X, DELP_Y, DELP_Z

      ! Domain length
      use geometry, only: XLENGTH, YLENGTH, ZLENGTH

      use discretelement, only: MAX_PIP, DES_EXPLICITLY_COUPLED

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO

      implicit none

      integer, intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(out  ) :: gradpg&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      integer         , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(c_real), intent(in   ) :: pvol(:)
      real(c_real), intent(in   ) :: des_pos_new(:,:)
      real(c_real), intent(inout) :: drag_fc(:,:)
      real(c_real), intent(in   ) :: dx, dy, dz

      integer         , intent(in   ) :: particle_state(:)

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, I, J, K
! mean pressure gradient for the case of periodic boundaries
      real(c_real) :: cPG(3)
! One over cell volume
      real(c_real) :: odx, ody, odz
!......................................................................!

! Calculate the gas phase pressure gradient. (dP/dx)
      CALL CALC_GRAD_DES(slo, shi, lo, hi, P_G, gradPg, flag, dx, dy, dz)

! Add in cyclic BC pressure drop.
      cPG(1) = merge(DELP_X/XLENGTH, ZERO, CYCLIC_X_PD)
      cPG(2) = merge(DELP_Y/YLENGTH, ZERO, CYCLIC_Y_PD)
      cPG(3) = merge(DELP_Z/ZLENGTH, ZERO, CYCLIC_Z_PD)

      DO K = slo(3),shi(3)
      DO J = slo(2),shi(2)
      DO I = slo(1),shi(1)

         gradPg(I,J,K,:) = cPG - gradPg(I,J,K,:)

      ENDDO
      ENDDO
      ENDDO

      IF(DES_EXPLICITLY_COUPLED) THEN

         odx = 1.0d0/dx
         ody = 1.0d0/dy
         odz = 1.0d0/dz

! Calculate the gas phase forces acting on each particle.
         DO NP=1,MAX_PIP

            if(nonexistent==particle_state(np) .or.        &
               entering_particle==particle_state(np) .or.  &
               entering_ghost==particle_state(np) .or.     &
               exiting_particle==particle_state(np)  .or.  &
               exiting_ghost==particle_state(np)) cycle

! Fluid cell containing the particle
            i = floor(des_pos_new(np,1)*odx) + 1
            j = floor(des_pos_new(np,2)*ody) + 1
            k = floor(des_pos_new(np,3)*odz) + 1

            if (.NOT.1.eq.flag(i,j,k,1)) CYCLE

! Include gas pressure and gas-solids drag
            DRAG_FC(NP,:) = DRAG_FC(NP,:) + gradPg(i,j,k,:)*PVOL(NP)

         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CALC_PG_GRAD

END MODULE CALC_PG_GRAD_MODULE
