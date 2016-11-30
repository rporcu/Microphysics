!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                         !
!  Subrourtine: DES_INIT_ARRAYS                                           !
!  Author: Jay Boyalakuntla                              Date: 12-Jun-04  !
!                                                                         !
!  Purpose: Initialize arrays at the start of the simulation. Note that   !
!  arrays based on the number of particles (MAX_PIP) should be added to   !
!  the DES_INIT_PARTICLE_ARRAYS as they need to be reinitialized after    !
!  particle arrays are grown (see PARTICLE_GROW).                         !
!                                                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_INIT_ARRAYS

      USE param
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE physprop
      USE des_bc
      USE run
      use desgrid
      use desmpi

      IMPLICIT NONE

!-----------------------------------------------

      PINC(:) = 0

      DES_ROP_S(:,:) = ZERO

      P_FORCE(:,:) = ZERO

      IF(allocated(DRAG_AM)) DRAG_AM = ZERO
      IF(allocated(DRAG_BM)) DRAG_BM = ZERO

      F_GDS = ZERO
      VXF_GDS = ZERO

      GRAV(:) = ZERO

      CALL DES_INIT_PARTICLE_ARRAYS(1,MAX_PIP)

      RETURN
      END SUBROUTINE DES_INIT_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                         !
!  Subrourtine: DES_INIT_PARTICLE_ARRAYS                                  !
!  Author: Jay Boyalakuntla                              Date: 12-Jun-04  !
!                                                                         !
!  Purpose: Initialize particle arrays. The upper and lower bounds are    !
!  passed so that after resizing particle arrays (see GROW_PARTICLE) the  !
!  new portions of the arrays can be initialized.                         !
!                                                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_INIT_PARTICLE_ARRAYS(LB,UB)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use desgrid
      use desmpi
      use discretelement
      use functions

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LB, UB
      INTEGER :: II

      IGLOBAL_ID(LB:UB) = 0
      PARTICLE_STATE(LB:UB) = NONEXISTENT

! Physical properties:
      DES_RADIUS(LB:UB) = ZERO
      RO_Sol(LB:UB) = ZERO
      PVOL(LB:UB) = ZERO
      PMASS(LB:UB) = ZERO
      OMOI(LB:UB) = ZERO

! Particle position, velocity, etc
      DES_POS_NEW(LB:UB,:) = ZERO
      DES_VEL_NEW(LB:UB,:) = ZERO
      OMEGA_NEW(LB:UB,:) = ZERO

! Particle state flag
      DO II = LB, UB
         call set_nonexistent(II)
      ENDDO
      NEIGHBOR_INDEX(:) = 0

! DES grid bin information
      DG_PIJK(LB:UB) = -1
      DG_PIJKPRV(LB:UB) = -1
      IGHOST_UPDATED(LB:UB) = .false.

! Fluid cell bin information
      PIJK(LB:UB,:) = 0

! Translation and rotational forces
      FC(LB:UB,:) = ZERO
      TOW(LB:UB,:) = ZERO

! Collision data
      WALL_COLLISION_FACET_ID(:,LB:UB) = -1
      WALL_COLLISION_PFT(:,:,LB:UB) = ZERO

! Initializing user defined array
      IF(DES_USR_VAR_SIZE > 0) &
         DES_USR_VAR(LB:UB,:) = ZERO

! Particle center drag coefficient and explicit drag force
      F_GP(LB:UB) = ZERO
      DRAG_FC(LB:UB,:) = ZERO

! Higher order time integration variables.
      IF (DO_OLD) THEN
         DES_ACC_OLD(LB:UB,:) = ZERO
         ROT_ACC_OLD(LB:UB,:) = ZERO
      ENDIF

      RETURN
      END SUBROUTINE DES_INIT_PARTICLE_ARRAYS
