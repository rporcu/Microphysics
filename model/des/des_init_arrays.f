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

      USE param1, only: zero
      USE discretelement, only: des_rop_s, pinc, f_gds, vxf_gds, grav, drag_am, drag_bm, max_pip

      IMPLICIT NONE

!-----------------------------------------------

      PINC(:,:,:) = 0

      DES_ROP_S(:,:,:,:) = ZERO

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

         use discretelement, only: des_radius, ro_sol, pmass, omoi, des_pos_new, des_vel_new, omega_new, particle_state, pvol
         use discretelement, only: dg_pijk, dg_pijkprv, ighost_updated, neighbor_index, fc, tow, wall_collision_facet_id, pijk
         use discretelement, only: rot_acc_old, des_usr_var_size
         use discretelement, only: wall_collision_pft, iglobal_id, drag_fc, des_acc_old, nonexistent, do_old, des_usr_var
         use param1, only: zero

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
         particle_state(II) = nonexistent
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
      DRAG_FC(LB:UB,:) = ZERO

! Higher order time integration variables.
      IF (DO_OLD) THEN
         DES_ACC_OLD(LB:UB,:) = ZERO
         ROT_ACC_OLD(LB:UB,:) = ZERO
      ENDIF

      RETURN
      END SUBROUTINE DES_INIT_PARTICLE_ARRAYS
