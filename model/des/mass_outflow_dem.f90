MODULE MASS_OUTFLOW_DEM_MODULE

      USE bc, only: bc_i_w, bc_j_s, bc_k_b
      USE bc, only: bc_type, bc_plane
      USE bc, only: bc_u_s, bc_v_s, bc_w_s
      USE des_bc, only: dem_bcmo, dem_bcmo_map, dem_bcmo_ijk, dem_bcmo_ijkstart, dem_bcmo_ijkend
      USE discretelement, only: normal_ghost, normal_particle
      USE discretelement, only: nonexistent, entering_ghost, exiting_ghost, exiting_particle, entering_particle
      USE discretelement, only: wall_collision_facet_id, pip
      USE discretelement, only: xe, yn, zt
      USE param1, only: zero

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MASS_OUTFLOW_DEM                                        !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entering the system.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MASS_OUTFLOW_DEM(FORCE_NSEARCH, particle_phase, iglobal_id, particle_state, &
         des_radius, omoi, pmass, pvol, ro_sol, &
         des_vel_new, des_pos_new, ppos, omega_new, fc, tow)

      implicit none

      LOGICAL, INTENT(INOUT) :: FORCE_NSEARCH

      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: des_radius, omoi, pmass, pvol, ro_sol
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: des_vel_new, des_pos_new, ppos, omega_new, fc, tow
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: iglobal_id
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase

      INTEGER :: IJK
      INTEGER :: LC, LP, NP, M
      INTEGER :: BCV, BCV_I, IDX

      DOUBLE PRECISION :: SGN
      DOUBLE PRECISION :: DIST

      LOGICAL :: FREEZE_VEL
      DOUBLE PRECISION :: FREEZE(3)

! HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
      return
! HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK


      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DELETE_PARTICLE                                         !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine is used to check if a new particle has fully !
!  entered the domain.  If so, the flag classifying the particle as new!
!  is removed, allowing the particle to respond to contact forces from !
!  walls and other particles.                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DELETE_PARTICLE(NP)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      iGLOBAL_ID(NP) = -1
      PARTICLE_STATE(NP) = NONEXISTENT

      DES_POS_NEW(NP,:) = ZERO
      DES_VEL_NEW(NP,:) = ZERO
      OMEGA_NEW(NP,:) = ZERO

      DES_RADIUS(NP) = ZERO
      PMASS(NP) = ZERO
      PVOL(NP) = ZERO
      RO_Sol(NP) = ZERO
      OMOI(NP) = ZERO

      FC(NP,:) = ZERO
      TOW(NP,:) = ZERO

      PPOS(NP,:) = ZERO

      WALL_COLLISION_FACET_ID(:,NP) = -1

      PIP = PIP - 1

      RETURN
      END SUBROUTINE DELETE_PARTICLE
   END SUBROUTINE MASS_OUTFLOW_DEM

END MODULE MASS_OUTFLOW_DEM_MODULE
