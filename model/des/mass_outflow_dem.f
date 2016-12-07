MODULE MASS_OUTFLOW_DEM_MODULE

      USE bc, only: bc_i_w, bc_j_s, bc_k_b
      USE bc, only: bc_type, bc_plane
      USE bc, only: bc_u_s, bc_v_s, bc_w_s
      USE des_bc, only: dem_bcmo, dem_bcmo_map, dem_bcmo_ijk, dem_bcmo_ijkstart, dem_bcmo_ijkend
      USE discretelement, only: normal_ghost, normal_particle
      USE discretelement, only: particle_state, nonexistent, entering_ghost, exiting_ghost, exiting_particle, entering_particle
      USE discretelement, only: wall_collision_facet_id, dg_pic, pip
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
      SUBROUTINE MASS_OUTFLOW_DEM(FORCE_NSEARCH, pijk, iglobal_id, particle_state, &
         des_radius, omoi, pmass, pvol, ro_sol, &
         des_vel_new, des_pos_new, ppos, omega_new, fc, tow)

      implicit none

      LOGICAL, INTENT(INOUT) :: FORCE_NSEARCH

      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: des_radius, omoi, pmass, pvol, ro_sol
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: des_vel_new, des_pos_new, ppos, omega_new, fc, tow
      INTEGER(KIND=1), DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: iglobal_id
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

      INTEGER :: IJK
      INTEGER :: LC, LP, NP, M
      INTEGER :: BCV, BCV_I, IDX

      DOUBLE PRECISION :: SGN
      DOUBLE PRECISION :: DIST

      LOGICAL :: FREEZE_VEL
      DOUBLE PRECISION :: FREEZE(3)

      DO BCV_I = 1, DEM_BCMO

         BCV = DEM_BCMO_MAP(BCV_I)

         FREEZE_VEL = (BC_TYPE(BCV) /= 'MASS_OUTFLOW')

         SELECT CASE (BC_PLANE(BCV))
         CASE('N'); FREEZE = (/0.0d0, 1.0d0, 0.0d0/); IDX=2; SGN=-1.0d0
         CASE('S'); FREEZE = (/0.0d0, 1.0d0, 0.0d0/); IDX=2; SGN= 1.0d0
         CASE('E'); FREEZE = (/1.0d0, 0.0d0, 0.0d0/); IDX=1; SGN=-1.0d0
         CASE('W'); FREEZE = (/1.0d0, 0.0d0, 0.0d0/); IDX=1; SGN= 1.0d0
         CASE('T'); FREEZE = (/0.0d0, 0.0d0, 1.0d0/); IDX=3; SGN=-1.0d0
         CASE('B'); FREEZE = (/0.0d0, 0.0d0, 1.0d0/); IDX=3; SGN= 1.0d0
         END SELECT

         DO LC=DEM_BCMO_IJKSTART(BCV_I), DEM_BCMO_IJKEND(BCV_I)
            IJK = DEM_BCMO_IJK(LC)
            DO LP= 1,DG_PIC(IJK)%ISIZE

               NP = DG_PIC(IJK)%P(LP)

               IF(NONEXISTENT==PARTICLE_STATE(NP)) CYCLE
               if ((PARTICLE_STATE(NP)==NORMAL_GHOST) .OR. &
                  (PARTICLE_STATE(NP)==ENTERING_GHOST) .OR. &
                  (PARTICLE_STATE(NP)==EXITING_GHOST)) CYCLE

               IF(ENTERING_PARTICLE==PARTICLE_STATE(NP)) CYCLE

               SELECT CASE (BC_PLANE(BCV))
               CASE('S'); DIST = YN(BC_J_s(BCV)-1) - DES_POS_NEW(NP,2)
               CASE('N'); DIST = DES_POS_NEW(NP,2) - YN(BC_J_s(BCV))
               CASE('W'); DIST = XE(BC_I_w(BCV)-1) - DES_POS_NEW(NP,1)
               CASE('E'); DIST = DES_POS_NEW(NP,1) - XE(BC_I_w(BCV))
               CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - DES_POS_NEW(NP,3)
               CASE('T'); DIST = DES_POS_NEW(NP,3) - ZT(BC_K_b(BCV))
               END SELECT

! The particle is still inside the domain
               IF(DIST > DES_RADIUS(NP)) THEN
                  particle_state(NP) = normal_particle

! Check if the particle is crossing over the outlet plane.
               ELSEIF(DIST > ZERO) THEN

! The velocity is 'frozen' normal to the outflow plane. This approach
! is strict because complex BCs (via STLs) can let particles pop through
! the wall along the outlet.
                  IF(FREEZE_VEL) THEN
! Only 'freeze' a particle's velocy if it has it moving out of the
! domain. Otherwise, particles flagged as exiting but moving away from
! the BC appear to moon-walk through the domain until it crashes.
                     IF(DES_VEL_NEW(NP,IDX)*SGN > 0.0d0) THEN
                        DES_VEL_NEW(NP,:) = DES_VEL_NEW(NP,:)*FREEZE(:)
! Set the flags for an exiting particle.
                        IF (NORMAL_GHOST==PARTICLE_STATE(NP)) THEN
                           PARTICLE_STATE(NP) = EXITING_GHOST
                        ELSE
                           PARTICLE_STATE(NP) = EXITING_PARTICLE
                        ENDIF
                     ENDIF

! The user specified velocity is applied to the exiting particle. This
! only applies to mass outflows where the speed at which particles
! exit needs to be controled.
                  ELSE
                     M = PIJK(NP,5)
                     DES_VEL_NEW(NP,1) = BC_U_s(BCV,M)
                     DES_VEL_NEW(NP,2) = BC_V_s(BCV,M)
                     DES_VEL_NEW(NP,3) = BC_W_s(BCV,M)
! Set the flags for an exiting particle.
                     IF (NORMAL_GHOST==PARTICLE_STATE(NP)) THEN
                        particle_state(NP) = EXITING_GHOST
                     ELSE
                        particle_state(NP) = EXITING_PARTICLE
                     ENDIF
                  ENDIF

! Ladies and gentlemen, the particle has left the building.
               ELSE
                  CALL DELETE_PARTICLE(NP)
                  FORCE_NSEARCH = .TRUE.
               ENDIF

            ENDDO
         ENDDO
      ENDDO

! Sync the search flag across all processes.
!      CALL GLOBAL_ALL_OR(FORCE_NSEARCH)

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
