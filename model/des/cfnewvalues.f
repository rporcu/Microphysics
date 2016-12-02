!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,
!           position, angular velocity etc
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES

      USE discretelement, only: des_acc_old, rot_acc_old, fc, tow, des_vel_new, des_pos_new, omega_new, ppos, do_nsearch, dtsolid
      USE discretelement, only: max_pip, intg_euler, omoi, intg_adams_bashforth, grav, des_radius, pmass, neighbor_search_rad_ratio
      USE discretelement, only: entering_particle, entering_ghost, nonexistent, exiting_particle, exiting_ghost, particle_state
      USE discretelement, only: normal_ghost
      USE param1, only: zero

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L
      DOUBLE PRECISION :: DD(3), NEIGHBOR_SEARCH_DIST
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

      DOUBLE PRECISION :: lVELo(3), lPOSo(3)

!-----------------------------------------------

! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE                       ! Only real particles
            IF(ENTERING_PARTICLE==PARTICLE_STATE(L).or.ENTERING_GHOST==PARTICLE_STATE(L)) CYCLE  ! Only non-entering
            IF(NORMAL_GHOST==PARTICLE_STATE(L)) CYCLE                             ! Skip ghost particles
            DES_ACC_OLD(L,:) = FC(L,:)/PMASS(L) + GRAV(:)
            ROT_ACC_OLD(L,:) = TOW(L,:)
         ENDDO
      ENDIF
      DO L = 1, MAX_PIP
! only process particles that exist
         IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE
! skip ghost particles
         IF(NORMAL_GHOST==PARTICLE_STATE(L).or.ENTERING_GHOST==PARTICLE_STATE(L).or.EXITING_GHOST==PARTICLE_STATE(L)) CYCLE

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.ENTERING_PARTICLE==PARTICLE_STATE(L) .AND. .NOT.ENTERING_GHOST==PARTICLE_STATE(L))THEN
            FC(L,:) = FC(L,:)/PMASS(L) + GRAV(:)
         ELSE
            FC(L,:) = ZERO
            TOW(L,:) = ZERO
         ENDIF

! Advance particle position, velocity
        IF (INTG_EULER) THEN
! first-order method
           DES_VEL_NEW(L,:) = DES_VEL_NEW(L,:) + FC(L,:)*DTSOLID
           DD(:) = DES_VEL_NEW(L,:)*DTSOLID
           DES_POS_NEW(L,:) = DES_POS_NEW(L,:) + DD(:)
            OMEGA_NEW(L,:)   = OMEGA_NEW(L,:) + TOW(L,:)*OMOI(L)*DTSOLID
         ELSEIF (INTG_ADAMS_BASHFORTH) THEN

            lVELo = DES_VEL_NEW(L,:)
            lPOSo = DES_POS_NEW(L,:)

! Second-order Adams-Bashforth/Trapezoidal scheme
            DES_VEL_NEW(L,:) = lVELo(:) + 0.5d0*&
               ( 3.d0*FC(L,:)-DES_ACC_OLD(L,:) )*DTSOLID

            OMEGA_NEW(L,:)   =  OMEGA_NEW(L,:) + 0.5d0*&
               ( 3.d0*TOW(L,:)*OMOI(L)-ROT_ACC_OLD(L,:) )*DTSOLID

            DD(:) = 0.5d0*( lVELo(:)+DES_VEL_NEW(L,:) )*DTSOLID

            DES_POS_NEW(L,:) = lPOSo(:) + DD(:)
            DES_ACC_OLD(L,:) = FC(L,:)
            ROT_ACC_OLD(L,:) = TOW(L,:)*OMOI(L)
         ENDIF


! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
         IF(.NOT.DO_NSEARCH) THEN
            DD(:) = DES_POS_NEW(L,:) - PPOS(L,:)
            NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*&
               DES_RADIUS(L)
            IF(dot_product(DD,DD).GE.NEIGHBOR_SEARCH_DIST**2) DO_NSEARCH = .TRUE.
         ENDIF

! Reset total contact force and torque
         FC(L,:) = ZERO
         TOW(L,:) = ZERO

      ENDDO

      FIRST_PASS = .FALSE.

      RETURN

      contains

        include 'functions.inc'

      END SUBROUTINE CFNEWVALUES
