MODULE CFNEWVALUES_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,
!           position, angular velocity etc
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES(particle_state, des_radius, pmass, omoi, &
         des_pos_new, des_vel_new, omega_new, fc, tow, &
         des_acc_old, rot_acc_old)

      USE discretelement, only: dtsolid
      USE discretelement, only: max_pip, intg_euler, intg_adams_bashforth
      USE discretelement, only: entering_particle, entering_ghost
      USE discretelement, only: nonexistent, exiting_ghost
      USE discretelement, only: normal_ghost
      USE param1, only: zero
      use constant, only: gravity

      IMPLICIT NONE

      INTEGER         , INTENT(IN   ) :: particle_state(:)
      real(c_real), INTENT(IN   ) :: des_radius(:)
      real(c_real), INTENT(IN   ) :: omoi(:)
      real(c_real), INTENT(IN   ) :: pmass(:)
      real(c_real), INTENT(INOUT) :: des_pos_new(:,:)
      real(c_real), INTENT(INOUT) :: des_vel_new(:,:)
      real(c_real), INTENT(INOUT) :: omega_new(:,:)
      real(c_real), INTENT(INOUT) :: fc(:,:)
      real(c_real), INTENT(INOUT) :: tow(:,:)
      real(c_real), INTENT(INOUT) :: des_acc_old(:,:)
      real(c_real), INTENT(INOUT) :: rot_acc_old(:,:)

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L
      real(c_real) :: DD(3)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

      real(c_real) :: lVELo(3), lPOSo(3)

!-----------------------------------------------

! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE
            IF(ENTERING_PARTICLE==PARTICLE_STATE(L).or.&
               ENTERING_GHOST==PARTICLE_STATE(L)) CYCLE  ! Only non-entering
            IF(NORMAL_GHOST==PARTICLE_STATE(L)) CYCLE
            DES_ACC_OLD(L,:) = FC(L,:)/PMASS(L) + GRAVITY(:)
            ROT_ACC_OLD(L,:) = TOW(L,:)
         ENDDO
      ENDIF
      DO L = 1, MAX_PIP
! only process particles that exist
         IF(NONEXISTENT==PARTICLE_STATE(L) .or.      &
            NORMAL_GHOST==PARTICLE_STATE(L).or.      &
            ENTERING_GHOST==PARTICLE_STATE(L).or.    &
            EXITING_GHOST==PARTICLE_STATE(L)) CYCLE

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.ENTERING_PARTICLE==PARTICLE_STATE(L) .AND. &
            .NOT.ENTERING_GHOST==PARTICLE_STATE(L))THEN
            FC(L,:) = FC(L,:)/PMASS(L) + GRAVITY(:)
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

! Reset total contact force and torque
         FC(L,:) = ZERO
         TOW(L,:) = ZERO

      ENDDO

      FIRST_PASS = .FALSE.

      RETURN

      contains

        include 'functions.inc'

      END SUBROUTINE CFNEWVALUES

END MODULE CFNEWVALUES_MODULE
