MODULE MASS_INFLOW_DEM_MODULE

      use bc, only: bc_i_w, bc_j_s, bc_k_b
      use bc, only: bc_plane
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use constant, only: PI
      use des_allocate, only: particle_grow
      use des_bc, only: dem_bc_poly_layout, dem_bcmi_map, dem_bcmi_ijk, dem_bcmi_ijkstart, dem_mi_time, dem_bcmi_ijkend
      use des_bc, only: dem_mi, dem_bcmi, numfrac_limit, pi_count, pi_factor
      use desgrid, only: dg_funijk
      use desgrid, only: dg_iend, dg_jend, dg_kend
      use desgrid, only: dg_iend2, dg_jend2, dg_kend2
      use desgrid, only: dg_istart, dg_jstart, dg_kstart
      use desgrid, only: dg_istart2, dg_jstart2, dg_kstart2
      use desgrid, only: iofpos, jofpos, kofpos
      use discretelement, only: des_explicitly_coupled, xe, yn, zt
      use discretelement, only: max_pip, dtsolid, imax_global_id
      use discretelement, only: dg_pic, s_time, pip, normal_particle
      use discretelement, only: entering_particle, entering_ghost, normal_ghost, exiting_particle, exiting_ghost, nonexistent
      use param1, only: half, zero
      use constant, only: d_p0, ro_s0

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_MASS_INLET                                          !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entering the system.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MASS_INFLOW_DEM(pijk, particle_phase, dg_pijk, iglobal_id, particle_state, &
         des_radius, omoi, pmass, pvol, ro_sol, &
         des_vel_new, des_pos_new, ppos, omega_new, drag_fc)

      implicit none

      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: des_radius, omoi, pmass, pvol, ro_sol
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: des_vel_new, des_pos_new, ppos, omega_new, drag_fc
      INTEGER(KIND=1), DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijk, iglobal_id
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase

      INTEGER :: IP, LS, M, NP, IJK, LC
      INTEGER :: BCV, BCV_I
      LOGICAL :: CHECK_FOR_ERRORS, OWNS

! I/J/K index of fluid cell containing the new particle.
      INTEGER :: IJKP(3)

      DOUBLE PRECISION :: DIST, POS(3)
! Random numbers -shared by all ranks-
      DOUBLE PRECISION :: RAND(3)

      CHECK_FOR_ERRORS = .FALSE.

      DO BCV_I = 1, DEM_BCMI
         BCV = DEM_BCMI_MAP(BCV_I)

         DO LC=DEM_BCMI_IJKSTART(BCV_I), DEM_BCMI_IJKEND(BCV_I)
           IJK = DEM_BCMI_IJK(LC)

            DO LS= 1,DG_PIC(IJK)%ISIZE
               NP = DG_PIC(IJK)%P(LS)
               IF(EXITING_PARTICLE==PARTICLE_STATE(NP) .or. EXITING_GHOST==PARTICLE_STATE(NP)) CYCLE
               SELECT CASE (BC_PLANE(BCV))
               CASE('N'); DIST = DES_POS_NEW(NP,2) - YN(BC_J_s(BCV))
               CASE('S'); DIST = YN(BC_J_s(BCV)-1) - DES_POS_NEW(NP,2)
               CASE('E'); DIST = DES_POS_NEW(NP,1) - XE(BC_I_w(BCV))
               CASE('W'); DIST = XE(BC_I_w(BCV)-1) - DES_POS_NEW(NP,1)
               CASE('T'); DIST = DES_POS_NEW(NP,3) - ZT(BC_K_b(BCV))
               CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - DES_POS_NEW(NP,3)
               CASE DEFAULT
                  STOP __LINE__
               END SELECT
! The particle is still inside the domain
               IF(DIST > DES_RADIUS(NP)) THEN
                  IF(ENTERING_PARTICLE==PARTICLE_STATE(NP)) PARTICLE_STATE(NP) = NORMAL_PARTICLE
                  IF(ENTERING_GHOST==PARTICLE_STATE(NP)) PARTICLE_STATE(NP) = NORMAL_GHOST
               ENDIF
            ENDDO
         ENDDO

! All ranks generate random numbers but PE_IO BCASTS its values so that
! the calculated position (with random wiggles) is the same on all ranks.
         CALL RANDOM_NUMBER(RAND)
         ! !CALL BCAST(RAND)

! Check if any particles need seeded this time step.
         IF(DEM_MI_TIME(BCV_I) > S_TIME) CYCLE

         LS = 1

! Loop over the particles being injected
         PLoop: DO IP = 1, PI_COUNT(BCV_I)

! Increment the global max particle ID (all ranks).
            iMAX_GLOBAL_ID = iMAX_GLOBAL_ID + 1

! Determine the location and phase of the incoming particle.
            CALL SEED_NEW_PARTICLE(BCV, BCV_I, RAND, M, POS, IJKP, OWNS)

! Only the rank receiving the new particle needs to continue
            IF(.NOT.OWNS) CYCLE PLoop

! Increment the number of particle on the processor by one. If the max
! number of particles is exceeded, set the error flag and cycle.
            PIP = PIP + 1
            CALL PARTICLE_GROW(PIP)

! Find the first free space in the particle existance array.
            NP_LP: DO NP = LS, MAX_PIP
               IF(NONEXISTENT==PARTICLE_STATE(NP)) THEN
                  LS = NP
                  EXIT NP_LP
               ENDIF
            ENDDO NP_LP

! Set the particle's global ID.
            iGLOBAL_ID(LS) = iMAX_GLOBAL_ID

! Set the properties of the new particle.
            CALL SET_NEW_PARTICLE_PROPS(BCV, M, LS, POS, IJKP, pijk, particle_phase, dg_pijk, particle_state, &
               des_radius, omoi, pmass, pvol, ro_sol, &
               des_vel_new, des_pos_new, ppos, omega_new, drag_fc)

         ENDDO PLoop

! Update the time for seeding the next particle.
         DEM_MI_TIME(BCV_I) = S_TIME + PI_FACTOR(BCV_I)*DTSOLID
! Set the flag for error checking.
         CHECK_FOR_ERRORS = .TRUE.
      ENDDO

      IF(CHECK_FOR_ERRORS) THEN
      ENDIF

      RETURN
      END SUBROUTINE MASS_INFLOW_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SEED_NEW_PARTICLE                                       !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose:  This routine uses the classification information to place !
!  a new particle in the proper location.                              !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SEED_NEW_PARTICLE(lBCV, lBCV_I, lRAND, lM, lPOS,      &
         lIJKP, lOWNS)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! The associated bc no.
      INTEGER, INTENT(IN) :: lBCV, lBCV_I
! Random numbers
      DOUBLE PRECISION, INTENT(IN) :: lRAND(3)
! Phase of incoming particle.
      INTEGER, INTENT(OUT) :: lM
! Position of incoming particle.
      DOUBLE PRECISION, INTENT(OUT) :: lPOS(3)
! I/J/K index of fluid cell containing the new particle.
      INTEGER, INTENT(OUT) :: lIJKP(3)
! Logical indicating that the current rank is the owner
      LOGICAL, INTENT(OUT) :: lOWNS

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! the associated bc no.
!      INTEGER :: BCV
! Random position
      DOUBLE PRECISION RAND1, RAND2
! Array index of vacant position
      INTEGER :: VACANCY
! Number of array position.
      INTEGER :: OCCUPANTS

      INTEGER :: RAND_I
      INTEGER :: lI, lJ, lK

      DOUBLE PRECISION :: WINDOW

!      IF(PARTICLE_PLCMNT(lBCV_I) == 'ORDR')THEN
         VACANCY = DEM_MI(lBCV_I)%VACANCY
         OCCUPANTS = DEM_MI(lBCV_I)%OCCUPANTS
         DEM_MI(lBCV_I)%VACANCY = MOD(VACANCY,OCCUPANTS) + 1
!      ELSE
!         !call bcast(lpar_rad)
!         !call bcast(lpar_pos)
!         !call bcast(m)
!      ENDIF


! Assign the phase of the incoming particle.
      IF(DEM_MI(lBCV_I)%POLYDISPERSE) THEN
         RAND_I = ceiling(dble(NUMFRAC_LIMIT)*lRAND(1))
         lM = DEM_BC_POLY_LAYOUT(lBCV_I, RAND_I)
      ELSE
         lM = DEM_BC_POLY_LAYOUT(lBCV_I,1)
      ENDIF

      WINDOW = DEM_MI(lBCV_I)%WINDOW
      RAND1 = HALF*D_p0(lM) + (WINDOW - D_p0(lM))*lRAND(1)
      RAND2 = HALF*D_p0(lM) + (WINDOW - D_p0(lM))*lRAND(2)


! Set the physical location and I/J/K location of the particle.
      SELECT CASE(BC_PLANE(lBCV))
      CASE('N','S')

         lPOS(1) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(3) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(2) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(1) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(2) = DEM_MI(lBCV_I)%L

      CASE('E','W')

         lPOS(2) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(3) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(1) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(2) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(1) = DEM_MI(lBCV_I)%L

      CASE('T','B')

         lPOS(1) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(2) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(3) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(1) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(2) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%L

      END SELECT

! Easier and cleaner to clear out the third component at the end.

      lI = IofPOS(lPOS(1))
      lJ = JofPOS(lPOS(2))
      lK = KofPOS(lPOS(3))

      lOWNS = ((DG_ISTART <= lI) .AND. (lI <= DG_IEND) .AND.           &
         (DG_JSTART <= lJ) .AND. (lJ <= DG_JEND))

      lOWNS = lOWNS .AND. (DG_KSTART<=lK) .AND. (lK<=DG_KEND)

      RETURN
      END SUBROUTINE SEED_NEW_PARTICLE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_NEW_PARTICLE_PROPS                                  !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose:  Set the properties of the new particle.                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_NEW_PARTICLE_PROPS(lBCV, lM, lNP, lPOS, lIJKP, pijk, particle_phase, dg_pijk, particle_state, &
         des_radius, omoi, pmass, pvol, ro_sol, &
         des_vel_new, des_pos_new, ppos, omega_new, drag_fc)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! The associated bc no.
      INTEGER, INTENT(IN) :: lBCV
! Phase of incoming particle.
      INTEGER, INTENT(IN) :: lM
! Index of new particle
      INTEGER, INTENT(IN) :: lNP
! Position of incoming particle.
      DOUBLE PRECISION, INTENT(IN) :: lPOS(3)
! I/J/K index of fluid cell containing the new particle.
      INTEGER, INTENT(IN) :: lIJKP(3)

      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: des_radius, omoi, pmass, pvol, ro_sol
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijk
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: des_vel_new, des_pos_new, ppos, omega_new, drag_fc
      INTEGER(KIND=1), DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase

! I/J/K index of DES grid cell
      INTEGER :: lI, lJ, lK

! Global phase index
      INTEGER :: BC_M

      BC_M = lM

! The particle exists and is entering, not exiting nor a ghost particle
      IF (NORMAL_GHOST==PARTICLE_STATE(lNP)) THEN
         PARTICLE_STATE(lNP) = ENTERING_GHOST
      ELSE
         PARTICLE_STATE(lNP) = ENTERING_PARTICLE
      ENDIF

! Set the initial position values based on mass inlet class
      DES_POS_NEW(lNP,:) = lPOS(:)
      PPOS(lNP,:) = lPOS(:)
      DES_VEL_NEW(lNP,1) = BC_U_s(lBCV,BC_M)
      DES_VEL_NEW(lNP,2) = BC_V_s(lBCV,BC_M)
      DES_VEL_NEW(lNP,3) = BC_W_s(lBCV,BC_M)

! Set the initial angular velocity values
      OMEGA_NEW(lNP,:) = 0

! Set the particle radius value
      DES_RADIUS(lNP) = HALF * D_P0(lM)

! Set the particle density value
      RO_Sol(lNP) = RO_S0(lM)

! Store the I/J/K indices of the particle.
      PIJK(lNP,1:3) = lIJKP(:)

! Set the particle mass phase
      particle_phase(lNP) = lM

! Calculate the DES grid cell indices.
      lI = min(DG_IEND2,max(DG_ISTART2,IOFPOS(DES_POS_NEW(lNP,1))))
      lJ = min(DG_JEND2,max(DG_JSTART2,JOFPOS(DES_POS_NEW(lNP,2))))
      lK = min(DG_KEND2,max(DG_KSTART2,KOFPOS(DES_POS_NEW(lNP,3))))
! Store the triple
      DG_PIJK(lNP) = DG_FUNIJK(lI,lJ,lK)

! Calculate the new particle's Volume, Mass, OMOI
      PVOL(lNP) = (4.0d0/3.0d0) * PI * DES_RADIUS(lNP)**3
      PMASS(lNP) = PVOL(lNP) * RO_Sol(lNP)
      OMOI(lNP) = 5.d0 / (2.d0 * PMASS(lNP) * DES_RADIUS(lNP)**2)

! Clear the drag force
      IF(DES_EXPLICITLY_COUPLED) DRAG_FC(lNP,:) = ZERO


      RETURN
      END SUBROUTINE SET_NEW_PARTICLE_PROPS

END MODULE MASS_INFLOW_DEM_MODULE
