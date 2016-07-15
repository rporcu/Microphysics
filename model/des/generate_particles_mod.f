      MODULE GENERATE_PARTICLES

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration based on maximum particle  C
!           radius and filling from top to bottom within specified     C
!           bounds                                                     C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 19-Mar-14  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GENERATE_PARTICLE_CONFIG

      use discretelement, only: PIP, PARTICLES
! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! Parameter for detecting unspecified values, zero, and one
      use param1, only: ONE
! Maximum number of initial conditions
      use param, only: DIMENSION_IC
! IC Region gas volume fraction.
      use ic, only: IC_EP_G

      use mpi_utility

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

      INTEGER :: ICV

! Initialize the error manager.
      CALL INIT_ERR_MSG("Generate_Particle_Config")

      DO ICV = 1, DIMENSION_IC

         IF(.NOT.IC_DEFINED(ICV)) CYCLE
         IF(IC_EP_G(ICV) == ONE) CYCLE

         CALL GENERATE_PARTICLE_CONFIG_DEM(ICV)

      ENDDO

      CALL GLOBAL_SUM(PIP,PARTICLES)

      WRITE(ERR_MSG, 1004) PARTICLES
 1004 FORMAT(/,'Total number of particles in the system: ',I15)

      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GENERATE_PARTICLE_CONFIG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
!  Authors: Rahul Garg                               Date: 21-Mar-2014 !
!                                                                      !
!  Purpose: Generate particle configuration for DEM solids for each IC !
!           region. Now using the particle linked lists for initial    !
!           build                                                      !
!           This routine will ultimately supersede the older rouine    !
!           that has not been deleted yet                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM(ICV)


! Global Variables:
! particle radius and density
      use discretelement, only: DES_RADIUS, RO_Sol
! particle position new and old
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! particle velocity new and old
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Simulation dimension (2D/3D)
      use discretelement, only: DIMN
! Number of particles in the system (current)
      use discretelement, only: PIP
! Number of DEM solids phases.
      use discretelement, only: DES_MMAX
! Flag to use _OLD variables
      use discretelement, only: DO_OLD
! Angular velocity
      use discretelement, only: OMEGA_OLD, OMEGA_NEW, PIJK
! DEM solid phase diameters and densities.
      use physprop, only: D_p0, RO_s0, SMAX
! IC Region solids volume fraction.
      use ic, only: IC_EP_S

! Constant: 3.14159...
      use constant, only: PI
! min and max physical co-ordinates of IC regions in each direction
      use ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
! initally specified velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s
! Flag to extend the lattice distribution in a given IC to available area
      use ic, only: IC_DES_FIT_TO_REGION
! Parameter for detecting unspecified values, zero, and one
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, Half
! Parameter for small and large numbers
      use param1, only: SMALL_NUMBER, LARGE_NUMBER

! to access random number generator subroutines
      use randomno
      use mpi_utility
      use functions, only: SET_NORMAL, FUNIJK, FLUID_AT

      use error_manager

! direction wise spans of the domain and grid spacing in each direction
      use geometry, only: xlength, ylength, zlength

      use cutcell, only : CARTESIAN_GRID, CUT_CELL_AT
      use STL_PREPROC_DES, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      use run, only: solids_model
      use des_allocate, only: PARTICLE_GROW

      use discretelement, only: MAX_RADIUS

      use discretelement, only: XE, YN, ZT

      use param, only: DIM_M, DIMENSION_I, DIMENSION_J, DIMENSION_K
      use functions, only: IS_ON_MYPE_WOBND

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICV

! Local variables
!---------------------------------------------------------------------//
! Starting positions in the axial directions
      DOUBLE PRECISION :: xINIT, yINIT, zINIT
! Fractor used to scale particle diameter
      DOUBLE PRECISION :: lFAC
! Particle position and velocity
      DOUBLE PRECISION :: POS(3), VEL(3)
! Number of particles in the lattice
      INTEGER :: SEED_X, SEED_Y, SEED_Z
! Loop indices phase/fluid cell
      INTEGER :: M, MM, I, J, K, IJK, LB, UB
! Loop indicies for seeding
      INTEGER :: II, JJ, KK
! Start and end bound for IC region.
      DOUBLE PRECISION :: IC_START(3), IC_END(3)
! Volume and lengths of the IC Region
      DOUBLE PRECISION :: DOM_VOL, DOML(3)
! Flag to skip the particle
      LOGICAL :: SKIP
! Diameter adjusted for space padding
      DOUBLE PRECISION :: ADJ_DIA
! Number of particles calculated from volume fracton
      INTEGER :: rPARTS(DIM_M), tPARTS
! Spacing between particles.
      DOUBLE PRECISION :: lDEL, lDX, lDY, lDZ
! Flag that the setup failed to fit the particles to the IC region
      LOGICAL :: FIT_FAILED
! Number of seeded particles
      INTEGER :: pCOUNT(DIM_M), tCOUNT

      DOUBLE PRECISION :: SOLIDS_DATA(0:DIM_M)

      LOGICAL :: VEL_FLUCT
      DOUBLE PRECISION :: VEL_SIG
      DOUBLE PRECISION, ALLOCATABLE :: randVEL(:,:)

!......................................................................!

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_DEM")

      WRITE(ERR_MSG,"(2/,'Generating initial particle configuration:')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      SOLIDS_DATA = ZERO
      CALL GET_IC_VOLUME(ICV, SOLIDS_DATA(0))

! setting particle seed spacing grid to be slightly greater than
! the maximum particle diameter. seed at ~particle radii
      lFAC = 1.05D0

! Setup local arrays with IC region bounds.
      IC_START(1)=IC_X_W(ICV);   IC_END(1)=IC_X_E(ICV)
      IC_START(2)=IC_Y_S(ICV);   IC_END(2)=IC_Y_N(ICV)
      IC_START(3)=IC_Z_B(ICV);   IC_END(3)=IC_Z_T(ICV)

      DOML = IC_END-IC_START
      IF(NO_K) DOML(3)=DZ(1)

! Volume of the IC region
      DOM_VOL = DOML(1)*DOML(2)*DOML(3)

      rPARTS=0
      VEL_FLUCT = .FALSE.
      DO M=1,SMAX+DES_MMAX
         IF(SOLIDS_MODEL(M) == 'DEM') THEN
! Number of particles for phase M
            rPARTS(M) = &
               floor((6.0d0*IC_EP_S(ICV,M)*DOM_VOL)/(PI*(D_P0(M)**3)))
         ENDIF
      ENDDO

! Total number of particles in this IC region.
      tPARTS = sum(rPARTS)
      IF(tPARTS == 0) RETURN

      ADJ_DIA = 2.0d0*MAX_RADIUS*lFAC


! Attempt to seed particle throughout the IC region
      FIT_FAILED=.FALSE.
      IF(IC_DES_FIT_TO_REGION(ICV)) THEN
         IF(NO_K) THEN
            lDEL = (DOML(1)-ADJ_DIA)*(DOML(2)-ADJ_DIA)
            lDEL = (lDEL/dble(tPARTS))**(1.0/2.0)
            SEED_X = max(1,ceiling((DOML(1)-ADJ_DIA)/lDEL))
            SEED_Y = max(1,ceiling((DOML(2)-ADJ_DIA)/lDEL))
            SEED_Z = 1
         ELSE
            lDEL = (DOML(1)-ADJ_DIA)*(DOML(2)-ADJ_DIA)*(DOML(3)-ADJ_DIA)
            lDEL = (lDEL/dble(tPARTS))**(1.0/3.0)
            SEED_X = max(1,ceiling((DOML(1)-ADJ_DIA)/lDEL))
            SEED_Y = max(1,ceiling((DOML(2)-ADJ_DIA)/lDEL))
            SEED_Z = max(1,ceiling((DOML(3)-ADJ_DIA)/lDEL))
         ENDIF
         FIT_FAILED=(dble(SEED_X*SEED_Y*SEED_Z) < tPARTS)
      ENDIF

! Generic filling
      IF(.NOT.IC_DES_FIT_TO_REGION(ICV) .OR. FIT_FAILED) THEN
         SEED_X = max(1,floor((IC_END(1)-IC_START(1)-ADJ_DIA)/ADJ_DIA))
         SEED_Y = max(1,floor((IC_END(2)-IC_START(2)-ADJ_DIA)/ADJ_DIA))
         SEED_Z = max(1,floor((IC_END(3)-IC_START(3)-ADJ_DIA)/ADJ_DIA))
      ENDIF

      lDX = DOML(1)/dble(SEED_X)
      lDY = DOML(2)/dble(SEED_Y)
      IF(DO_K) THEN
         lDZ = DOML(3)/dble(SEED_Z)
      ELSE
         lDZ = 0.0d0
      ENDIF

      xINIT = IC_START(1)+HALF*lDX
      yINIT = IC_START(2)+HALF*lDY
      zINIT = IC_START(3)+HALF*lDZ


      M=1
      pCOUNT = 0
      tCOUNT = 0
      JJ_LP: DO  JJ=1, SEED_Y
      KK_LP: DO  KK=1, SEED_Z
      II_LP: DO  II=1, SEED_X

! Exit if all particles were seeded.
         IF(tCOUNT > int(tPARTS)) THEN
            EXIT JJ_LP
! Find the next phase that needs to be seeded
         ELSEIF(pCOUNT(M) > int(rPARTS(M))) THEN
            MM_LP: DO MM=M+1,SMAX+DES_MMAX
               IF(rPARTS(MM) > 0.0) THEN
                  M=MM
                  EXIT MM_LP
               ENDIF
            ENDDO MM_LP
            IF(MM > SMAX+DES_MMAX) THEN
               EXIT JJ_LP
            ENDIF
         ENDIF

         pCOUNT(M) = pCOUNT(M) + 1
         tCOUNT = tCOUNT + 1

! Position the particle
         POS(1) = xINIT + (II-1)*lDX
         POS(2) = YINIT + (JJ-1)*lDY
         POS(3) = ZINIT + (KK-1)*lDZ

! Bin the parcel to the fuild grid.
         K=1
         IF(DO_K) CALL PIC_SEARCH(K, POS(3), ZT, DIMENSION_K, KMIN2, KMAX2)
         CALL PIC_SEARCH(J, POS(2), YN, DIMENSION_J, JMIN2, JMAX2)
         CALL PIC_SEARCH(I, POS(1), XE, DIMENSION_I, IMIN2, IMAX2)

! Skip cells that are not part of the local fuild domain.
         IF(.NOT.IS_ON_MYPE_WOBND(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

         IJK = FUNIJK(I,J,K)
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IF(CARTESIAN_GRID) THEN
            CALL CHECK_IF_PARTICLE_OVERLAPS_STL(POS, I, J, K, SKIP)
            IF(SKIP) CYCLE
         ENDIF

         PIP = PIP + 1
         CALL PARTICLE_GROW(PIP)

         CALL SET_NORMAL(PIP)

         IF(VEL_FLUCT) THEN
            VEL(1) = randVEL(tCOUNT,1)
            VEL(2) = randVEL(tCOUNT,2)
            VEL(3) = randVEL(tCOUNT,3)
         ELSE
            VEL(1) = IC_U_s(ICV,M)
            VEL(2) = IC_V_s(ICV,M)
            VEL(3) = IC_W_s(ICV,M)
         ENDIF
         IF(NO_K) VEL(3) = 0.0d0


         DES_POS_NEW(:,PIP) = POS(:)
         DES_VEL_NEW(:,PIP) = VEL(:)
         OMEGA_NEW(:,PIP) = 0.0d0

         DES_RADIUS(PIP) = D_P0(M)*HALF
         RO_SOL(PIP) =  RO_S0(M)

         PIJK(PIP,1) = I
         PIJK(PIP,2) = J
         PIJK(PIP,3) = K
         PIJK(PIP,4) = IJK
         PIJK(PIP,5) = M

         IF(DO_OLD) THEN
            DES_VEL_OLD(:,PIP) = DES_VEL_NEW(:,PIP)
            DES_POS_OLD(:,PIP) = DES_POS_NEW(:,PIP)
            OMEGA_OLD(:,PIP) = ZERO
         ENDIF

         SOLIDS_DATA(M) = SOLIDS_DATA(M) + 1.0

      ENDDO II_LP
      ENDDO KK_LP
      ENDDO JJ_LP

! Collect the data
      CALL GLOBAL_ALL_SUM(SOLIDS_DATA)

! Verify that the IC region volume is not zero.
      IF(SOLIDS_DATA(0) <= 0.0d0) THEN
         WRITE(ERR_MSG,1000) ICV, SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 
1000 FORMAT('Error 1000: Invalid IC region volume: IC=',I3,' VOL=',&
         ES15.4,/'Please correct the mfix.dat file.')

      WRITE(ERR_MSG,2000) ICV
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      DO M=1, SMAX+DES_MMAX
         IF(SOLIDS_DATA(M) < SMALL_NUMBER) CYCLE
         WRITE(ERR_MSG,2010) M, int(SOLIDS_DATA(M)), IC_EP_S(ICV,M),   &
            (dble(SOLIDS_DATA(M))*(Pi/6.0d0)*D_P0(M)**3)/SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO

      IF(allocated(randVEL)) deallocate(randVEL)

      CALL FINL_ERR_MSG

      RETURN

 2000 FORMAT(/2x,'|',43('-'),'|',/2x,'| IC Region: ',I3,28x,'|',/2x,   &
         '|',43('-'),'|',/2x,'| Phase | Number of |    EPs    |    EP',&
         's    |',/2x,'|   ID  | Particles | Specified |   Actual  |', &
         /2x,'|-------|',3(11('-'),'|'))

 2010 FORMAT(2x,'|  ',I3,'  |',1x,I9,1x,'|',2(1x,ES9.2,1x,'|'),/2x,    &
         '|-------|',3(11('-'),'|'))

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GET_IC_VOLUME                                           !
!  Author: J.Musser                                 Date: 26-Aug-2015  !
!                                                                      !
!  Purpose: Calculate the actual volume of the IC region.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GET_IC_VOLUME(ICV, IC_VOL)

! IC region index bounds
      use ic, only: IC_I_w, IC_I_e
      use ic, only: IC_J_s, IC_J_n
      use ic, only: IC_K_b, IC_K_t

! Volume of computational cells.
      use geometry, only: VOL

      use functions

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Index of IC region
      INTEGER, INTENT(IN) :: ICV
! Total calculated volume of IC region
      DOUBLE PRECISION, INTENT(OUT) :: IC_VOL

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I, J, K, IJK
!......................................................................!


      IC_VOL = 0.0d0
      DO K = IC_K_B(ICV), IC_K_T(ICV)
      DO J = IC_J_S(ICV), IC_J_N(ICV)
      DO I = IC_I_W(ICV), IC_I_E(ICV)

         IF(.NOT.IS_ON_MYPE_WOBND(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

         IJK = FUNIJK(I,J,K)
         IF(FLUID_AT(IJK)) IC_VOL = IC_VOL + VOL(IJK)

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GET_IC_VOLUME

      END MODULE GENERATE_PARTICLES
