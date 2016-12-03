!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_ALLOCATE                                           C
!                                                                      C
!  Purpose: subroutines to allocate all DEM arrays                     C
!                                                                      C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE DES_ALLOCATE

      USE compar, only: iend3, jend3, kend3
      USE compar, only: istart3, jstart3, kstart3
      USE compar, only: numpes
      USE des_bc, only: numfrac_limit
      USE des_bc, only: pi_factor, pi_count, dem_mi_time, dem_mi, dem_bc_poly_layout, dem_bcmi_ijkstart, dem_bcmi_ijkend, dem_bcmi
      USE discretelement
      USE error_manager, only: err_msg, ival, flush_err_msg, finl_err_msg, init_err_msg
      USE param1, only: undefined_i
      USE particle_filter, only: DES_INTERP_DPVM
      USE particle_filter, only: DES_INTERP_GARG
      USE particle_filter, only: DES_INTERP_GAUSS
      USE particle_filter, only: DES_INTERP_SCHEME_ENUM
      USE physprop, only: mmax

  PUBLIC:: DES_ALLOCATE_ARRAYS, ADD_PAIR, PARTICLE_GROW, ALLOCATE_DEM_MI

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_ALLOCATE_ARRAYS                                     C
!  Purpose: Original allocate arrays subroutines for DES               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_ALLOCATE_ARRAYS

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K
!-----------------------------------------------

      CALL INIT_ERR_MSG("DES_ALLOCATE_ARRAYS")

! For parallel processing the array size required should be either
! specified by the USEr or could be determined from total particles
! with some factor.
      MAX_PIP = merge(0, PARTICLES/numPEs, PARTICLES==UNDEFINED_I)
      MAX_PIP = MAX(MAX_PIP,4)

      WRITE(ERR_MSG,1000) trim(iVal(MAX_PIP))
      CALL FLUSH_ERR_MSG(HEADER = .FALSE., FOOTER = .FALSE.)

 1000 FORMAT('Initial DES Particle array size: ',A)

! DES Allocatable arrays
!-----------------------------------------------
! Dynamic particle info including another index for parallel
! processing for ghost
      ALLOCATE( PARTICLE_STATE (MAX_PIP) )
      ALLOCATE (iglobal_id(max_pip))

! Particle attributes
! Radius, density, mass, moment of inertia
      Allocate(  DES_RADIUS (MAX_PIP) )
      Allocate(  RO_Sol (MAX_PIP) )
      Allocate(  PVOL (MAX_PIP) )
      Allocate(  PMASS (MAX_PIP) )
      Allocate(  OMOI (MAX_PIP) )

! Old and new particle positions, velocities (translational and
! rotational)
      Allocate(  DES_POS_NEW (MAX_PIP,DIMN) )
      Allocate(  DES_VEL_NEW (MAX_PIP,DIMN) )
      Allocate(  OMEGA_NEW (MAX_PIP,DIMN) )

      IF (DO_OLD) THEN
         Allocate(  DES_ACC_OLD (MAX_PIP,DIMN) )
         Allocate(  ROT_ACC_OLD (MAX_PIP,DIMN))
      ENDIF

! Allocating user defined array
      IF(DES_USR_VAR_SIZE > 0) &
         Allocate( DES_USR_VAR(MAX_PIP,DES_USR_VAR_SIZE) )

! Particle positions at the last call neighbor search algorithm call
      Allocate(  PPOS (MAX_PIP,DIMN) )

! Total, normal and tangetial forces
      Allocate(  FC (MAX_PIP,DIMN) )

! Torque
      Allocate(  TOW (MAX_PIP,DIMN) )


! allocate variable for des grid binning
      allocate(dg_pijk(max_pip)); dg_pijk=0
      allocate(dg_pijkprv(max_pip)); dg_pijkprv=0

! allocate variables related to ghost particles
      allocate(ighost_updated(max_pip))



      Allocate(  wall_collision_facet_id (COLLISION_ARRAY_MAX, MAX_PIP) )
      wall_collision_facet_id(:,:) = -1
      Allocate(  wall_collision_PFT (DIMN, COLLISION_ARRAY_MAX, MAX_PIP) )

      NEIGH_MAX = MAX_PIP

      Allocate(  NEIGHBOR_INDEX (MAX_PIP) )
      Allocate(  NEIGHBOR_INDEX_OLD (MAX_PIP) )
      Allocate(  NEIGHBORS (NEIGH_MAX) )
      NEIGHBORS(:) = 0
      Allocate(  NEIGHBORS_OLD (NEIGH_MAX) )
      Allocate(  PFT_NEIGHBOR (3,NEIGH_MAX) )
      Allocate(  PFT_NEIGHBOR_OLD (3,NEIGH_MAX) )

! Variable that stores the particle in cell information (ID) on the
! computational fluid grid defined by imax, jmax and kmax in mfix.dat
      ALLOCATE(PIC(istart3:iend3, jstart3:jend3, kstart3:kend3))
      DO K = kstart3, kend3
         DO J = jstart3, jend3
            DO I = istart3, iend3
               NULLIFY(pic(i,j,k)%p)
            ENDDO
         ENDDO
      ENDDO

! Particles in a computational fluid cell (for volume fraction)
      Allocate(  PINC (istart3:iend3, jstart3:jend3, kstart3:kend3))

! For each particle track its i,j,k location on computational fluid grid
! defined by imax, jmax and kmax in mfix.dat and phase no.
      Allocate(  PIJK (MAX_PIP,5) )

      ALLOCATE(DRAG_AM(istart3:iend3, jstart3:jend3, kstart3:kend3))
      ALLOCATE(DRAG_BM(istart3:iend3, jstart3:jend3, kstart3:kend3, DIMN))

! Explicit drag force acting on a particle.
      Allocate(DRAG_FC (MAX_PIP,DIMN) )

! Volume of nodes
      ALLOCATE(DES_VOL_NODE(istart3:iend3, jstart3:jend3, kstart3:kend3))

      ALLOCATE(F_GDS(istart3:iend3, jstart3:jend3, kstart3:kend3))

      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_GARG)
         ALLOCATE(DES_ROPS_NODE(istart3:iend3, jstart3:jend3, kstart3:kend3, MMAX))
         ALLOCATE(DES_VEL_NODE(istart3:iend3, jstart3:jend3, kstart3:kend3, DIMN, MMAX))
      END SELECT

! Bulk density in a computational fluid cell / for communication with
! MFIX continuum
      ALLOCATE( DES_ROP_S(istart3:iend3, jstart3:jend3, kstart3:kend3, MMAX))

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE DES_ALLOCATE_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_DEM_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_DEM_MI

      IMPLICIT NONE
!-----------------------------------------------

! Particle injection factor
      Allocate( PI_FACTOR (DEM_BCMI) )
! Particle injection count (injection number)
      Allocate( PI_COUNT (DEM_BCMI) )
! Particle injection time scale
      Allocate( DEM_MI_TIME (DEM_BCMI) )
! Array used for polydisperse inlets: stores the particle number
! distribution of an inlet scaled with numfrac_limit
      Allocate( DEM_BC_POLY_LAYOUT( DEM_BCMI, NUMFRAC_LIMIT ) )
! Data structure for storing BC data.
      Allocate( DEM_MI(DEM_BCMI) )

! Initializiation
! Integer arrays
      PI_FACTOR(:) = -1
      PI_COUNT(:) = -1
      DEM_BC_POLY_LAYOUT(:,:) = -1
! Double precision arrays
      DEM_MI_TIME(:) = UNDEFINED

      allocate( DEM_BCMI_IJKSTART(DEM_BCMI) )
      allocate( DEM_BCMI_IJKEND(DEM_BCMI) )

      DEM_BCMI_IJKSTART = -1
      DEM_BCMI_IJKEND   = -1

! Boundary classification
!         Allocate( PARTICLE_PLCMNT (DES_BCMI) )
! Character precision arrays
!         PARTICLE_PLCMNT(:) = UNDEFINED_C

      RETURN
      END SUBROUTINE ALLOCATE_DEM_MI


!``````````````````````````````````````````````````````````````````````!
! Subroutine: ADD_PAIR                                                 !
!                                                                      !
! Purpose: Adds a neighbor pair to the pairs array.                    !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      INTEGER FUNCTION add_pair(ii,jj)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ii,jj

      if (NEIGHBOR_INDEX(ii) > NEIGH_MAX) then
         NEIGH_MAX = 2*NEIGH_MAX
         CALL NEIGHBOR_GROW(NEIGH_MAX)
      endif

      NEIGHBORS(NEIGHBOR_INDEX(ii)) = jj
      NEIGHBOR_INDEX(ii) = NEIGHBOR_INDEX(ii) + 1
      add_pair = NEIGHBOR_INDEX(ii)

      RETURN
      END FUNCTION add_pair

!``````````````````````````````````````````````````````````````````````!
! Subroutine: NEIGHBOR_GROW                                            !
!                                                                      !
! Purpose: Grow neighbors arrays to new_neigh_max. Note that neighbor      !
! max should be increased before calling this routine. Also, no        !
! assumption to the previous array size is made as needed for restarts.!
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE NEIGHBOR_GROW(new_neigh_max)
        IMPLICIT NONE

        integer, intent(in) :: new_neigh_max

        INTEGER :: lSIZE1
        INTEGER, DIMENSION(:), ALLOCATABLE :: neigh_tmp
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pf_tmp

        lSIZE1 = size(neighbors,1)

        allocate(neigh_tmp(new_neigh_max))
        neigh_tmp(1:lSIZE1) = neighbors(1:lSIZE1)
        neigh_tmp(lSIZE1+1:) = 0
        call move_alloc(neigh_tmp,neighbors)

        allocate(neigh_tmp(new_neigh_max))
        neigh_tmp(1:lSIZE1) = neighbors_old(1:lSIZE1)
        neigh_tmp(lSIZE1+1:) = 0
        call move_alloc(neigh_tmp,neighbors_old)

        allocate(pf_tmp(3,new_neigh_max))
        pf_tmp(:,1:lSIZE1) = pft_neighbor(:,1:lSIZE1)
        pf_tmp(:,lSIZE1+1:) = 0
        call move_alloc(pf_tmp,pft_neighbor)

        allocate(pf_tmp(3,new_neigh_max))
        pf_tmp(:,1:lSIZE1) = pft_neighbor_old(:,1:lSIZE1)
        pf_tmp(:,lSIZE1+1:) = 0
        call move_alloc(pf_tmp,pft_neighbor_old)


      END SUBROUTINE NEIGHBOR_GROW

!``````````````````````````````````````````````````````````````````````!
! Subroutine: PARTICLE_GROW                                            !
!                                                                      !
! Purpose: Grow particle arrays to new_max_pip. Note that pair         !
! max should be increased before calling this routine. Also, no        !
! assumption to the previous array size is made as needed for restarts.!
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE PARTICLE_GROW(new_max_pip)

        IMPLICIT NONE

        integer, intent(in) :: new_max_pip

        DO WHILE (MAX_PIP < new_max_pip)
           MAX_PIP = MAX_PIP*2

           call real_grow(des_radius,MAX_PIP)
           call real_grow(RO_Sol,MAX_PIP)
           call real_grow(PVOL,MAX_PIP)
           call real_grow(PMASS,MAX_PIP)
           call real_grow(OMOI,MAX_PIP)
           call real_grow2_reverse(DES_POS_NEW,MAX_PIP)
           call real_grow2_reverse(DES_VEL_NEW,MAX_PIP)
           call real_grow2_reverse(OMEGA_NEW,MAX_PIP)
           call real_grow2_reverse(PPOS,MAX_PIP)
           call byte_grow(PARTICLE_STATE,MAX_PIP)
           call integer_grow(iglobal_id,MAX_PIP)
           call integer_grow2_reverse(pijk,MAX_PIP)
           call integer_grow(dg_pijk,MAX_PIP)
           call integer_grow(dg_pijkprv,MAX_PIP)
           call logical_grow(ighost_updated,MAX_PIP)
           call real_grow2_reverse(FC,MAX_PIP)
           call real_grow2_reverse(TOW,MAX_PIP)
           call integer_grow2(WALL_COLLISION_FACET_ID,MAX_PIP)
           call real_grow3(WALL_COLLISION_PFT,MAX_PIP)
           call real_grow2_reverse(DRAG_FC,MAX_PIP)

           call integer_grow(NEIGHBOR_INDEX,MAX_PIP)
           call integer_grow(NEIGHBOR_INDEX_OLD,MAX_PIP)


           IF (DO_OLD) THEN
              call real_grow2_reverse(DES_ACC_OLD,MAX_PIP)
              call real_grow2_reverse(ROT_ACC_OLD,MAX_PIP)
           ENDIF

           IF(DES_USR_VAR_SIZE > 0) &
              call real_grow2_reverse(DES_USR_VAR,MAX_PIP)

           CALL DES_INIT_PARTICLE_ARRAYS(MAX_PIP/2+1,MAX_PIP)

        ENDDO

      RETURN

      END SUBROUTINE PARTICLE_GROW

      SUBROUTINE BYTE_GROW(byte_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: byte_array
        INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE :: byte_tmp
        INTEGER lSIZE

        lSIZE = size(byte_array,1)
        allocate(byte_tmp(new_size))
        byte_tmp(1:lSIZE) = byte_array(1:lSIZE)
        call move_alloc(byte_tmp,byte_array)

      END SUBROUTINE BYTE_GROW

      SUBROUTINE INTEGER_GROW(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE

        lSIZE = size(integer_array,1)
        allocate(integer_tmp(new_size))
        integer_tmp(1:lSIZE) = integer_array(1:lSIZE)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW

      SUBROUTINE INTEGER_GROW2_reverse(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        allocate(integer_tmp(new_size,lSIZE2))
        integer_tmp(1:lSIZE,:) = integer_array(1:lSIZE,:)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW2_reverse

      SUBROUTINE INTEGER_GROW2(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        allocate(integer_tmp(lSIZE,new_size))
        integer_tmp(:,1:lSIZE2) = integer_array(:,1:lSIZE2)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW2

      SUBROUTINE LOGICAL_GROW(logical_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:), ALLOCATABLE :: logical_tmp
        INTEGER lSIZE

        lSIZE = size(logical_array,1)
        allocate(logical_tmp(new_size))
        logical_tmp(1:lSIZE) = logical_array(1:lSIZE)
        call move_alloc(logical_tmp,logical_array)

      END SUBROUTINE LOGICAL_GROW

      SUBROUTINE LOGICAL_GROW2(logical_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: logical_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(logical_array,1)
        lSIZE2 = size(logical_array,2)
        allocate(logical_tmp(lSIZE,new_size))
        logical_tmp(:,1:lSIZE2) = logical_array(:,1:lSIZE2)
        call move_alloc(logical_tmp,logical_array)

      END SUBROUTINE LOGICAL_GROW2

      SUBROUTINE REAL_GROW(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE

        lSIZE = size(real_array,1)
        allocate(real_tmp(new_size))
        real_tmp(1:lSIZE) = real_array(1:lSIZE)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW

      SUBROUTINE REAL_GROW2(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(lSIZE,new_size))
        real_tmp(:,1:lSIZE2) = real_array(:,1:lSIZE2)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW2

      SUBROUTINE REAL_GROW2_reverse(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(new_size,lSIZE2))
        real_tmp(1:lSIZE,:) = real_array(1:lSIZE,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW2_REVERSE

      SUBROUTINE REAL_GROW3(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2, lSIZE3

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lSIZE3 = size(real_array,3)
        allocate(real_tmp(lSIZE,lSIZE2,new_size))
        real_tmp(:,:,1:lSIZE3) = real_array(:,:,1:lSIZE3)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW3

      SUBROUTINE LOGICAL_GROW2_REVERSE(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(new_size,lSIZE2))
        real_tmp(1:lSIZE,:) = real_array(1:lSIZE,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE LOGICAL_GROW2_REVERSE

    END MODULE DES_ALLOCATE
