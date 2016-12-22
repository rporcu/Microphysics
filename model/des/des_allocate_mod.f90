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
      USE constant, only: mmax

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

      WRITE(ERR_MSG,1000) trim(iVal(MAX_PIP))
      CALL FLUSH_ERR_MSG(HEADER = .FALSE., FOOTER = .FALSE.)

 1000 FORMAT('Initial DES Particle array size: ',A)

! DES Allocatable arrays
!-----------------------------------------------

! allocate variables related to ghost particles
      allocate(ighost_updated(max_pip))

      Allocate(  wall_collision_facet_id (COLLISION_ARRAY_MAX, MAX_PIP) )
      wall_collision_facet_id(:,:) = -1

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

! Explicit drag force acting on a particle.
      Allocate(DRAG_FC (MAX_PIP,DIMN) )

! Volume of nodes
      ALLOCATE(DES_VOL_NODE(istart3:iend3, jstart3:jend3, kstart3:kend3))


      CALL FINL_ERR_MSG

      RETURN
   CONTAINS



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
         stop 17414
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
      ! SUBROUTINE NEIGHBOR_GROW(new_neigh_max)
      !   IMPLICIT NONE

      !   integer, intent(in) :: new_neigh_max

      !   INTEGER :: lSIZE1
      !   INTEGER, DIMENSION(:), ALLOCATABLE :: neigh_tmp
      !   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pf_tmp

      !   lSIZE1 = size(neighbors,1)

      !   allocate(neigh_tmp(new_neigh_max))
      !   neigh_tmp(1:lSIZE1) = neighbors(1:lSIZE1)
      !   neigh_tmp(lSIZE1+1:) = 0
      !   call move_alloc(neigh_tmp,neighbors)

      !   allocate(neigh_tmp(new_neigh_max))
      !   neigh_tmp(1:lSIZE1) = neighbors_old(1:lSIZE1)
      !   neigh_tmp(lSIZE1+1:) = 0
      !   call move_alloc(neigh_tmp,neighbors_old)

      !   allocate(pf_tmp(3,new_neigh_max))
      !   pf_tmp(:,1:lSIZE1) = pft_neighbor(:,1:lSIZE1)
      !   pf_tmp(:,lSIZE1+1:) = 0
      !   call move_alloc(pf_tmp,pft_neighbor)

      !   allocate(pf_tmp(3,new_neigh_max))
      !   pf_tmp(:,1:lSIZE1) = pft_neighbor_old(:,1:lSIZE1)
      !   pf_tmp(:,lSIZE1+1:) = 0
      !   call move_alloc(pf_tmp,pft_neighbor_old)

      ! END SUBROUTINE NEIGHBOR_GROW

!``````````````````````````````````````````````````````````````````````!
! Subroutine: PARTICLE_GROW                                            !
!                                                                      !
! Purpose: Grow particle arrays to new_max_pip. Note that pair         !
! max should be increased before calling this routine. Also, no        !
! assumption to the previous array size is made as needed for restarts.!
!``````````````````````````````````````````````````````````````````````!
      ! SUBROUTINE PARTICLE_GROW(new_max_pip)

      !   IMPLICIT NONE

      !   integer, intent(in) :: new_max_pip
      !   write(*,*) 'Death in particle grow'
      !   stop 887

      ! RETURN

      ! CONTAINS
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
      ! SUBROUTINE DES_INIT_PARTICLE_ARRAYS(LB,UB)

      !    use discretelement, only: ighost_updated, neighbor_index, wall_collision_facet_id
      !    use discretelement, only: drag_fc, nonexistent
      !    use param1, only: zero

      ! IMPLICIT NONE

      ! INTEGER, INTENT(IN) :: LB, UB

      ! NEIGHBOR_INDEX(:) = 0

      ! ! DES grid bin information
      ! IGHOST_UPDATED(LB:UB) = .false.

      ! ! Collision data
      ! WALL_COLLISION_FACET_ID(:,LB:UB) = -1

      ! ! Particle center drag coefficient and explicit drag force
      ! DRAG_FC(LB:UB,:) = ZERO

      ! RETURN
      ! END SUBROUTINE DES_INIT_PARTICLE_ARRAYS

      ! END SUBROUTINE PARTICLE_GROW

      SUBROUTINE BYTE_GROW(byte_array,new_size,ival)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: byte_array
        INTEGER, DIMENSION(:), ALLOCATABLE :: byte_tmp
        INTEGER :: ival ! initial value
        INTEGER lSIZE

        lSIZE = size(byte_array,1)
        allocate(byte_tmp(new_size))
        byte_tmp(1:lSIZE) = byte_array(1:lSIZE)
        call move_alloc(byte_tmp,byte_array)
        byte_array(new_size+1:) = ival

      END SUBROUTINE BYTE_GROW

      SUBROUTINE INTEGER_GROW(integer_array,new_size,ival)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:), ALLOCATABLE :: integer_tmp
        INTEGER :: ival
        INTEGER lSIZE

        lSIZE = size(integer_array,1)
        allocate(integer_tmp(new_size))
        integer_tmp(1:lSIZE) = integer_array(1:lSIZE)
        call move_alloc(integer_tmp,integer_array)
        integer_array(new_size+1:) = ival

      END SUBROUTINE INTEGER_GROW

      SUBROUTINE INTEGER_GROW2_reverse(integer_array,new_size,ival)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER :: ival
        INTEGER lSIZE, lSIZE2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        allocate(integer_tmp(new_size,lSIZE2))
        integer_tmp(1:lSIZE,:) = integer_array(1:lSIZE,:)
        call move_alloc(integer_tmp,integer_array)
        integer_array(new_size+1:,:) = ival

      END SUBROUTINE INTEGER_GROW2_reverse

      SUBROUTINE INTEGER_GROW2(integer_array,new_size,ival)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER :: ival
        INTEGER lSIZE, lSIZE2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        allocate(integer_tmp(lSIZE,new_size))
        integer_tmp(:,1:lSIZE2) = integer_array(:,1:lSIZE2)
        call move_alloc(integer_tmp,integer_array)
        integer_array(:,new_size+1:) = ival

      END SUBROUTINE INTEGER_GROW2

      SUBROUTINE LOGICAL_GROW(logical_array,new_size,ival)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:), ALLOCATABLE :: logical_tmp
        LOGICAL :: ival
        INTEGER lSIZE

        lSIZE = size(logical_array,1)
        allocate(logical_tmp(new_size))
        logical_tmp(1:lSIZE) = logical_array(1:lSIZE)
        call move_alloc(logical_tmp,logical_array)
        logical_array(new_size+1:) = ival

      END SUBROUTINE LOGICAL_GROW

      SUBROUTINE LOGICAL_GROW2(logical_array,new_size,ival)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: logical_tmp
        LOGICAL :: ival
        INTEGER lSIZE, lSIZE2

        lSIZE = size(logical_array,1)
        lSIZE2 = size(logical_array,2)
        allocate(logical_tmp(lSIZE,new_size))
        logical_tmp(:,1:lSIZE2) = logical_array(:,1:lSIZE2)
        call move_alloc(logical_tmp,logical_array)
        logical_array(:,new_size+1:) = ival

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
        real_array(new_size+1:) = ZERO

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
        real_array(:,new_size+1:) = ZERO

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
        real_array(new_size+1:,:) = ZERO

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
        real_array(:,:,new_size+1:) = ZERO

      END SUBROUTINE REAL_GROW3

      SUBROUTINE LOGICAL_GROW2_REVERSE(logical_array,new_size,ival)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: logical_tmp
        LOGICAL :: ival
        INTEGER lSIZE, lSIZE2

        lSIZE = size(logical_array,1)
        lSIZE2 = size(logical_array,2)
        allocate(logical_tmp(new_size,lSIZE2))
        logical_tmp(1:lSIZE,:) = logical_array(1:lSIZE,:)
        call move_alloc(logical_tmp,logical_array)
        logical_array(new_size+1:,:) = ival

      END SUBROUTINE LOGICAL_GROW2_REVERSE

    END MODULE DES_ALLOCATE
