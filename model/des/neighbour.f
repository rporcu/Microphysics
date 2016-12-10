MODULE NEIGHBOUR_MODULE
   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NEIGHBOUR                                              C
!  Purpose: DES - Neighbors search;
!           N-Square,
!           Quadtree(2D)/Octree(3D)  (use at own risk)
!           Cell linked
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NEIGHBOUR

      USE discretelement, only: ppos, neighbor_index_old, pft_neighbor_old, neighbor_index, do_nsearch, max_pip, neighbors_old
      USE discretelement, only: des_pos_new, pft_neighbor, neighbors, des_radius, dg_pijk, particle_state
      USE desgrid, only: desgrid_neigh_build

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

INTEGER :: cc,ll,cc_start,cc_end,cc_start_old,cc_end_old,cc_old
LOGICAL :: found

!-----------------------------------------------
! Reset PPOS and NEIGHBOURS back to initialized values
      PPOS(:,:) = DES_POS_NEW(:,:)
      neighbor_index_old(:) = neighbor_index(:)

      do cc=1, size(neighbors)
         neighbors_old(cc) = neighbors(cc)
         pft_neighbor_old(:,cc) = pft_neighbor(:,cc)
      enddo

      NEIGHBOR_INDEX(:) = 0

      CALL DESGRID_NEIGH_BUILD(des_pos_new, dg_pijk, particle_state, des_radius)

      do ll = 1, max_pip

         CC_START = 1
         IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
         CC_END   = NEIGHBOR_INDEX(LL)

         CC_START_OLD = 1
         IF (LL.gt.1) CC_START_OLD = NEIGHBOR_INDEX_OLD(LL-1)
         CC_END_OLD   = NEIGHBOR_INDEX_OLD(LL)

         DO CC = CC_START, CC_END-1
            found = .false.
            DO CC_OLD = CC_START_OLD, CC_END_OLD-1
               if (neighbors(cc) .eq. neighbors_old(cc_old)) then
                  pft_neighbor(:,cc) = pft_neighbor_old(:,cc_old)
                  found = .true.
                  exit
               endif
            enddo

            if (.not.found) pft_neighbor(:,cc) = 0.0
         enddo
      enddo

! resetting do_nsearch to false here since neighbor search will have
! just been invoked
      DO_NSEARCH = .FALSE.

      RETURN
      END SUBROUTINE NEIGHBOUR
END MODULE NEIGHBOUR_MODULE
