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

      SUBROUTINE NEIGHBOUR(pijk, pinc, particle_state, des_radius, &
         des_pos_new, ppos, neighbor_index, neighbor_index_old)

      USE discretelement, only: do_nsearch, max_pip, neighbors_old
      USE discretelement, only: neighbors
      USE desgrid, only: desgrid_neigh_build

      IMPLICIT NONE

      integer, intent(in) :: pijk(:,:)
      integer, intent(in) :: pinc(:,:,:)
      integer, intent(inout) :: particle_state(:)
      integer, intent(inout) :: neighbor_index(:)
      integer, intent(inout) :: neighbor_index_old(:)

      double precision, intent(inout) :: des_radius(:)
      double precision, intent(in   ) :: des_pos_new(:,:)
      double precision, intent(inout) :: ppos(:,:)

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      INTEGER :: cc,ll,cc_start,cc_end,cc_start_old,cc_end_old,cc_old
      LOGICAL :: found

!-----------------------------------------------
! Reset PPOS and NEIGHBOURS back to initialized values
      PPOS(:,:) = DES_POS_NEW(:,:)
      neighbor_index_old(:) = neighbor_index(:)

      return

      do cc=1, size(neighbors)
         neighbors_old(cc) = neighbors(cc)
      enddo

      NEIGHBOR_INDEX(:) = 0

      call desgrid_neigh_build(des_pos_new, pijk, pinc, particle_state, &
         des_radius, neighbor_index)

      do ll = 1, max_pip

         CC_START = 1
         IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
         CC_END   = NEIGHBOR_INDEX(LL)

         CC_START_OLD = 1
         IF (LL.gt.1) CC_START_OLD = NEIGHBOR_INDEX_OLD(LL-1)
         CC_END_OLD   = NEIGHBOR_INDEX_OLD(LL)

      enddo

! resetting do_nsearch to false here since neighbor search will have
! just been invoked
      DO_NSEARCH = .FALSE.

      RETURN
      END SUBROUTINE NEIGHBOUR
END MODULE NEIGHBOUR_MODULE
