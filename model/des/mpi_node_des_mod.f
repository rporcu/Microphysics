!------------------------------------------------------------------------
! Module           : desmpi
! Purpose          : Contains wrapper class for mpi communications- send,recv
!
! Author           : Pradeep.G
!
! Purpose          : Module contains subroutines and variables related to
!                    des mpi communication.
!
! Comments         : do_nsearch flag should be set to true before calling
!                    des_par_exchange; when do_nsearch is true ghost particles of the
!                    system will be updated, which will be later used to generate
!                    neighbour list.
!
! Contains following subroutines:
!    des_addnodevalues, des_addnodevalues2, des_addnodevalues_mean_fields
!------------------------------------------------------------------------
      module mpi_node_des

      use compar, only: iend2, jend2, kend2
      use compar, only: iend3, jend3, kend3
      use compar, only: istart2, jstart2, kstart2
      use compar, only: istart3, jstart3, kstart3
      use compar, only: nodesi, nodesj, nodesk
      use desgrid, only: des_periodic_walls_x, des_periodic_walls_y, des_periodic_walls_z
      use discretelement, only: des_rops_node, des_vel_node
      use functions, only: imax1, jmax1, kmax1
      use param, only: DIMENSION_N_s

      contains

!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues_mean_fields
! Purpose          : This routine is specially used for computing mean
!                    fields by backward interpolation.
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues_mean_fields()

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: li, lj, lk
!-----------------------------------------------

! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            des_rops_node(1,lj,lk,:)  = des_rops_node(1,lj,lk,:)+des_rops_node(imax1,lj,lk,:)
            des_vel_node(1,lj,lk,:,:) = des_vel_node(1,lj,lk,:,:)+des_vel_node(imax1,lj,lk,:,:)
            des_rops_node(imax1,lj,lk,:)  = des_rops_node(1,lj,lk,:)
            des_vel_node(imax1,lj,lk,:,:) = des_vel_node(1,lj,lk,:,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            des_rops_node(li,1,lk,:)  = des_rops_node(li,1,lk,:)+des_rops_node(li,jmax1,lk,:)
            des_vel_node(li,1,lk,:,:) = des_vel_node(li,1,lk,:,:)+des_vel_node(li,jmax1,lk,:,:)
            des_rops_node(li,jmax1,lk,:)  = des_rops_node(li,1,lk,:)
            des_vel_node(li,jmax1,lk,:,:) = des_vel_node(li,1,lk,:,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            des_rops_node(li,lj,1,:)  = des_rops_node(li,lj,1,:)+des_rops_node(li,lj,kmax1,:)
            des_vel_node(li,lj,1,:,:) = des_vel_node(li,lj,1,:,:)+des_vel_node(li,lj,kmax1,:,:)
            des_rops_node(li,lj,kmax1,:)  = des_rops_node(li,lj,1,:)
            des_vel_node(li,lj,kmax1,:,:) = des_vel_node(li,lj,1,:,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues_mean_fields



!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues
! Purpose          : This routine is specially used for des_drag_gs
!                    The backward interpolation in des_drag_gs computes
!                    the grid node values of drag_am and drag_bm
!                    node values are from istart2 to iend1;
!                    hence a separate module is created to exchange
!                    node values
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues(drag_am, drag_bm)

!-----------------------------------------------
      implicit none
!-----------------------------------------------
      DOUBLE PRECISION, INTENT(INOUT) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

! local variables
!-----------------------------------------------
      integer :: li, lj, lk
!-----------------------------------------------


! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            drag_am(1,lj,lk) = drag_am(1,lj,lk)+drag_am(imax1,lj,lk)
            drag_bm(1,lj,lk,:) = drag_bm(1,lj,lk,:)+drag_bm(imax1,lj,lk,:)
            drag_am(imax1,lj,lk) = drag_am(1,lj,lk)
            drag_bm(imax1,lj,lk,:) = drag_bm(1,lj,lk,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            drag_am(li,1,lk) = drag_am(li,1,lk)+drag_am(li,jmax1,lk)
            drag_bm(li,1,lk,:) = drag_bm(li,1,lk,:)+drag_bm(li,jmax1,lk,:)
            drag_am(li,jmax1,lk) = drag_am(li,1,lk)
            drag_bm(li,jmax1,lk,:) = drag_bm(li,1,lk,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            drag_am(li,lj,1) = drag_am(li,lj,1)+drag_am(li,lj,kmax1)
            drag_bm(li,lj,1,:) = drag_bm(li,lj,1,:)+drag_bm(li,lj,kmax1,:)
            drag_am(li,lj,kmax1) = drag_am(li,lj,1)
            drag_bm(li,lj,kmax1,:) = drag_bm(li,lj,1,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues


!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues2
! Purpose          : This routine is specially used for calc_des_rop_s
!                    The backward interpolation in calc_des_rop_s computes
!                    the grid node values of des_rops_node
!                    node values are from istart2 to iend1;
!                    hence a separate module is created to exchange
!                    node values
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues2()

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: li, lj, lk
!-----------------------------------------------

! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            des_rops_node(1,lj,lk,:) = des_rops_node(1,lj,lk,:)+des_rops_node(imax1,lj,lk,:)
            des_rops_node(imax1,lj,lk,:) = des_rops_node(1,lj,lk,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            des_rops_node(li,1,lk,:) = des_rops_node(li,1,lk,:)+des_rops_node(li,jmax1,lk,:)
            des_rops_node(li,jmax1,lk,:) = des_rops_node(li,1,lk,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            des_rops_node(li,lj,1,:) = des_rops_node(li,lj,1,:)+des_rops_node(li,lj,kmax1,:)
            des_rops_node(li,lj,kmax1,:) = des_rops_node(li,lj,1,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues2



      end module mpi_node_des
