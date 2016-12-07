MODULE CFASSIGN_MODULE

   use compar, only: imap_c, jmap_c, kmap_c
   use compar, only: istart1,iend1,jstart1,jend1,kstart1,kend1
   use compar, only: istart2,iend2,jstart2,jend2,kstart2,kend2
   use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
   use des_allocate, only: des_vol_node
   use functions, only: fluid_at
   use geometry, only: dx, dy, dz
   use param1, only: zero

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! subroutine: compute_volume_of_nodes                                      C
! Author: Rahul Garg                                                       C
! Dare: Dec 20, 2011                                                       C
! Purpose:  Due to the staggered grid, the interpolaiton of mean fields    C
! is always first done at the nodes (of the scalar cell) for any quantity. C
! For a quantity like drag force or ep_s, one needs to assign a geometric  C
! volume to a node. In the past this was done on-the-fly in drag_fgs.f.    C
! VCELL was the variable that was used and it used to be incorrecty set to C
! the volume of the scalar cell that the particle belongs to.              C
! This will be ok for uniform grid and will not hold for non-uniform grids.C
! In addition, at the edges (in 2-D, the node volume is only half of       C
! the standard cell volume. And for corner (in 2-D), it is only one fourth.C
! This was also done on-the-fly earlier                                    C
! But again the volume of the cell was used, which was not correct         C
! also not extendable to cut-cell. So this routine computes the geoemetric C
! volume of the nodes.                                                     C
!                                                                          C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE compute_volume_of_nodes

      implicit none
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! ijk index of fluid grid and corresponding i, j, k indices
      integer :: iraw, jraw, kraw
! i, j, k and (i+1, j+1, k+1) indices corrected for any
! cyclic ghost boundaries on fluid grid
      INTEGER :: I, J, K, ip, jp, kp
! volume of indicated grid
      double precision :: vol_ijk, vol_ipjk, vol_ijpk, vol_ipjpk
      double precision :: vol_ijkp, vol_ipjkp, vol_ijpkp, vol_ipjpkp
      double precision :: vol_node_count, vol_node_actual_count
! weighting factor
      double precision :: avg_factor
! not used?
      double precision :: vol_node_uncorr
!-----------------------------------------------

      avg_factor = 0.125d0

! compute the volume at the grid nodes
! grid nodes start from istart2 to iend1
      vol_node_count = 8.

        DO Kraw = kstart3, kend3
        DO Jraw = jstart3, jend3
        DO Iraw = istart3, iend3

         des_vol_node(iraw,jraw,kraw) = zero

! start at 1 (ghost cell) and go to last fluid cell. why start at a
! ghost cell and not a fluid cell?
! See below

! Since we are interested in nodes making up the interpolation stencil,
! their numbering goes from 1 to iend1.
! Think of a case with IMAX = 3. Here the nodes where the interpolation will be
! done will run from 1 (=istart2) to 4 (=iend1).
         IF(iraw.LT.istart2 .OR. iraw.GT.iend1) CYCLE
         IF(jraw.LT.jstart2 .OR. jraw.GT.jend1) CYCLE
         IF(kraw.LT.kstart2 .OR. kraw.GT.kend1) CYCLE

! this will force indices of ghost cells on cyclic borders back to
! the corresponding fluid cells. since we are using i, j, k indices and
! not just a composite ijk index we need these to be shifted to account
! for periodicity
         I = imap_c(iraw)
         J = jmap_c(jraw)
         K = kmap_c(kraw)
         IP = imap_c(iraw+1)
         JP = jmap_c(jraw+1)

! the existing variable vol(ijk) is not used here for cut-cell reasons
! see r. garg for discussion
         vol_ijk   = dx*dy*dz
         vol_ipjk  = dx*dy*dz
         vol_ijpk  = dx*dy*dz
         vol_ipjpk = dx*dy*dz

         vol_node_uncorr = avg_factor*(vol_ijk + vol_ipjk + vol_ijpk + vol_ipjpk)
         vol_node_actual_count = vol_node_count

         IF(.NOT.fluid_at(iraw,jraw,kraw)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijk  = zero
         ENDIF

         IF(.NOT.fluid_at(ip,j,k)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjk  = zero
         ENDIF

         IF(.NOT.fluid_at(i,jp,k)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijpk  = zero
         ENDIF

         IF(.NOT.fluid_at(ip,jp,k)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjpk = zero
         ENDIF

! this will have non-zero values for non-fluid cells at the
! west/south/bottom borders but not for east/north/top borders?
         des_vol_node(iraw,jraw,kraw) = avg_factor*(vol_ijk + vol_ipjk + &
                                         vol_ijpk + vol_ipjpk)

         KP     = kmap_c(kraw+1)

         vol_ijkp   = dx*dy*dz
         vol_ipjkp  = dx*dy*dz
         vol_ijpkp  = dx*dy*dz
         vol_ipjpkp = dx*dy*dz

         vol_node_uncorr = avg_factor*(vol_node_uncorr + vol_ijkp + &
            vol_ipjkp + vol_ijpkp + vol_ipjpkp)

         IF(.NOT.fluid_at(i,j,kp)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijkp   = zero
         ENDIF

         IF(.NOT.fluid_at(ip,j,kp)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjkp  = zero
         ENDIF

         IF(.NOT.fluid_at(i,jp,kp)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijpkp  = zero
         ENDIF

         IF(.NOT.fluid_at(ip,jp,kp)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjpkp = zero
         ENDIF

         des_vol_node(iraw,jraw,kraw) = des_vol_node(iraw,jraw,kraw) + avg_factor*&
            (vol_ijkp + vol_ipjpkp + vol_ijpkp + vol_ipjkp)


      ENDDO
      ENDDO
      ENDDO

      RETURN
      RETURN
      END SUBROUTINE compute_volume_of_nodes

END MODULE CFASSIGN_MODULE
