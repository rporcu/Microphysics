!-------------------------------------------------------------------!
!                                                                   !
! Purpose:                                                          !
!    Variables to be declared for parallel information. Removed     !
!    from geometry_mod and put in COMPAR module with some           !
!    additional variables used by AEOLUS                            !
!                                                                   !
! Added by Ed and Sreekanth on 06/22/99.                            !
!-------------------------------------------------------------------!

      MODULE compar

!-----------------------------------------------
! Modules
!-----------------------------------------------

! myPE - my processor id (it varies from 0 to nproc-1)
! numPEs - total number of nodes
      integer, parameter :: myPE = 0
      integer, parameter :: numPEs = 1

! specify the rank of the PE to be used for I/O
      INTEGER :: PE_IO = 0

! nodesi, nodesj and nodesk represent the number of nodes
! in i, j, k directions respectively.
! nodesj = 1 (No decomposition along j-direction)
! For 1-D decomposition, nodesk = nproc for a 3d problem and
! nodesi = nproc for a 2D problem.
      integer :: nodesi, nodesj, nodesk

! root represents the 'root' processor. For now it is defaulted to
! zero
      integer :: root=0

! nlayers_bicgs - Number of layers for send_recv in bicgs
      integer :: nlayers_bicgs = 1

! Variables used for mapping i, j, k to ii, jj, kk to take care of
! of cyclic conditions...
      integer, allocatable,dimension(:) :: imap, jmap, kmap

      integer :: &
                istart3, iend3, jstart3, jend3, &
                kstart3, kend3, istart2, iend2, jstart2, jend2, &
                kstart2, kend2, istart1, iend1, jstart1, jend1, &
                kstart1, kend1

      integer :: istart, iend, jstart, jend, kstart, kend

      integer :: c0, c1, c2

      END MODULE compar
