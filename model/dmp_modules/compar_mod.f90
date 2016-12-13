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

! mpierr - used by AEOLUS for error checking
      INTEGER :: mpierr

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


! -istart1_all contains the starting i value for all the processors
!  excluding the ghost regions. istart2_all is for one extra ghost
!  layer and istart3_all is for two ghost layers. Similarly
!  iend1_all, iend2_all and iend3_all contain the ending values.
!  Similarly for j and k, jstart..., kstart.... are  prescribed.

! -All the variables without the '_all' represent that processor
!  values. So ijkstart3 denotes the starting value of ijk, which
!  belongs to the processor = funijk(istart3_all(myid),
!  jstart3_all(myid), kstart3_all(myid) for a 1-d decompostion of a
!  3D problem. For more details see gridmap_mod.f90. Similarly the
!  end values are denoted by ijkend3_all

! -displs has the necessary shift information to position the buffer
!  in the scatterv and gatherv routines.

! -ijksize3 is the size of the element owned by each processor plus
!  the ghost regions.
! -'_all' has information about all the processor mapping based on
!  above convention
      integer, allocatable,dimension(:) ::  &
                ijkstart3_all,ijkend3_all,    &
                istart_all,istart1_all,istart2_all,istart3_all, &
                jstart_all,jstart1_all,jstart2_all,jstart3_all, &
                kstart_all,kstart1_all,kstart2_all,kstart3_all, &
                iend_all,iend1_all,iend2_all,iend3_all, &
                jend_all,jend1_all,jend2_all,jend3_all, &
                kend_all,kend1_all,kend2_all,kend3_all, &
                ijksize3_all, displs

! Variables used for mapping i, j, k to ii, jj, kk to take care of
! of cyclic conditions...
      integer, allocatable,dimension(:) :: imap, jmap, kmap
      integer, allocatable,dimension(:) :: imap_c, jmap_c, kmap_c

      integer :: &
                ijksize3, ijkstart3,ijkend3, &
                istart3, iend3, jstart3, jend3, &
                kstart3, kend3, istart2, iend2, jstart2, jend2, &
                kstart2, kend2, istart1, iend1, jstart1, jend1, &
                kstart1, kend1

      integer :: istart, iend, jstart, jend, kstart, kend


! declaration for storing filebasename, e.g. mfix00000.dat
      CHARACTER(len=5) :: fbname
      INTEGER :: idbg = 1

! Funijk coefficients
      integer :: c0, c1, c2


      END MODULE compar
