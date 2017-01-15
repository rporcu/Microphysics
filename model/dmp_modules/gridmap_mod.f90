!--------------------------------------------------------------------
! Purpose:
! Contains following subroutines:
!    partition, gridmap_init
!--------------------------------------------------------------------
module gridmap

   implicit none

contains

!----------------------------------------------------------------------!
! Purpose: Initializing all the variables from the information         !
! obtained in the above routine.                                       !
!----------------------------------------------------------------------!
   SUBROUTINE GRIDMAP_INIT

      use compar, only: iend, jend, kend
      use compar, only: imap, jmap, kmap
      use compar, only: istart, jstart, kstart
      use compar, only: mype, numpes
      use geometry, only: cyclic_x, cyclic_y, cyclic_z
      use geometry, only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd

      use geometry, only: domlo, domhi
      use compar, only: nodesi, nodesj, nodesk

      implicit none

! Local variables
!---------------------------------------------------------------------//
! Loop indicies
      integer :: iproc, ii, jj, kk
! Local flags for cyclic boundarys
      LOGICAL :: CYC_XL, CYC_YL, CYC_ZL
! Amount of load imbalance
      INTEGER :: IMBALANCE
! Theoritical speedup (based on load imbalance)
      CHARACTER(len=32) :: AMDAHL_SPEEDUP

!......................................................................!

! Set local flags for cyclic boundaries.
      CYC_XL = (CYCLIC_X .OR. CYCLIC_X_PD)
      CYC_YL = (CYCLIC_Y .OR. CYCLIC_Y_PD)
      CYC_ZL = (CYCLIC_Z .OR. CYCLIC_Z_PD)

! Setup mapping to take care of cyclic boundary conditions
! ---------------------------------------------------------------->>>
! consider cyclic boundary condition using the imap(:),jmap(:),kmap(:)
! indirection arrays

      allocate( imap( domlo(1)-1:domhi(1)+1 ) )
      allocate( jmap( domlo(2)-1:domhi(2)+1 ) )
      allocate( kmap( domlo(3)-1:domhi(3)+1 ) )

      do kk=domlo(3)-1,domhi(3)+1
        kmap(kk) = kk
      enddo

      do jj=domlo(2)-1,domhi(2)+1
        jmap(jj) = jj
      enddo

      do ii=domlo(1)-1,domhi(1)+1
        imap(ii) = ii
      enddo

      if (CYC_ZL) then
         kmap( domhi(3)+1 ) = domlo(3)
         kmap( domlo(3)-1 ) = domhi(3)
      endif

      if (CYC_YL) then
         jmap( domhi(2)+1 ) = domlo(2)
         jmap( domlo(2)-1 ) = domhi(2)
      endif

      if (CYC_XL) then
         imap( domhi(1)+1 ) = domlo(1)
         imap( domlo(1)-1 ) = domhi(1)
      endif

! End setup mapping to take care of cyclic boundary conditions
! ----------------------------------------------------------------<<<

! Defining new set of varaibles to define upper and lower bound of the
! indices to include actual physical boundaries of the problem
      istart = domlo(1)
      iend   = domhi(1)
      jstart = domlo(2)
      jend   = domhi(2)
      kstart = domlo(3)
      kend   = domhi(3)

      istart = domlo(1)-1
      iend   = domhi(1)-1
      jstart = domlo(2)-1
      jend   = domhi(2)+1
      kstart = domlo(3)+1
      kend   = domhi(3)+1

 1000 FORMAT('Parallel load balancing statistics:',2/,13x,'Comp. cells',&
      4X,'Processor',/3X,'maximum   ',I11,4X,I9,/3X,'minimum   ',I11,&
      4X,I9,/3X,'average   ',I11,6X,'-N/A-',2/,3X,'Maximum speedup ',&
      '(Amdahls Law) = ',A)

      END SUBROUTINE GRIDMAP_INIT

      END MODULE GRIDMAP
