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

      use compar, only: c0, c1, c2
      use compar, only: iend, jend, kend
      use compar, only: iend1, jend1, kend1
      use compar, only: iend2, jend2, kend2
      use compar, only: iend3, jend3, kend3
      use compar, only: imap, jmap, kmap
      use compar, only: istart, jstart, kstart
      use compar, only: istart1, jstart1, kstart1
      use compar, only: istart2, jstart2, kstart2
      use compar, only: istart3, jstart3, kstart3
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

      istart1   = domlo(1)
      iend1     = domhi(1)
      jstart1   = domlo(2)
      jend1     = domhi(2)
      kstart1   = domlo(3)
      kend1     = domhi(3)

      istart2   = domlo(1)-1
      iend2     = domhi(1)+1
      jstart2   = domlo(2)-1
      jend2     = domhi(2)+1
      kstart2   = domlo(3)-1
      kend2     = domhi(3)+1

      istart3   = domlo(1)-1
      iend3     = domhi(1)+1
      jstart3   = domlo(2)-1
      jend3     = domhi(2)+1
      kstart3   = domlo(3)-1
      kend3     = domhi(3)+1


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
      istart = istart1
      iend = iend1
      jstart = jstart1
      jend = jend1
      kstart = kstart1
      kend = kend1

      if(istart2.eq.domlo(1)-1) istart = istart2
      if(iend2  .eq.domhi(1)+1) iend = iend2
      if(jstart2.eq.domlo(2)-1) jstart = jstart2
      if(jend2  .eq.domhi(2)+1) jend = jend2
      if(kstart2.eq.domlo(3)-1) kstart = kstart2
      if(kend2  .eq.domhi(3)+1) kend = kend2

! Setup coefficients of FUINIJK
      c0 = 1 - istart3
      c1 = (iend3-istart3+1)
      c2 = (iend3-istart3+1)* (jend3-jstart3+1)
      c0 =  c0  - c1*jstart3 - c2*kstart3


! Call to sendrecv_init to set all the communication pattern

      RETURN

 1000 FORMAT('Parallel load balancing statistics:',2/,13x,'Comp. cells',&
      4X,'Processor',/3X,'maximum   ',I11,4X,I9,/3X,'minimum   ',I11,&
      4X,I9,/3X,'average   ',I11,6X,'-N/A-',2/,3X,'Maximum speedup ',&
      '(Amdahls Law) = ',A)

      END SUBROUTINE GRIDMAP_INIT

      END MODULE GRIDMAP
