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
   subroutine gridmap_init

      use compar, only: iend, jend, kend
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

! End setup mapping to take care of cyclic boundary conditions
! ----------------------------------------------------------------<<<

      ! Defining new set of varaibles to define upper and lower bound of the
      ! indices to include actual physical boundaries of the problem
      istart = domlo(1)-1
      iend   = domhi(1)-1
      jstart = domlo(2)-1
      jend   = domhi(2)+1
      kstart = domlo(3)+1
      kend   = domhi(3)+1

      end subroutine gridmap_init

      end module gridmap
