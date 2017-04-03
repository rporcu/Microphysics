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
      integer :: PE_IO = 0

      end module compar
