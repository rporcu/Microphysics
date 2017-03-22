module functions

! iminus, iplus, jminus, jplus, kminus, kplus:
! Functions for calculating indicated directional shift in given IJK
! index. This will generally return the ijk index of the computational
! cell corresponding to the indicated shift regardless of the wall
! status of that computational cell. It may not return corner cells
! unless the ijk cell itself is a corner cell.
!---------------------------------------------------------------------//
! Additional functions
!---------------------------------------------------------------------//
contains

  include 'functions.inc'

end module functions
