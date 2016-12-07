MODULE functions

  USE geometry, only: imax1, jmax1, kmax1

! Functions for generating IJK indices for indicated basis:
!---------------------------------------------------------------------//
! INTEGER :: FUNIJK


! Logical functions to determine whether index is on my PE's domain or
! indicated subset:
!---------------------------------------------------------------------//
! LOGICAL :: IS_ON_myPE_plus1layer
! LOGICAL :: IS_ON_myPE_plus2layers


! ieast, iwest, jsouth, jnorth, kbot, ktop:
! Functions for calculating indicated directional shift in given IJK
! index. This will return the ijk index of the computational cell
! corresponding to the indicated shift when that computational cell
! is NOT a wall cell. If the computational cell is a wall cell then
! this will return its own ijk index. For example, east_of will return
! IPJK when IPJK is a fluid or flow at cell. However, if IPJK is a
! wall cell east_of will return IJK.
!---------------------------------------------------------------------//

! iminus, iplus, jminus, jplus, kminus, kplus:
! Functions for calculating indicated directional shift in given IJK
! index. This will generally return the ijk index of the computational
! cell corresponding to the indicated shift regardless of the wall
! status of that computational cell. It may not return corner cells
! unless the ijk cell itself is a corner cell.
!---------------------------------------------------------------------//

! Logical functions to identify indicated condition:
!---------------------------------------------------------------------//
! logical function to identify various fluid/flow cells
! LOGICAL :: fluid_at
! LOGICAL :: FLOW_AT

! logical function to identify various wall cells
! LOGICAL :: wall_at

! logical function to identify different flow at boundaries
! LOGICAL :: FLOW_AT_E, FLOW_AT_N, FLOW_AT_T
! LOGICAL :: MFLOW_AT_E, MFLOW_AT_N, MFLOW_AT_T

! Logical functions to identify different impermeable/semipermeable
! surface boundaries (specific type of internal surface)
! LOGICAL :: ip_at_e, ip_at_n, ip_at_t


! Additional functions
!---------------------------------------------------------------------//
! DOUBLE PRECISION :: ZMAX
! INTEGER FUNCTION :: FUNLM


CONTAINS

  INCLUDE 'functions.inc'

END MODULE functions
