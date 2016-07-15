!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: fldvar                                                      C
!  Purpose: Common block containing field variables data               C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE fldvar

! Specified constant gas density
      DOUBLE PRECISION :: RO_g0
! Specified constant gas viscosity
      DOUBLE PRECISION :: MU_g0
! Average molecular weight of gas
      DOUBLE PRECISION :: MW_AVG

! Void fraction (current and previous time-step values)
      DOUBLE PRECISION, ALLOCATABLE ::  EP_g(:), EP_go(:)
! Gas pressure (current and previous time-step values)
      DOUBLE PRECISION, ALLOCATABLE ::  P_g(:), P_gO(:)
! Gas density (current and previous time-step values)
      DOUBLE PRECISION, ALLOCATABLE ::  RO_g(:), RO_go(:)
! Macroscopic gas density (current and previous time-step values)
      DOUBLE PRECISION, ALLOCATABLE ::  ROP_g(:), ROP_go(:)
! x-component of gas velocity (current and previous time-step values)
      DOUBLE PRECISION, ALLOCATABLE ::  U_g(:), U_go(:)
! y-component of gas velocity (current and previous time-step values)
      DOUBLE PRECISION, ALLOCATABLE ::  V_g(:), V_go(:)
! z-component of gas velocity (current and previous time-step values)
      DOUBLE PRECISION, ALLOCATABLE ::  W_g(:), W_go(:)

! Gas viscosity
      DOUBLE PRECISION, ALLOCATABLE ::  MU_g(:)
! Second coefficient of viscosity
      DOUBLE PRECISION, ALLOCATABLE ::  LAMBDA_G(:)
! trace of D_g at i, j, k
      DOUBLE PRECISION, ALLOCATABLE ::  trD_g(:)

! diffusive component of conv-dif
! stores diffusive term in matrix for x, y, z momentum cell
      DOUBLE PRECISION, ALLOCATABLE :: DF_gu(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: DF_gv(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: DF_gw(:,:)

! cross terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_U_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_W_g


!--- Physical Properties ---------------------------------------------//


!--- Pressure correction ---------------------------------------------//

      DOUBLE PRECISION, ALLOCATABLE ::  d_e(:)
      DOUBLE PRECISION, ALLOCATABLE ::  d_n(:)
      DOUBLE PRECISION, ALLOCATABLE ::  d_t(:)
!
      DOUBLE PRECISION, ALLOCATABLE ::  Pp_g(:)

      END MODULE fldvar
