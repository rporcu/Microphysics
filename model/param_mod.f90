      MODULE param

! Parameters describing problem size: (set from user input)
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Number of computational cells
      INTEGER :: DIMENSION_3  !

! Maximum number of species.
      INTEGER :: DIMENSION_N_g ! Gas
      INTEGER :: DIMENSION_N_s ! Solids


! Parameters limiting user-specifed input.
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Maximum number of reactions defined in data file
      INTEGER, PARAMETER :: DIMENSION_RXN = 100
! Number of user defined constants
      INTEGER, PARAMETER :: DIMENSION_C = 500
! Maximum number of items for specifying initial conditions
      INTEGER, PARAMETER :: DIMENSION_IC = 500
! Maximum number of items for specifying boundary conditions
      INTEGER, PARAMETER :: DIMENSION_BC = 500
! Maximum number of items for specifying point sources
      INTEGER, PARAMETER :: DIMENSION_PS = 5000
! Maximum number of solids phases
      INTEGER, PARAMETER :: DIM_M = 10
! Maximum number of gas species
      INTEGER, PARAMETER :: DIM_N_g = 100
! Maximum number of solids species per phase.
      INTEGER, PARAMETER :: DIM_N_s = 100
! Maximum of DIM_N_g and DIM_N_s
      INTEGER, PARAMETER :: DIM_N = max(DIM_N_g, DIM_N_s)
! Maximum number of species.
      INTEGER, PARAMETER :: DIM_N_ALL = 2*DIM_N
! Maximum of the number of cells in the x direction.
      INTEGER, PARAMETER :: DIM_I = 5000
! Maximum of the number of cells in the y direction.
      INTEGER, PARAMETER :: DIM_J = 5000
! Maximum of the number of cells in the z direction.
      INTEGER, PARAMETER :: DIM_K = 5000
! Maximum number of user-defined output files
      INTEGER, PARAMETER :: DIMENSION_USR = 5
! Number of Equation types:
!  1) Gas pressure
!  2) Solids volume fraction
!  3) Gas and solids U-Momentum equation
!  4) Gas and solids V-Momentum equation
!  5) Gas and solids W-Momentum equation
!  6) Temperature
!  7) Species Mass Fractions
!  8) Granular Temperature
! 10) DES Diffusion
      INTEGER, PARAMETER :: DIM_EQS = 10

      END MODULE param
