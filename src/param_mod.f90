      MODULE param

! Parameters describing problem size: (set from user input)
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

! Maximum number of species.
      integer :: DIMENSION_N_g ! Gas
      integer :: DIMENSION_N_s ! Solids


! Parameters limiting user-specifed input.
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Maximum number of reactions defined in data file
      integer, PARAMETER :: DIMENSION_RXN = 100
! Number of user defined constants
      integer, PARAMETER :: DIMENSION_C = 500
! Maximum number of items for specifying initial conditions
      integer, PARAMETER :: DIMENSION_IC = 500
! Maximum number of items for specifying boundary conditions
      integer, PARAMETER :: DIM_BC = 500
! Maximum number of items for specifying point sources
      integer, PARAMETER :: DIMENSION_PS = 5000
! Maximum number of solids phases
      integer, PARAMETER :: DIM_M = 10
! Maximum number of gas species
      integer, PARAMETER :: DIM_N_g = 100
! Maximum number of solids species per phase.
      integer, PARAMETER :: DIM_N_s = 100
! Maximum of DIM_N_g and DIM_N_s
      integer, PARAMETER :: DIM_N = max(DIM_N_g, DIM_N_s)
! Maximum number of species.
      integer, PARAMETER :: DIM_N_ALL = 2*DIM_N
! Maximum of the number of cells in the x direction.
      integer, PARAMETER :: DIM_I = 5000
! Maximum of the number of cells in the y direction.
      integer, PARAMETER :: DIM_J = 5000
! Maximum of the number of cells in the z direction.
      integer, PARAMETER :: DIM_K = 5000
! Maximum number of user-defined output files
      integer, PARAMETER :: DIMENSION_USR = 5
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
      integer, PARAMETER :: DIM_EQS = 10

      END MODULE param
