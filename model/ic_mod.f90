!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ic                                                          !
!  Author: M. Syamlal                                 Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: Global initial conditions variables.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE ic

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

! Maximum number of IC regions.
      use param, only: DIMENSION_IC
! Maximum number of solids phases.
      use param, only: DIM_M
! Maximum number of gas phase species
      use param, only: DIM_N_g
! Maximum number of solids phase species
      use param, only: DIM_N_s

! IC region West face, X-coordinate
      real(c_real) :: IC_X_w (DIMENSION_IC)

! IC region East face, X-coordinate
      real(c_real) :: IC_X_e (DIMENSION_IC)

! IC region South face, Y-coordinate
      real(c_real) :: IC_Y_s (DIMENSION_IC)

! IC region North face, Y-coordinate
      real(c_real) :: IC_Y_n (DIMENSION_IC)

! IC region Bottom face, Z-coordinate
      real(c_real) :: IC_Z_b (DIMENSION_IC)

! IC region Top face, Z-coordinate
      real(c_real) :: IC_Z_t (DIMENSION_IC)

! IC region, West face, I Index
      INTEGER :: IC_I_w (DIMENSION_IC)

! IC region, East face, I Index
      INTEGER :: IC_I_e (DIMENSION_IC)

! IC region, South face, J Index
      INTEGER :: IC_J_s (DIMENSION_IC)

! IC region, North face, J Index
      INTEGER :: IC_J_n (DIMENSION_IC)

! IC region, Bottom face, K Index
      INTEGER :: IC_K_b (DIMENSION_IC)

! IC region, Top face, K Index
      INTEGER :: IC_K_t (DIMENSION_IC)

! Type of initial condition: PATCH
      CHARACTER(LEN=16) :: IC_TYPE(DIMENSION_IC)

! Initial gas phase volume fraction
      real(c_real) :: IC_EP_g (DIMENSION_IC)

! Initial gas pressure
      real(c_real) :: IC_P_g (DIMENSION_IC)

! Initial macroscopic density of solids phases
      real(c_real) :: IC_ROP_s(DIMENSION_IC, DIM_M)

! Initial solids phase volume fraction
      real(c_real) :: IC_EP_s (DIMENSION_IC, DIM_M)

! Initial gas phase temperature
      real(c_real) :: IC_T_g(DIMENSION_IC)

! Initial solids phase temperature
      real(c_real) :: IC_T_s(DIMENSION_IC, DIM_M)

! Initial x-component of gas velocity
      real(c_real) :: IC_U_g(DIMENSION_IC)

! Initial x-component of solids phase velocity
      real(c_real) :: IC_U_s(DIMENSION_IC, DIM_M)

! Initial y-component of gas velocity
      real(c_real) :: IC_V_g(DIMENSION_IC)

! Initial y-component of solids phase velocity
      real(c_real) :: IC_V_s(DIMENSION_IC, DIM_M)

! Initial z-component of gas velocity
      real(c_real) :: IC_W_g(DIMENSION_IC)

! Initial z-component of solids phase velocity
      real(c_real) :: IC_W_s(DIMENSION_IC, DIM_M)

! Logical variable to determine whether an ic is defined
      LOGICAL :: IC_DEFINED (DIMENSION_IC)

! Initial gas species mass fractions
      real(c_real) :: IC_X_g(DIMENSION_IC, DIM_N_g)

! Initial solids species mass fractions
      real(c_real) :: IC_X_s(DIMENSION_IC, DIM_M, DIM_N_s)

! Flag to extend the lattice distribution in a given IC to available area
      LOGICAL :: IC_DES_FIT_TO_REGION (DIMENSION_IC)

!  Cell flag definitions
!
!     BC_TYPE        FLAG       Cell type
!  ------------------------------------------------------
!     UNDEF            0    Undefined
!     -                1    Cell containing gas or solids or both
!   P_INFLOW          10    Specified pressure inflow cell
!   P_OUTFLOW         11    Specified pressure outflow cell
!   MASS_INFLOW       20    Specified mass flux inflow cell
!   MASS_OUTFLOW      21    Specified mass flux outflow cell
!   OUTFLOW           31    outflow cell
!   NO_SLIP_WALL     100    Internal/external wall with no-slip b.c.
!   FREE_SLIP_WALL   101    Internal/external wall with free-slip
!   PAR_SLIP_WALL    102    Internal/external wall with partial-slip b.c.
!   CYCLIC           106    Cyclic b.c.
!   CYCLIC_PD        107    Cyclic b.c. with pressure drop

   integer, parameter :: UNDEF_CELL =   0
   integer, parameter :: FLUID_     =   1
   integer, parameter :: PINF_      =  10
   integer, parameter :: POUT_      =  11
   integer, parameter :: MINF_      =  20
   integer, parameter :: MOUT_      =  21
   integer, parameter :: OUTF_      =  31
   integer, parameter :: NSW_       = 100
   integer, parameter :: FSW_       = 101
   integer, parameter :: PSW_       = 102
   integer, parameter :: CYCL_      = 106
   integer, parameter :: CYCP_      = 107

      END MODULE ic
